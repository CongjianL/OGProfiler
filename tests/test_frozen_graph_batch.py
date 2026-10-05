import json

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from benchmarks.og_extraction.frozen_graph_batch import run_batch
from ogprofiler.core.manifest import sha256_file


@pytest.fixture
def frozen(tmp_path):
    run, ref = tmp_path / "run", tmp_path / "ref"
    ref.mkdir()
    h = run / "hierarchy/components/component=00000000"
    e = run / "components/edges/component=00000000/part.parquet"
    h.mkdir(parents=True)

    def write(path, rows):
        path.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.Table.from_pylist(rows), path)

    write(
        run / "input/proteins.parquet",
        [dict(protein_id=i, species_id=i, original_id=f"p{i}") for i in range(3)],
    )
    write(
        run / "input/species.parquet", [dict(species_id=i, species_name=f"s{i}") for i in range(3)]
    )
    write(
        run / "components/index.parquet",
        [dict(protein_id=i, component_id=0 if i < 2 else 1) for i in range(3)],
    )
    (run / "run.yaml").write_text("{}")
    write(
        h / "nodes.parquet",
        [
            dict(cluster_id=0, parent_id=None, n_genes=2),
            dict(cluster_id=1, parent_id=0, n_genes=1),
            dict(cluster_id=2, parent_id=0, n_genes=1),
        ],
    )
    write(
        h / "members.parquet",
        [dict(protein_id=0, terminal_cluster_id=1), dict(protein_id=1, terminal_cluster_id=2)],
    )
    write(e, [dict(u=0, v=1, weight=2.0)])
    (h / "hierarchy-manifest.json").write_text(
        json.dumps(
            dict(
                hierarchy_status="RESOLVED",
                structural_validation_passed=True,
                input_checksums={str(e.relative_to(run)): sha256_file(e)},
                output_checksums={
                    name: sha256_file(h / name) for name in ["nodes.parquet", "members.parquet"]
                },
            )
        )
    )
    (ref / "Orthogroups.tsv").write_text("Orthogroup\ts0\ts1\ts2\nOG0\tp0\tp1\t\n")
    (ref / "Orthogroups_UnassignedGenes.tsv").write_text("Orthogroup\ts0\ts1\ts2\nU0\t\t\tp2\n")
    return run, ref, e


def test_frozen_batch_predictions_then_evaluation(frozen, tmp_path):
    run, ref, _ = frozen
    out = tmp_path / "out"
    r = run_batch(run, ref, out, [1.0, 3.0])
    assert r["completed"]
    assert r["candidates"][0]["primary"]["pair_f1"] == 1
    assert r["candidates"][1]["primary"]["pair_recall"] == 0
    assert len((out / "candidate-00/members.tsv").read_text().splitlines()) == 4
    assert r["candidates"][0]["strata"]["cross_species"]["recall"] == 1
    # Reference perturbation affects evaluation, never predictions.
    (ref / "Orthogroups.tsv").write_text("Orthogroup\ts0\ts1\ts2\nOG0\tp0\t\t\nOG1\t\tp1\t\n")
    out2 = tmp_path / "out2"
    run_batch(run, ref, out2, [1.0, 3.0])
    for i in range(2):
        assert (out / f"candidate-{i:02d}/members.tsv").read_bytes() == (
            out2 / f"candidate-{i:02d}/members.tsv"
        ).read_bytes()
    with pytest.raises(FileExistsError):
        run_batch(run, ref, out, [1.0])


def test_corrupt_or_extra_edges_rejected(frozen, tmp_path):
    run, ref, edge = frozen
    extra = edge.with_name("extra.parquet")
    extra.write_bytes(edge.read_bytes())
    with pytest.raises(ValueError, match="Edge partition"):
        run_batch(run, ref, tmp_path / "extra", [1.0])
    extra.unlink()
    edge.write_bytes(b"bad")
    with pytest.raises(ValueError, match="Checksum"):
        run_batch(run, ref, tmp_path / "corrupt", [1.0])


@pytest.mark.parametrize("penalties", [[], [1.0, 1.0], [-1.0], [float("nan")]])
def test_bad_protocol_before_output(frozen, tmp_path, penalties):
    run, ref, _ = frozen
    with pytest.raises(ValueError, match="penalties"):
        run_batch(run, ref, tmp_path / "bad", penalties)
    assert not (tmp_path / "bad").exists()


def test_strength_batch_fixed_objective_and_reference_independence(frozen, tmp_path):
    run, ref, _ = frozen
    r = run_batch(run, ref, tmp_path / "strength", [None], strength_null=True)
    assert r["protocol"]["algorithm"] == "fixed-tree-component-strength-null-v1"
    assert r["protocol"]["resolution"] == 1
    assert r["candidates"][0]["primary"]["pair_f1"] == 1
    (ref / "Orthogroups.tsv").write_text("Orthogroup\ts0\ts1\ts2\nOG0\tp0\t\t\nOG1\t\tp1\t\n")
    run_batch(run, ref, tmp_path / "perturbed", [None], strength_null=True)
    assert (tmp_path / "strength/candidate-00/members.tsv").read_bytes() == (
        tmp_path / "perturbed/candidate-00/members.tsv"
    ).read_bytes()
    with pytest.raises(ValueError, match="fixed resolution"):
        run_batch(run, ref, tmp_path / "bad-strength", [0.1], strength_null=True)
