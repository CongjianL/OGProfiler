import csv
import json

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from benchmarks.qfo.evaluate import evaluate, of_assigned, partition
from ogprofiler.core.manifest import sha256_file


def tsv(path, fields, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(fields)
        w.writerows(data)


def fixture(tmp_path):
    names = ["sp|A|ONE", "sp|B|TWO", "sp|C|THREE", "sp|D|FOUR"]
    roots = []
    for method in ("v2", "of3", "v1"):
        r = tmp_path / method
        roots.append(r)
        inp = r / "campaign/input"
        inp.mkdir(parents=True)
        for i in range(2):
            (inp / f"s{i}.fasta").write_text(
                "".join(f">{n}\nACDE\n" for n in names[i * 2 : (i + 1) * 2])
            )
        hashes = {p.name: sha256_file(p) for p in inp.glob("*.fasta")}
        (r / "campaign/campaign-input.json").write_text(
            json.dumps(
                dict(
                    collection="bacteria",
                    mode="full",
                    n_proteins=4,
                    n_species=2,
                    input_sha256=hashes,
                    dataset_sha256="fixture",
                )
            )
        )
        (r / "qfo-methods-completion.json").write_text(
            json.dumps(dict(methods_completed=True, mode="full", method=method))
        )
        tsv(
            r / f"campaign/{method}-timing/metadata.tsv",
            ["key", "value"],
            [
                ("exit_code", "0"),
                ("start_time_utc", "2026-10-07T00:00:00Z"),
                ("end_time_utc", "2026-10-07T00:01:00Z"),
            ],
        )
    v2, of3, v1 = roots
    p = v2 / "campaign/v2/input/proteins.parquet"
    p.parent.mkdir(parents=True)
    pq.write_table(
        pa.Table.from_pylist([dict(protein_id=i, original_id=n) for i, n in enumerate(names)]), p
    )
    p = v2 / "campaign/v2/components/index.parquet"
    p.parent.mkdir(parents=True)
    pq.write_table(pa.Table.from_pylist([dict(protein_id=i, component_id=0) for i in range(4)]), p)
    tsv(
        v2 / "campaign/v2/results/members.tsv",
        ["original_id", "family_id"],
        zip(names, ["P1", "P1", "P2", "P2"], strict=True),
    )
    tsv(
        v1 / "campaign/v1-groups.tsv",
        ["protein_id", "group_id"],
        zip(names, ["P1", "P2", "P1", "P3"], strict=True),
    )
    tsv(
        v1 / "campaign/v1-id-map.tsv",
        ["adapted_id", "original_id"],
        zip([f"x{i}" for i in range(4)], names, strict=True),
    )
    tsv(
        of3 / "campaign/of3-groups.tsv",
        ["protein_id", "group_id"],
        zip(names, ["R1", "R1", "R2", "UNASSIGNED_D"], strict=True),
    )
    primary = of3 / "primary.tsv"
    tsv(
        primary,
        ["Orthogroup", "s0", "s1"],
        [("R1", ", ".join(names[:2]), ""), ("R2", "", names[2])],
    )
    (of3 / "campaign/of3-assigned-primary.txt").write_text(str(primary))
    return roots, names


def test_scores_assigned_only_and_reports_unassigned_separately(tmp_path):
    roots, _ = fixture(tmp_path)
    result = evaluate(*roots, tmp_path / "metrics")
    assert result["of3_assigned"] == 3 and result["of3_unassigned"] == 1
    assert result["metrics"]["V2"]["pair_f1"] == 1
    assert result["metrics"]["V1"]["pair_f1"] == 0
    assert result["metrics"]["V1"]["bcubed_f1"] == pytest.approx(2 / 3)
    assert result["coverage"]["V2"]["of3_unassigned_attached_to_assigned"] == 1
    assert not result["official_qfo_scores_produced"]
    assert result["campaign_method_span_seconds"] == 60


def test_rejects_changed_input_across_methods(tmp_path):
    roots, _ = fixture(tmp_path)
    p = roots[1] / "campaign/input/s0.fasta"
    p.write_text(p.read_text().replace("ACDE", "ACDF"))
    with pytest.raises(ValueError, match="changed"):
        evaluate(*roots, tmp_path / "metrics")


def test_rejects_incomplete_or_duplicate_partition(tmp_path):
    p = tmp_path / "p.tsv"
    tsv(p, ["protein_id", "group_id"], [("a", "G1"), ("a", "G2")])
    with pytest.raises(ValueError, match="Duplicate"):
        partition(p, "protein_id", "group_id", {"a", "b"})
    tsv(p, ["protein_id", "group_id"], [("a", "G1")])
    with pytest.raises(ValueError, match="Incomplete"):
        partition(p, "protein_id", "group_id", {"a", "b"})


def test_rejects_of_unknown_and_duplicate_ids(tmp_path):
    p = tmp_path / "of.tsv"
    for ids in ("a, a", "a, unknown"):
        tsv(p, ["Orthogroup", "s0"], [("R1", ids)])
        with pytest.raises(ValueError, match="Unknown/duplicate"):
            of_assigned(p, {"a"}, ["s0"])
