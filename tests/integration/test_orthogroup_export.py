from __future__ import annotations

import csv
import json
from pathlib import Path

import pyarrow.parquet as pq
import pytest
from test_orthogroup_stage import make_run, write_table

from ogprofiler.cli import main
from ogprofiler.core.manifest import sha256_file
from ogprofiler.evolution.stage import run_network_annotation_stage
from ogprofiler.exceptions import ExportError, PhylogenyError
from ogprofiler.input.fasta import parse_fasta
from ogprofiler.orthogroups.models import OrthogroupConfig
from ogprofiler.orthogroups.stage import run_orthogroup_stage
from ogprofiler.orthology.stage import run_orthology_stage
from ogprofiler.output.stage import TABLE_FIELDS, run_export_stage
from ogprofiler.phylogeny.selection import select_refinement_families


def tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def export_fixture(root, workers=1, add_unassigned=False):
    make_run(root)
    hierarchy = root / "hierarchy/components/component=00000000"
    nodes = pq.ParquetFile(hierarchy / "nodes.parquet").read().to_pylist()
    for row in nodes:
        row["child_count"] = 2 if row["cluster_id"] == 0 else 0
        row["terminal_reason"] = None if row["cluster_id"] == 0 else "MIN_SIZE"
    write_table(hierarchy / "nodes.parquet", nodes)
    if add_unassigned:
        path = root / "input/proteins.parquet"
        write_table(
            path,
            pq.ParquetFile(path).read().to_pylist()
            + [
                dict(protein_id=3, species_id=2, original_id="unselected"),
            ],
        )
        path = root / "components/index.parquet"
        write_table(
            path, pq.ParquetFile(path).read().to_pylist() + [dict(protein_id=3, component_id=2)]
        )
        directory = root / "hierarchy/components/component=00000002"
        write_table(
            directory / "nodes.parquet",
            [
                dict(
                    component_id=2,
                    cluster_id=0,
                    parent_id=None,
                    depth=0,
                    n_genes=1,
                    n_species=1,
                    child_count=0,
                    terminal_reason="MIN_SIZE",
                )
            ],
        )
        write_table(directory / "members.parquet", [dict(protein_id=3, terminal_cluster_id=0)])
    (root / "input/proteins.faa").write_text(
        ">OGP2P000000000000\nAAAA\n>OGP2P000000000001\nCCCC\n>OGP2P000000000002\nGGGG\n"
        + (">OGP2P000000000003\nTTTT\n" if add_unassigned else "")
    )
    run_network_annotation_stage(root, 0.0, [])
    run_orthogroup_stage(root, OrthogroupConfig(), [], workers=workers)
    return root


def test_default_cli_exports_internal_node_og_fasta_and_separate_terminals(tmp_path):
    root = export_fixture(tmp_path)
    assert main(["export", "--run", str(root), "--all-family-fasta"]) == 0
    results = root / "results"
    groups = tsv(results / "families.tsv")
    assert [
        (g["family_id"], g["n_genes"], g["selection_type"], g["source_cluster_id"]) for g in groups
    ] == [
        ("OG000000000", "1", "SSN_ISOLATE", ""),
        ("OG000000001", "2", "EVENT_I", "0"),
    ]
    assert {
        m["family_id"] for m in tsv(results / "members.tsv") if m["protein_id"] in {"0", "1"}
    } == {"OG000000001"}
    terminal = tsv(results / "terminal_families.tsv")
    assert [(g["family_id"], g["n_genes"]) for g in terminal] == [
        ("TF000000000", "1"),
        ("TF000000001", "1"),
        ("TF000000002", "1"),
    ]
    fasta = list(parse_fasta(results / "fasta/OG000000001.faa", "error"))
    assert [r.identifier for r in fasta] == ["OGP2P000000000000", "OGP2P000000000001"]
    assert [r.sequence for r in fasta] == ["AAAA", "CCCC"]
    assert (results / "fasta/OG000000001.faa").read_text().splitlines()[::2] == [
        ">OGP2P000000000000 original_id=same protein_id=0 species_id=0",
        ">OGP2P000000000001 original_id=same protein_id=1 species_id=1",
    ]
    stats = {row["metric"]: int(row["value"]) for row in tsv(results / "statistics.tsv")}
    assert stats == dict(
        orthogroups=2,
        terminal_families=3,
        input_proteins=3,
        assigned_proteins=3,
        unassigned_proteins=0,
        duplicate_assignments=0,
    )
    manifest = json.loads((results / "export-manifest.json").read_text())
    assert manifest["counts"]["families"] == manifest["counts"]["orthogroups"] == 2
    assert all(
        g["membership_hash"] == manifest["family_membership_sha256"][g["family_id"]] for g in groups
    )


def test_unassigned_is_preserved_and_not_forced_to_singleton(tmp_path):
    root = export_fixture(tmp_path, add_unassigned=True)
    run_export_stage(root, [], all_family_fasta=True)
    results = root / "results"
    assert len(tsv(results / "families.tsv")) == 2
    assert {m["protein_id"] for m in tsv(results / "members.tsv")} == {"0", "1", "2"}
    assert tsv(results / "unassigned.tsv") == [
        dict(
            component_id="2",
            protein_id="3",
            species_id="2",
            original_id="unselected",
            terminal_cluster_id="0",
            reason="NO_V1_SELECTION",
        )
    ]
    assert len(tsv(results / "terminal_members.tsv")) == 4
    assert "TTTT" not in "".join(p.read_text() for p in (results / "fasta").glob("*.faa"))


def test_global_ids_and_all_exchange_tables_match_serial_parallel(tmp_path):
    serial = export_fixture(tmp_path / "serial", workers=1)
    parallel = export_fixture(tmp_path / "parallel", workers=2)
    run_export_stage(serial, [])
    run_export_stage(parallel, [])
    assert all(
        (serial / "results" / name).read_bytes() == (parallel / "results" / name).read_bytes()
        for name in TABLE_FIELDS
    )


@pytest.mark.parametrize(
    "name", list(TABLE_FIELDS) + ["fasta/OG000000001.faa", "export-manifest.json"]
)
def test_export_resume_repairs_every_managed_artifact(tmp_path, name):
    root = export_fixture(tmp_path)
    _, reused, count = run_export_stage(root, [], all_family_fasta=True)
    assert (reused, count) == (False, 2)
    expected = (root / "results" / name).read_bytes()
    assert run_export_stage(root, [], all_family_fasta=True)[1]
    (root / "results" / name).write_bytes(b"corrupt")
    assert not run_export_stage(root, [], all_family_fasta=True)[1]
    if name != "export-manifest.json":
        assert (root / "results" / name).read_bytes() == expected


def test_changed_upstream_or_og_config_is_rejected_until_og_stage_rebuilt(tmp_path):
    root = export_fixture(tmp_path)
    run_export_stage(root, [])
    path = root / "hierarchy/components/component=00000000/members.parquet"
    pq.write_table(pq.ParquetFile(path).read(), path, compression="gzip")
    with pytest.raises(ExportError, match="stale or corrupt"):
        run_export_stage(root, [])
    run_orthogroup_stage(root, OrthogroupConfig(), [])
    assert not run_export_stage(root, [])[1]
    with pytest.raises(ExportError, match="configuration changed"):
        run_export_stage(root, [], config=OrthogroupConfig(species_overlap_count=1))


@pytest.mark.parametrize(
    "artifact",
    [
        "groups.parquet",
        "members.parquet",
        "unassigned.parquet",
        "og-manifest.json",
        "root_manifest",
    ],
)
def test_export_checks_og_inputs_before_using_export_cache(tmp_path, artifact):
    root = export_fixture(tmp_path)
    run_export_stage(root, [])
    path = (
        root / "orthogroups/og-manifest.json"
        if artifact == "root_manifest"
        else root / "orthogroups/components/component=00000000" / artifact
    )
    path.write_bytes(b"corrupt")
    with pytest.raises(ExportError):
        run_export_stage(root, [])


def test_explicit_terminal_diagnostic_does_not_replace_default_og_results(tmp_path):
    root = export_fixture(tmp_path)
    run_export_stage(root, [])
    before = (root / "results/families.tsv").read_bytes()
    path, reused, count = run_export_stage(root, [], strategy="terminal")
    assert path == root / "results/terminal-families/export-manifest.json"
    assert (reused, count) == (False, 3)
    assert (root / "results/families.tsv").read_bytes() == before
    assert main(["export", "--run", str(root), "--strategy", "terminal"]) == 0


def test_old_export_cache_is_replaced_and_obsolete_managed_fasta_archived(tmp_path):
    root = export_fixture(tmp_path)
    manifest, _, _ = run_export_stage(root, [], all_family_fasta=True)
    stale = root / "results/fasta/OG000000000.faa"
    original = stale.read_bytes()
    value = json.loads(manifest.read_text())
    value["algorithm_version"] = "terminal-family-export-v1"
    manifest.write_text(json.dumps(value))
    assert not run_export_stage(root, [], fasta_families=("OG000000001",))[1]
    assert not stale.exists()
    archived = list((root / "results/previous-fasta").rglob("OG000000000.faa"))
    assert len(archived) == 1 and archived[0].read_bytes() == original


def test_pairwise_policy_is_independent_of_og_grouping_and_cache(tmp_path):
    root = export_fixture(tmp_path)
    path, reused, pairs = run_orthology_stage(root, [], chunk_size=10)
    assert not reused and pairs == 1
    policy = json.loads(path.read_text())["parameters"]
    assert policy["grouping_dependency"] == "none"
    assert policy["event_source"] == "network_event"
    (root / "orthogroups/og-manifest.json").write_text("invalid")
    assert run_orthology_stage(root, [], chunk_size=10)[1]


def test_phylogeny_selects_terminal_diagnostics_not_internal_og(tmp_path):
    root = export_fixture(tmp_path)
    run_export_stage(root, [])
    selected = select_refinement_families(
        root,
        explicit_family_ids=("TF000000001",),
        selection_events=set(),
        large_family_size=100,
        max_families=1,
    )
    assert len(selected) == 1
    assert selected[0].cluster_id == 1 and selected[0].n_genes == 1
    with pytest.raises(PhylogenyError, match="Unknown family"):
        select_refinement_families(
            root,
            explicit_family_ids=("OG000000001",),
            selection_events=set(),
            large_family_size=100,
            max_families=1,
        )


def test_cli_end_to_end_stage_range_reuses_hierarchy(tmp_path):
    root = export_fixture(tmp_path)
    nodes = root / "hierarchy/components/component=00000000/nodes.parquet"
    checksum = sha256_file(nodes)
    assert (
        main(
            [
                "run",
                "--out",
                str(root),
                "--from-stage",
                "annotate-network",
                "--until-stage",
                "export",
            ]
        )
        == 0
    )
    assert sha256_file(nodes) == checksum
    assert (root / "results/statistics.tsv").is_file()


def test_protein_and_component_reindexing_preserves_global_og_ids(tmp_path):
    root = export_fixture(tmp_path)
    run_export_stage(root, [])
    before = [
        (row["family_id"], row["membership_hash"]) for row in tsv(root / "results/families.tsv")
    ]
    mapping = {0: 10, 1: 11, 2: 12}
    components = {0: 7, 1: 5}
    paths = [
        root / "input/proteins.parquet",
        root / "components/index.parquet",
        root / "components/singleton_terminal_families.parquet",
        root / "hierarchy/components/component=00000000/nodes.parquet",
        root / "hierarchy/components/component=00000000/members.parquet",
    ]
    for path in paths:
        rows = pq.ParquetFile(path).read().to_pylist()
        for row in rows:
            if "protein_id" in row:
                row["protein_id"] = mapping[row["protein_id"]]
            if "component_id" in row:
                row["component_id"] = components[row["component_id"]]
        write_table(path, list(reversed(rows)))
    (root / "hierarchy/components/component=00000000").rename(
        root / "hierarchy/components/component=00000007"
    )
    (root / "input/proteins.faa").write_text(
        ">OGP2P000000000010\nAAAA\n>OGP2P000000000011\nCCCC\n>OGP2P000000000012\nGGGG\n"
    )
    run_network_annotation_stage(root, 0, [])
    run_orthogroup_stage(root, OrthogroupConfig(), [])
    run_export_stage(root, [], all_family_fasta=True)
    assert [
        (row["family_id"], row["membership_hash"]) for row in tsv(root / "results/families.tsv")
    ] == before
    assert {row["protein_id"] for row in tsv(root / "results/members.tsv")} == {"10", "11", "12"}


def test_optional_pairwise_pipeline_is_explicit_and_uses_event_engine(tmp_path):
    root = export_fixture(tmp_path)
    assert (
        main(
            [
                "run",
                "--out",
                str(root),
                "--from-stage",
                "orthogroups",
                "--until-stage",
                "export",
                "--set",
                "output.emit_pairwise_orthologs=true",
            ]
        )
        == 0
    )
    manifest = json.loads((root / "results/ortholog-manifest.json").read_text())
    assert manifest["counts"]["pairs"] == 1
    assert manifest["parameters"]["grouping_dependency"] == "none"


def test_export_partial_write_failure_has_no_completion_marker_and_recovers(tmp_path, monkeypatch):
    import ogprofiler.output.stage as stage

    root = export_fixture(tmp_path)
    manifest, _, _ = run_export_stage(root, [])
    (root / "results/members.tsv").write_bytes(b"corrupt")
    writer = stage.write_tsv

    def fail_members(path, fields, rows):
        if path.name == "members.tsv":
            raise OSError("injected disk failure")
        return writer(path, fields, rows)

    monkeypatch.setattr(stage, "write_tsv", fail_members)
    with pytest.raises(OSError, match="disk failure"):
        run_export_stage(root, [])
    assert not manifest.exists()
    monkeypatch.setattr(stage, "write_tsv", writer)
    assert not run_export_stage(root, [])[1]
    assert manifest.exists()


def test_semantically_inconsistent_member_hash_is_not_exported(tmp_path):
    root = export_fixture(tmp_path)
    directory = root / "orthogroups/components/component=00000000"
    path = directory / "groups.parquet"
    rows = pq.ParquetFile(path).read().to_pylist()
    rows[0]["membership_hash"] = "invalid"
    write_table(path, rows, pq.ParquetFile(path).schema_arrow)
    component_manifest = directory / "og-manifest.json"
    value = json.loads(component_manifest.read_text())
    value["output_checksums"]["groups.parquet"] = sha256_file(path)
    component_manifest.write_text(json.dumps(value))
    root_manifest = root / "orthogroups/og-manifest.json"
    value = json.loads(root_manifest.read_text())
    value["component_manifests"]["0"] = sha256_file(component_manifest)
    root_manifest.write_text(json.dumps(value))
    with pytest.raises(ExportError, match="membership hash mismatch"):
        run_export_stage(root, [])


def test_downstream_hierarchy_metrics_use_terminal_diagnostics(tmp_path):
    from ogprofiler.benchmark.metrics import evolution_metrics
    from ogprofiler.benchmark.synthetic_metrics import _predicted_clades

    root = export_fixture(tmp_path)
    path = root / "input/proteins.parquet"
    rows = pq.ParquetFile(path).read().to_pylist()
    rows[0]["original_id"], rows[1]["original_id"] = "alpha", "beta"
    write_table(path, rows)
    run_orthogroup_stage(root, OrthogroupConfig(), [])
    run_export_stage(root, [])
    evolution = evolution_metrics(root)
    assert evolution["evaluated_network_nodes"] == 1
    assert evolution["species_overlap_consistency"] == 1.0
    clades, events = _predicted_clades(root)
    assert clades == {frozenset({"alpha", "beta"})}
    assert events[frozenset({"alpha", "beta"})] == "SPECIATION_LIKE"


def test_cross_species_members_within_og_do_not_imply_pairwise_orthologs(tmp_path):
    import pyarrow as pa

    root = tmp_path
    write_table(
        root / "input/proteins.parquet",
        [dict(protein_id=p, species_id=p % 2, original_id=f"g{p}") for p in range(4)],
    )
    write_table(root / "input/species.parquet", [dict(species_id=s) for s in range(2)])
    write_table(
        root / "components/index.parquet", [dict(protein_id=p, component_id=0) for p in range(4)]
    )
    write_table(
        root / "components/singleton_terminal_families.parquet",
        [],
        pa.schema([("protein_id", pa.int64()), ("component_id", pa.int64())]),
    )
    directory = root / "hierarchy/components/component=00000000"
    write_table(
        directory / "nodes.parquet",
        [
            dict(
                component_id=0,
                cluster_id=0,
                parent_id=None,
                depth=0,
                n_genes=4,
                n_species=2,
                child_count=2,
                terminal_reason=None,
            ),
            dict(
                component_id=0,
                cluster_id=1,
                parent_id=0,
                depth=1,
                n_genes=2,
                n_species=2,
                child_count=0,
                terminal_reason="MIN_SIZE",
            ),
            dict(
                component_id=0,
                cluster_id=2,
                parent_id=0,
                depth=1,
                n_genes=2,
                n_species=2,
                child_count=0,
                terminal_reason="MIN_SIZE",
            ),
        ],
    )
    write_table(
        directory / "members.parquet",
        [dict(protein_id=p, terminal_cluster_id=1 + p // 2) for p in range(4)],
    )
    (root / "input/proteins.faa").write_text("".join(f">OGP2P{p:012d}\nAAAA\n" for p in range(4)))
    run_network_annotation_stage(root, 0, [])
    run_orthogroup_stage(root, OrthogroupConfig(), [])
    run_export_stage(root, [])
    groups = tsv(root / "results/families.tsv")
    assert len(groups) == 2 and all(g["n_species"] == "2" for g in groups)
    assert run_orthology_stage(root, [], chunk_size=10)[2] == 0


def test_source_cli_subprocess_runs_parallel_og_and_default_export(tmp_path):
    import os
    import subprocess
    import sys

    root = export_fixture(tmp_path)
    repo = Path(__file__).parents[2]
    environment = {**os.environ, "PYTHONPATH": str(repo / "src")}
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "ogprofiler",
            "run",
            "--out",
            str(root),
            "--from-stage",
            "orthogroups",
            "--until-stage",
            "export",
            "--set",
            "runtime.workers=2",
            "--set",
            "orthogroups.species_overlap_count=1",
        ],
        cwd=repo,
        env=environment,
        text=True,
        capture_output=True,
        timeout=30,
    )
    assert completed.returncode == 0, completed.stderr
    manifest = json.loads((root / "orthogroups/og-manifest.json").read_text())
    assert manifest["parameters"]["species_overlap_count"] == 1
    assert manifest["runtime"]["workers"] == 2
    assert len(tsv(root / "results/families.tsv")) == 2
    assert "pipeline complete" in completed.stdout
