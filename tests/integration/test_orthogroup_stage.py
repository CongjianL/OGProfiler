from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from ogprofiler.cli import main
from ogprofiler.config import load_config
from ogprofiler.core.checkpoint import CheckpointStore
from ogprofiler.exceptions import HierarchyError, InputError
from ogprofiler.orthogroups.stage import SCHEMAS, OrthogroupConfig, run_orthogroup_stage
from ogprofiler.ux import pipeline_commands, run_status


def write_table(path, rows, schema=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows, schema=schema), path)


def make_run(root):
    write_table(
        root / "input/proteins.parquet",
        [
            dict(protein_id=0, species_id=0, original_id="same"),
            dict(protein_id=1, species_id=1, original_id="same"),
            dict(protein_id=2, species_id=0, original_id="isolate"),
        ],
    )
    write_table(root / "input/species.parquet", [dict(species_id=i) for i in range(3)])
    write_table(
        root / "components/index.parquet",
        [
            dict(protein_id=0, component_id=0),
            dict(protein_id=1, component_id=0),
            dict(protein_id=2, component_id=1),
        ],
    )
    write_table(
        root / "components/singleton_terminal_families.parquet",
        [
            dict(protein_id=2, component_id=1, terminal_reason="SINGLETON"),
        ],
    )
    hierarchy = root / "hierarchy/components/component=00000000"
    write_table(
        hierarchy / "nodes.parquet",
        [
            dict(component_id=0, cluster_id=0, parent_id=None, depth=0, n_genes=2, n_species=2),
            dict(component_id=0, cluster_id=1, parent_id=0, depth=1, n_genes=1, n_species=1),
            dict(component_id=0, cluster_id=2, parent_id=0, depth=1, n_genes=1, n_species=1),
        ],
    )
    write_table(
        hierarchy / "members.parquet",
        [
            dict(protein_id=0, terminal_cluster_id=1),
            dict(protein_id=1, terminal_cluster_id=2),
        ],
    )
    write_table(
        root / "evolution/components/component=00000000/events.parquet",
        [
            dict(component_id=0, cluster_id=0, network_event="SPECIATION"),
        ],
    )
    return root


def read_component(root, component=0):
    directory = root / "orthogroups/components" / f"component={component:08d}"
    return {name: pq.ParquetFile(directory / name).read().to_pylist() for name in SCHEMAS}


def test_roundtrip_and_manifest_verified_resume(tmp_path):
    root = make_run(tmp_path)
    config = OrthogroupConfig()
    path, components, reused = run_orthogroup_stage(root, config, ["test"])
    assert (components, reused) == (2, 0)
    rows = read_component(root)
    assert [(r["local_group_id"], r["protein_id"]) for r in rows["members.parquet"]] == [
        (0, 0),
        (0, 1),
    ]
    group = rows["groups.parquet"][0]
    assert (
        group["source_cluster_id"],
        group["selection_type"],
        group["processing_level"],
        group["n_genes"],
        group["n_species"],
    ) == (0, "EVENT_I", 2, 2, 2)
    assert rows["unassigned.parquet"] == []
    assert rows["v1_events.parquet"][0]["v1_event"] == "I"
    isolate = read_component(root, 1)
    assert isolate["groups.parquet"][0]["selection_type"] == "SSN_ISOLATE"
    assert isolate["members.parquet"][0]["protein_id"] == 2
    assert not (root / "hierarchy/components/component=00000001").exists()
    before = {
        p: p.read_bytes() for p in (root / "orthogroups/components").rglob("*") if p.is_file()
    }
    assert run_orthogroup_stage(root, config, ["resume"])[2] == 2
    assert all(p.read_bytes() == data for p, data in before.items())
    assert json.loads(path.read_text())["status"] == "DONE"


@pytest.mark.parametrize("name", list(SCHEMAS) + ["og-manifest.json"])
def test_corrupted_artifact_rebuilds_only_component(tmp_path, name):
    root = make_run(tmp_path)
    run_orthogroup_stage(root, OrthogroupConfig(), [])
    expected = read_component(root)
    path = root / "orthogroups/components/component=00000000" / name
    path.write_bytes(b"corrupt")
    assert run_orthogroup_stage(root, OrthogroupConfig(), [])[2] == 1
    assert read_component(root) == expected


@pytest.mark.parametrize(
    "change", ["overlap", "protein", "species", "events", "members", "isolate"]
)
def test_configuration_and_input_changes_invalidate_cache(tmp_path, change):
    root = make_run(tmp_path)
    run_orthogroup_stage(root, OrthogroupConfig(), [])
    config = OrthogroupConfig()
    if change == "overlap":
        config = OrthogroupConfig(species_overlap_count=1)
    else:
        relative = {
            "protein": "input/proteins.parquet",
            "species": "input/species.parquet",
            "events": "evolution/components/component=00000000/events.parquet",
            "members": "hierarchy/components/component=00000000/members.parquet",
            "isolate": "components/singleton_terminal_families.parquet",
        }[change]
        path = root / relative
        table = pq.ParquetFile(path).read()
        # Equivalent content, changed physical artifact identity.
        pq.write_table(table, path, compression="gzip")
    reused = run_orthogroup_stage(root, config, [])[2]
    assert reused == (1 if change in {"events", "members"} else 0)


def test_missing_input_partial_failure_and_recovery(tmp_path):
    root = make_run(tmp_path)
    events = root / "evolution/components/component=00000000/events.parquet"
    data = events.read_bytes()
    events.unlink()
    with pytest.raises(HierarchyError, match="components \\[0\\]"):
        run_orthogroup_stage(root, OrthogroupConfig(), [], retries=0)
    manifest = root / "orthogroups/og-manifest.json"
    assert json.loads(manifest.read_text())["status"] == "FAILED"
    store = CheckpointStore(root / "run.db")
    assert store.get("orthogroups", "0").status == "FAILED"
    assert store.get("orthogroups", "1").status == "DONE"
    assert not (root / "orthogroups/components/component=00000000/og-manifest.json").exists()
    assert (
        next(row for row in run_status(root)["stages"] if row["stage"] == "orthogroups")["status"]
        == "FAILED"
    )
    events.write_bytes(data)
    assert run_orthogroup_stage(root, OrthogroupConfig(), [])[2] == 1
    assert store.get("orthogroups", "0").status == "DONE"


def test_interrupted_publication_is_not_reused_and_retries(tmp_path, monkeypatch):
    import ogprofiler.orthogroups.stage as stage

    root = make_run(tmp_path)
    original = stage._write_rows
    calls = []

    def fail_once(path, schema, rows):
        calls.append(path.name)
        if len(calls) == 2:
            raise OSError("injected disk write failure")
        return original(path, schema, rows)

    monkeypatch.setattr(stage, "_write_rows", fail_once)
    with pytest.raises(HierarchyError):
        run_orthogroup_stage(root, OrthogroupConfig(), [], retries=0)
    assert not (root / "orthogroups/components/component=00000000/og-manifest.json").exists()
    assert not list(root.rglob(".og-build-*"))
    assert run_orthogroup_stage(root, OrthogroupConfig(), [], retries=1)[2] == 1
    assert CheckpointStore(root / "run.db").get("orthogroups", "0").attempts == 2


def test_serial_parallel_and_stale_checkpoint_resume(tmp_path):
    serial = make_run(tmp_path / "serial")
    parallel = make_run(tmp_path / "parallel")
    run_orthogroup_stage(serial, OrthogroupConfig(), [], workers=1)
    run_orthogroup_stage(parallel, OrthogroupConfig(), [], workers=2)
    for component in (0, 1):
        assert read_component(serial, component) == read_component(parallel, component)
    store = CheckpointStore(parallel / "run.db")
    store.invalidate("orthogroups", "0", "simulate interruption")
    store.start("orthogroups", "0")
    assert run_orthogroup_stage(parallel, OrthogroupConfig(), [], workers=2)[2] == 2
    assert store.get("orthogroups", "0").status == "DONE"


@pytest.mark.parametrize(
    "override",
    [
        "orthogroups.species_overlap_count=-1",
        "orthogroups.species_overlap_count=0.1",
        "orthogroups.species_overlap_count=true",
        "orthogroups.strategy=terminal",
        "orthogroups.refinement=true",
        "orthogroups.refinement=1",
    ],
)
def test_invalid_orthogroup_configuration(override):
    with pytest.raises(InputError):
        load_config(overrides=[override])


def test_cli_and_pipeline_stage_range(tmp_path):
    root = make_run(tmp_path)
    assert load_config()["orthogroups"] == {
        "strategy": "v1_compatible",
        "species_overlap_count": 0,
        "refinement": False,
    }
    assert main(["orthogroups", "--run", str(root)]) == 0
    assert main(["orthogroups", "--run", str(root), "--set", "orthogroups.refinement=true"]) == 2
    commands = pipeline_commands(
        root,
        None,
        config=None,
        overrides=(),
        from_stage="annotate-network",
        until_stage="orthogroups",
    )
    assert [c[0] for c in commands] == ["annotate-network", "orthogroups"]
    assert (
        main(
            [
                "run",
                "--out",
                str(root),
                "--from-stage",
                "orthogroups",
                "--until-stage",
                "orthogroups",
            ]
        )
        == 0
    )


def test_partial_atomic_replacement_recovers_without_completion_marker(tmp_path, monkeypatch):
    import ogprofiler.orthogroups.stage as stage

    root = make_run(tmp_path)
    run_orthogroup_stage(root, OrthogroupConfig(), [])
    expected = read_component(root)
    path = root / "orthogroups/components/component=00000000/groups.parquet"
    path.write_bytes(b"invalidates cache")
    replace = stage.os.replace

    def fail_publication(source, target):
        if Path(target).name == "members.parquet" and "component=00000000" in str(target):
            raise OSError("injected partial publication")
        return replace(source, target)

    monkeypatch.setattr(stage.os, "replace", fail_publication)
    with pytest.raises(HierarchyError):
        run_orthogroup_stage(root, OrthogroupConfig(), [], retries=0)
    directory = path.parent
    assert not (directory / "og-manifest.json").exists()
    assert json.loads((directory / "og-failure.json").read_text())["status"] == "FAILED"
    monkeypatch.setattr(stage.os, "replace", replace)
    assert run_orthogroup_stage(root, OrthogroupConfig(), [])[2] == 1
    assert read_component(root) == expected
    assert not (directory / "og-failure.json").exists()


def test_transient_failure_retries_in_same_run(tmp_path, monkeypatch):
    import ogprofiler.orthogroups.stage as stage

    root = make_run(tmp_path)
    write = stage._write_rows
    failed = False

    def fail_once(path, schema, rows):
        nonlocal failed
        if not failed:
            failed = True
            raise OSError("transient")
        return write(path, schema, rows)

    monkeypatch.setattr(stage, "_write_rows", fail_once)
    run_orthogroup_stage(root, OrthogroupConfig(), [], retries=1)
    store = CheckpointStore(root / "run.db")
    assert store.get("orthogroups", "0").attempts == 2
    assert store.get("orthogroups", "0").status == "DONE"


@pytest.mark.parametrize(
    "change", ["protein_partition", "unknown_species", "original_identity", "orphan_isolate"]
)
def test_invalid_upstream_metadata_is_rejected(tmp_path, change):
    root = make_run(tmp_path)
    if change == "orphan_isolate":
        write_table(
            root / "components/singleton_terminal_families.parquet",
            [
                dict(protein_id=99, component_id=99),
            ],
        )
    else:
        path = root / "input/proteins.parquet"
        rows = pq.ParquetFile(path).read().to_pylist()
        if change == "protein_partition":
            rows[0]["protein_id"] = 99
        elif change == "unknown_species":
            rows[0]["species_id"] = 99
        else:
            rows[2]["original_id"] = "same"
        write_table(path, rows)
    with pytest.raises(HierarchyError):
        run_orthogroup_stage(root, OrthogroupConfig(), [])
    assert not (root / "orthogroups/og-manifest.json").exists()


def test_corrupted_input_returns_cli_error(tmp_path):
    root = make_run(tmp_path)
    (root / "input/proteins.parquet").write_bytes(b"corrupt")
    assert main(["orthogroups", "--run", str(root)]) == 2


def test_stale_running_without_complete_artifact_rebuilds(tmp_path):
    root = make_run(tmp_path)
    run_orthogroup_stage(root, OrthogroupConfig(), [])
    store = CheckpointStore(root / "run.db")
    store.invalidate("orthogroups", "0", "simulate crash")
    store.start("orthogroups", "0")
    (root / "orthogroups/components/component=00000000/og-manifest.json").unlink()
    assert run_orthogroup_stage(root, OrthogroupConfig(), [])[2] == 1
    assert store.get("orthogroups", "0").status == "DONE"
    assert store.get("orthogroups", "0").attempts == 3
