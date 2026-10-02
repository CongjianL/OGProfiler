"""Production unresolved/resume/publication contracts (small local fixtures)."""

from __future__ import annotations

import json
from pathlib import Path

import pyarrow.parquet as pq
import pytest
from test_hierarchy_scheduler import _prepare_multicomponent_run

from ogprofiler.cli import main
from ogprofiler.config import hierarchy_config, load_config
from ogprofiler.core.checkpoint import CheckpointStore
from ogprofiler.exceptions import HierarchyError
from ogprofiler.hierarchy.gate import require_resolved_hierarchy
from ogprofiler.hierarchy.stage import run_hierarchy_component_stage
from ogprofiler.orthogroups.models import OrthogroupConfig
from ogprofiler.orthogroups.stage import run_orthogroup_stage, verified_orthogroup_inputs


@pytest.mark.parametrize("topology", ["kway_v1", "soft_binary_24_v2"])
def test_h4_acceptance_gates_on_complete_serial_parallel_replay(tmp_path: Path, topology):
    import shutil

    import yaml

    from benchmarks.og_extraction.hierarchy_regression import validate_h4

    run = _prepare_multicomponent_run(tmp_path)
    baseline = load_config(str(run / "run.yaml"))
    baseline["hierarchy"]["subtree_release_size"] = 1
    baseline["hierarchy"]["topology_policy"] = topology
    (run / "run.yaml").write_text(yaml.safe_dump(baseline))
    parallel = tmp_path / "parallel"
    parallel.mkdir()
    for name in ("input", "edges", "components"):
        (parallel / name).symlink_to(run / name, target_is_directory=True)
    shutil.copy2(run / "run.yaml", parallel / "run.yaml")
    assert main(["hierarchy", "--run", str(run), "--component-id", "0"]) == 0
    assert (
        main(
            [
                "hierarchy",
                "--run",
                str(parallel),
                "--component-id",
                "0",
                "--set",
                "hierarchy.subtree_workers=2",
            ]
        )
        == 0
    )
    acceptance = validate_h4(run, parallel)
    assert acceptance["passed"]
    if topology == "soft_binary_24_v2":
        assert acceptance["checks"]["fallback_evidence"]
        assert acceptance["hierarchy"]["fallback_count"] > 0
    folder = parallel / "hierarchy/components/component=00000000"
    members = pq.ParquetFile(folder / "members.parquet").read()
    pq.write_table(members.slice(1), folder / "members.parquet")
    report = validate_h4(run, parallel)
    assert not report["passed"]
    assert not report["checks"]["parallel_manifest"]
    assert not report["checks"]["serial_parallel_equal"]


def test_iteration_budget_is_recorded_and_invalidates_resume(tmp_path: Path):
    from dataclasses import replace

    from ogprofiler.hierarchy.stage import hierarchy_component_is_verified

    run = _prepare_multicomponent_run(tmp_path)
    config = hierarchy_config(load_config()["hierarchy"])
    old = replace(config, leiden_iterations=2)
    output, resumed = run_hierarchy_component_stage(run, 0, old, [])
    assert not resumed
    assert hierarchy_component_is_verified(run, 0, old)
    assert not hierarchy_component_is_verified(run, 0, config)
    output, resumed = run_hierarchy_component_stage(run, 0, config, [])
    assert not resumed
    manifest_path = output / "hierarchy-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["parameters"]["hierarchy"]["leiden_iterations"] == 10
    assert manifest["algorithm_version"] == "hierarchical-leiden-v4"
    assert manifest["schema_version"] == 3
    assert run_hierarchy_component_stage(run, 0, config, [])[1]
    manifest["algorithm_version"] = "hierarchical-leiden-v2"
    manifest_path.write_text(json.dumps(manifest))
    assert not hierarchy_component_is_verified(run, 0, config)


def test_unresolved_diagnostics_resume_and_explicit_retry(tmp_path: Path):
    run = _prepare_multicomponent_run(tmp_path)
    command = ["hierarchy-all", "--run", str(run), "--set", "hierarchy.gamma_max=.1"]
    assert main(command) == 2
    manifest_path = run / "hierarchy/scheduler-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["counts"]["unresolved"] == 2
    assert manifest["counts"]["failed"] == 0
    store = CheckpointStore(run / "run.db")
    attempts = store.get("hierarchy", "0").attempts
    assert store.get("hierarchy", "0").status == "UNRESOLVED"
    output = run / "hierarchy/components/component=00000000"
    assert pq.ParquetFile(output / "members.parquet").read().num_rows == 6
    candidates = pq.ParquetFile(output / "candidates.parquet").read()
    assert candidates.num_rows <= 24
    assert "violations" in candidates.column_names
    assert main(command) == 2
    assert store.get("hierarchy", "0").attempts == attempts
    assert json.loads(manifest_path.read_text())["counts"]["scheduled"] == 0
    for consumer in (
        lambda: require_resolved_hierarchy(run),
        lambda: run_orthogroup_stage(run, OrthogroupConfig(), []),
        lambda: verified_orthogroup_inputs(run),
    ):
        with pytest.raises(HierarchyError, match="UNRESOLVED"):
            consumer()
    assert main(["export", "--run", str(run)]) == 2
    assert not (run / "results/export-manifest.json").is_file()
    assert main([*command, "--retry-unresolved"]) == 2
    assert store.get("hierarchy", "0").attempts == attempts + 1


def test_gate_uses_recorded_effective_config_not_prepare_defaults(tmp_path: Path):
    run = _prepare_multicomponent_run(tmp_path)
    config = hierarchy_config(
        load_config(overrides=["hierarchy.recursion_stop_size=10"])["hierarchy"]
    )
    for component in (0, 1):
        run_hierarchy_component_stage(run, component, config, [])
    require_resolved_hierarchy(run)


def test_checkpoint_v1_migration_preserves_state(tmp_path: Path):
    import sqlite3

    from ogprofiler.core.workspace import Workspace

    run = tmp_path / "run"
    run.mkdir()
    with sqlite3.connect(run / "run.db") as db:
        db.execute("""CREATE TABLE tasks (
            stage TEXT NOT NULL, task_id TEXT NOT NULL,
            status TEXT CHECK(status IN ('PENDING','RUNNING','DONE','FAILED','INVALID')),
            started_at TEXT, completed_at TEXT, input_hash TEXT, output_path TEXT,
            error TEXT, attempts INTEGER NOT NULL DEFAULT 0, algorithm_version TEXT,
            PRIMARY KEY(stage,task_id))""")
        db.execute(
            "INSERT INTO tasks(stage,task_id,status,attempts,algorithm_version) "
            "VALUES('hierarchy','0','RUNNING',3,'v1')"
        )
    Workspace.create(run)
    store = CheckpointStore(run / "run.db")
    assert store.get("hierarchy", "0").attempts == 3
    store.unresolved("hierarchy", "0", "diagnostics")
    Workspace.create(run)
    assert store.get("hierarchy", "0").status == "UNRESOLVED"
    assert store.get("hierarchy", "0").attempts == 3


def test_gate_invalidates_changed_baseline_config(tmp_path: Path):
    run = _prepare_multicomponent_run(tmp_path)
    config = hierarchy_config(
        load_config(overrides=["hierarchy.recursion_stop_size=10"])["hierarchy"]
    )
    for component in (0, 1):
        run_hierarchy_component_stage(run, component, config, [])
    require_resolved_hierarchy(run)
    with (run / "run.yaml").open("a") as handle:
        handle.write("\n# changed config provenance\n")
    with pytest.raises(HierarchyError, match="input checksum mismatch"):
        require_resolved_hierarchy(run)


def test_soft_policy_parquet_resume_identity_and_schema(tmp_path: Path):
    from dataclasses import replace

    from ogprofiler.hierarchy.stage import hierarchy_component_is_verified

    run = _prepare_multicomponent_run(tmp_path)
    base = hierarchy_config(load_config()["hierarchy"])
    soft = replace(base, resolution=replace(base.resolution, topology_policy="soft_binary_24_v2"))
    run_hierarchy_component_stage(run, 0, base, [])
    assert not hierarchy_component_is_verified(run, 0, soft)
    output, resumed = run_hierarchy_component_stage(run, 0, soft, [])
    assert not resumed
    manifest_path = output / "hierarchy-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    assert (
        manifest["parameters"]["hierarchy"]["resolution"]["topology_policy"] == "soft_binary_24_v2"
    )
    nodes = pq.ParquetFile(output / "nodes.parquet").read().to_pylist()
    candidates = pq.ParquetFile(output / "candidates.parquet").read().to_pylist()
    fallback = [n for n in nodes if n["selection_kind"] == "FALLBACK_KWAY"]
    assert fallback
    assert manifest["fallback_count"] == len(fallback)
    for node in fallback:
        chosen = next(
            c for c in candidates if c["cluster_id"] == node["cluster_id"] and c["selected"]
        )
        assert chosen["kway_eligible"] and not chosen["binary_eligible"]
        assert chosen["original_violations"] == []
        assert chosen["selection_kind"] == node["selection_kind"]
    assert run_hierarchy_component_stage(run, 0, soft, [])[1]
    assert not hierarchy_component_is_verified(run, 0, base)
    manifest["schema_version"] = 2
    manifest_path.write_text(json.dumps(manifest))
    assert not hierarchy_component_is_verified(run, 0, soft)
