from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pyarrow.parquet as pq

from ogprofiler.cli import main
from ogprofiler.core.checkpoint import CheckpointStore
from ogprofiler.similarity.io import write_retained_edges
from ogprofiler.similarity.models import RetainedEdge


def _prepare_multicomponent_run(tmp_path: Path) -> Path:
    proteomes = tmp_path / "proteomes"
    proteomes.mkdir()
    for species in range(3):
        records = "".join(f">g{index}\n{'A' * 20}\n" for index in range(6))
        (proteomes / f"species_{species}.faa").write_text(records, encoding="utf-8")
    run = tmp_path / "run"
    assert main(["prepare", "--proteomes", str(proteomes), "--out", str(run)]) == 0
    species_by_id = {protein_id: protein_id // 6 for protein_id in range(18)}
    groups = [(0, 1, 6, 7, 12, 13), (2, 3, 8, 9)]
    edges: list[RetainedEdge] = []
    for group in groups:
        for offset, left in enumerate(group):
            for right in group[offset + 1 :]:
                edges.append(
                    RetainedEdge(
                        left,
                        right,
                        species_by_id[left],
                        species_by_id[right],
                        1.0,
                        1.0,
                        1.0,
                        100.0,
                        "TEST",
                    )
                )
    write_retained_edges(run / "edges" / "retained_edges.parquet", edges)
    assert main(["components", "--run", str(run)]) == 0
    return run


def _command(run: Path, seed: int = 42) -> list[str]:
    return [
        "hierarchy-all",
        "--run",
        str(run),
        "--set",
        "runtime.workers=2",
        "--set",
        "runtime.component_retries=0",
        "--set",
        "hierarchy.stability_mode=fast",
        "--set",
        "hierarchy.gamma_max=0.1",
        "--set",
        f"hierarchy.seed={seed}",
    ]


def test_scheduler_parallel_resume_corruption_and_config_invalidation(tmp_path: Path) -> None:
    run = _prepare_multicomponent_run(tmp_path)
    command = _command(run)
    assert main(command) == 0
    scheduler_path = run / "hierarchy" / "scheduler-manifest.json"
    first = json.loads(scheduler_path.read_text(encoding="utf-8"))
    assert first["counts"] == {
        "completed": 2,
        "failed": 0,
        "scheduled": 2,
        "singleton_components": 8,
        "skipped": 0,
    }
    statistics = pq.read_table(run / "components" / "statistics.parquet").to_pylist()
    assert [row["n_vertices"] for row in statistics[:2]] == [6, 4]

    assert main(command) == 0
    resumed = json.loads(scheduler_path.read_text(encoding="utf-8"))
    assert resumed["counts"]["scheduled"] == 0
    assert resumed["counts"]["skipped"] == 2

    store = CheckpointStore(run / "run.db")
    store.invalidate("hierarchy", "0", "simulate scheduler interruption")
    store.start("hierarchy", "0")
    assert main(command) == 0
    recovered = json.loads(scheduler_path.read_text(encoding="utf-8"))
    assert recovered["counts"]["scheduled"] == 0
    assert recovered["counts"]["skipped"] == 2
    assert store.get("hierarchy", "0").status == "DONE"  # type: ignore[union-attr]

    damaged = run / "hierarchy" / "components" / "component=00000000" / "nodes.parquet"
    damaged.write_bytes(b"corrupt")
    assert main(command) == 0
    rebuilt = json.loads(scheduler_path.read_text(encoding="utf-8"))
    assert rebuilt["counts"]["scheduled"] == 1
    assert rebuilt["counts"]["skipped"] == 1
    assert pq.read_table(damaged).num_rows >= 1

    assert main(_command(run, seed=99)) == 0
    changed = json.loads(scheduler_path.read_text(encoding="utf-8"))
    assert changed["counts"]["scheduled"] == 2
    with sqlite3.connect(run / "run.db") as database:
        invalid_events = database.execute(
            "SELECT COUNT(*) FROM task_events WHERE stage='hierarchy' AND status='INVALID'"
        ).fetchone()[0]
    assert invalid_events >= 2


def test_worker_failure_isolated_and_failed_only_resume(tmp_path: Path) -> None:
    run = _prepare_multicomponent_run(tmp_path)
    broken = next((run / "components" / "edges" / "component=00000001").glob("*.parquet"))
    original = broken.read_bytes()
    broken.write_bytes(b"corrupt")
    command = _command(run)
    assert main(command) == 2
    scheduler_path = run / "hierarchy" / "scheduler-manifest.json"
    failed = json.loads(scheduler_path.read_text(encoding="utf-8"))
    assert failed["counts"]["completed"] == 1
    assert failed["counts"]["failed"] == 1
    assert failed["failed_component_ids"] == [1]
    assert (run / "hierarchy" / "components" / "component=00000000").is_dir()

    broken.write_bytes(original)
    retry = [*_command(run), "--failed-only"]
    assert main(retry) == 0
    recovered = json.loads(scheduler_path.read_text(encoding="utf-8"))
    assert recovered["counts"]["scheduled"] == 1
    assert recovered["counts"]["completed"] == 1
    assert recovered["counts"]["failed"] == 0
    assert recovered["counts"]["skipped"] == 1
