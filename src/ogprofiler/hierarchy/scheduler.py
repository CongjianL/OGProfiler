"""Largest-first process scheduler with SQLite-owned resume semantics."""

from __future__ import annotations

import json
import multiprocessing
import traceback
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path

import pyarrow.parquet as pq

from ogprofiler.core.checkpoint import CheckpointStore, TaskRecord
from ogprofiler.core.manifest import sha256_json, write_json
from ogprofiler.core.workspace import Workspace
from ogprofiler.graph.partition import non_singleton_component_ids
from ogprofiler.hierarchy.engine import HierarchyConfig
from ogprofiler.hierarchy.stage import (
    HIERARCHY_ALGORITHM_VERSION,
    hierarchy_component_identity,
    hierarchy_component_is_verified,
    run_hierarchy_component_stage,
)

HIERARCHY_STAGE = "hierarchy"


@dataclass(frozen=True, slots=True)
class WorkerSummary:
    component_id: int
    status: str
    output_path: str | None
    reused: bool
    hierarchy_nodes: int
    terminal_families: int
    runtime_seconds: float
    error: str | None = None


@dataclass(frozen=True, slots=True)
class SchedulerResult:
    scheduled: int
    completed: int
    skipped: int
    failed: int
    singleton_components: int
    summaries: tuple[WorkerSummary, ...]


_WORKER_RUN_ROOT: Path | None = None
_WORKER_CONFIG: HierarchyConfig | None = None
_WORKER_COMMAND: list[str] | None = None


def _initialize_worker(run_root: str, config: HierarchyConfig, command: list[str]) -> None:
    global _WORKER_RUN_ROOT, _WORKER_CONFIG, _WORKER_COMMAND
    _WORKER_RUN_ROOT = Path(run_root)
    _WORKER_CONFIG = config
    _WORKER_COMMAND = command


def _run_worker(component_id: int) -> WorkerSummary:
    if _WORKER_RUN_ROOT is None or _WORKER_CONFIG is None or _WORKER_COMMAND is None:
        raise RuntimeError("Hierarchy worker was not initialized")
    try:
        output, reused = run_hierarchy_component_stage(
            _WORKER_RUN_ROOT, component_id, _WORKER_CONFIG, _WORKER_COMMAND
        )
        metrics = json.loads((output / "metrics.json").read_text(encoding="utf-8"))
        return WorkerSummary(
            component_id=component_id,
            status="DONE",
            output_path=str(output),
            reused=reused,
            hierarchy_nodes=int(metrics["hierarchy_node_count"]),
            terminal_families=int(metrics["terminal_family_count"]),
            runtime_seconds=float(metrics["runtime_seconds"]),
        )
    except Exception as error:  # worker boundary intentionally captures component failure
        detail = "".join(traceback.format_exception_only(type(error), error)).strip()
        return WorkerSummary(component_id, "FAILED", None, False, 0, 0, 0.0, detail)


def _component_sizes(run_root: Path) -> dict[int, int]:
    rows = pq.read_table(
        run_root / "components" / "statistics.parquet",
        columns=["component_id", "n_vertices"],
    ).to_pylist()
    return {int(row["component_id"]): int(row["n_vertices"]) for row in rows}


def _runnable(record: TaskRecord, retries: int, failed_only: bool) -> bool:
    if failed_only and record.status not in {"FAILED", "INVALID"}:
        return False
    if record.status in {"PENDING", "INVALID"}:
        return True
    return record.status == "FAILED" and record.attempts <= retries


def run_hierarchy_scheduler(
    run_root: Path,
    config: HierarchyConfig,
    command: list[str],
    *,
    workers: int,
    retries: int,
    failed_only: bool = False,
) -> SchedulerResult:
    workspace = Workspace.create(run_root)
    store = CheckpointStore(workspace.database_path)
    sizes = _component_sizes(run_root)
    component_ids = non_singleton_component_ids(run_root / "components")
    singleton_count = sum(1 for size in sizes.values() if size == 1)
    pending: list[int] = []
    skipped = 0
    for component_id in component_ids:
        identity, _, _ = hierarchy_component_identity(run_root, component_id, config)
        record = store.register(
            HIERARCHY_STAGE, str(component_id), identity, HIERARCHY_ALGORITHM_VERSION
        )
        if record.status == "DONE":
            if hierarchy_component_is_verified(run_root, component_id, config):
                skipped += 1
                continue
            store.invalidate(HIERARCHY_STAGE, str(component_id), "output checksum mismatch")
            record = store.get(HIERARCHY_STAGE, str(component_id)) or record
        elif hierarchy_component_is_verified(run_root, component_id, config):
            output = run_root / "hierarchy" / "components" / f"component={component_id:08d}"
            store.finish(HIERARCHY_STAGE, str(component_id), str(output))
            skipped += 1
            continue
        if _runnable(record, retries, failed_only):
            pending.append(component_id)

    pending.sort(key=lambda component_id: (-sizes[component_id], component_id))
    summaries: list[WorkerSummary] = []
    remaining = pending
    while remaining:
        current = list(remaining)
        remaining = []
        for component_id in current:
            store.start(HIERARCHY_STAGE, str(component_id))
        context = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(
            max_workers=workers,
            mp_context=context,
            initializer=_initialize_worker,
            initargs=(str(run_root), config, command),
        ) as executor:
            futures: dict[Future[WorkerSummary], int] = {
                executor.submit(_run_worker, component_id): component_id
                for component_id in current
            }
            for future in as_completed(futures):
                component_id = futures[future]
                try:
                    summary = future.result()
                except Exception as error:
                    detail = "".join(traceback.format_exception_only(type(error), error)).strip()
                    summary = WorkerSummary(
                        component_id, "FAILED", None, False, 0, 0, 0.0, detail
                    )
                if summary.status == "DONE" and summary.output_path is not None:
                    store.finish(HIERARCHY_STAGE, str(component_id), summary.output_path)
                    summaries.append(summary)
                else:
                    store.fail(
                        HIERARCHY_STAGE,
                        str(component_id),
                        summary.error or "unknown worker failure",
                    )
                    retry_record = store.get(HIERARCHY_STAGE, str(component_id))
                    if retry_record is not None and retry_record.attempts <= retries:
                        remaining.append(component_id)
                    else:
                        summaries.append(summary)
        remaining.sort(key=lambda component_id: (-sizes[component_id], component_id))

    records = store.list(HIERARCHY_STAGE)
    current_ids = set(component_ids)
    failed_ids = sorted(
        int(record.task_id)
        for record in records
        if record.status == "FAILED" and int(record.task_id) in current_ids
    )
    result = SchedulerResult(
        scheduled=len(pending),
        completed=sum(summary.status == "DONE" for summary in summaries),
        skipped=skipped,
        failed=len(failed_ids),
        singleton_components=singleton_count,
        summaries=tuple(sorted(summaries, key=lambda summary: summary.component_id)),
    )
    manifest = {
        "algorithm_version": "hierarchy-scheduler-v1",
        "command": command,
        "parameters": {
            "workers": workers,
            "retries": retries,
            "failed_only": failed_only,
            "hierarchy": asdict(config),
        },
        "parameters_sha256": sha256_json(
            {"workers": workers, "retries": retries, "hierarchy": asdict(config)}
        ),
        "counts": {
            "scheduled": result.scheduled,
            "completed": result.completed,
            "skipped": result.skipped,
            "failed": result.failed,
            "singleton_components": singleton_count,
        },
        "failed_component_ids": failed_ids,
        "summaries": [asdict(summary) for summary in result.summaries],
    }
    write_json(run_root / "hierarchy" / "scheduler-manifest.json", manifest)
    return result
