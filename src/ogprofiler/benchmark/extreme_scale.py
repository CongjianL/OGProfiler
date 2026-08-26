"""Phase 14 protocols and measurements for extreme-scale execution."""

from __future__ import annotations

import json
import random
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from ogprofiler.core.manifest import sha256_file, write_json


@dataclass(frozen=True, slots=True)
class ComponentTask:
    array_index: int
    component_id: int
    n_vertices: int
    n_edges: int
    estimated_cost: int


def build_component_task_manifest(
    statistics_path: Path, output_path: Path, *, include_singletons: bool = False
) -> tuple[ComponentTask, ...]:
    """Write a stable largest-first Slurm array manifest."""
    rows = pq.read_table(
        statistics_path,
        columns=["component_id", "n_vertices", "n_edges"],
        memory_map=True,
    ).to_pylist()
    eligible = [row for row in rows if include_singletons or int(row["n_vertices"]) > 1]
    eligible.sort(
        key=lambda row: (
            -(int(row["n_vertices"]) + int(row["n_edges"])),
            int(row["component_id"]),
        )
    )
    tasks = tuple(
        ComponentTask(
            array_index=index,
            component_id=int(row["component_id"]),
            n_vertices=int(row["n_vertices"]),
            n_edges=int(row["n_edges"]),
            estimated_cost=int(row["n_vertices"]) + int(row["n_edges"]),
        )
        for index, row in enumerate(eligible)
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pylist([asdict(task) for task in tasks])
    pq.write_table(table, output_path, compression="zstd", row_group_size=65_536)
    write_json(
        output_path.with_suffix(".manifest.json"),
        {
            "algorithm_version": "phase14-component-task-manifest-v1",
            "ordering": "descending(n_vertices+n_edges),ascending(component_id)",
            "include_singletons": include_singletons,
            "task_count": len(tasks),
            "statistics_sha256": sha256_file(statistics_path),
            "output_sha256": sha256_file(output_path),
        },
    )
    return tasks


def merge_component_task_results(
    task_manifest: Path, run_root: Path, output_path: Path
) -> pa.Table:
    """Verify every task artifact and publish the downstream completeness barrier."""
    tasks = pq.read_table(task_manifest, memory_map=True).to_pylist()
    merged: list[dict[str, Any]] = []
    for task in tasks:
        component_id = int(task["component_id"])
        component_root = (
            run_root / "hierarchy" / "components" / f"component={component_id:08d}"
        )
        manifest_path = component_root / "hierarchy-manifest.json"
        if not manifest_path.is_file():
            raise ValueError(f"Missing hierarchy result for component {component_id}")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        checksums = manifest.get("output_checksums")
        if not isinstance(checksums, dict):
            raise ValueError(f"Invalid hierarchy manifest for component {component_id}")
        for name, expected in checksums.items():
            artifact = component_root / str(name)
            if not artifact.is_file() or sha256_file(artifact) != expected:
                raise ValueError(f"Hierarchy checksum mismatch: component {component_id}/{name}")
        metrics = manifest.get("metrics", {})
        merged.append(
            {
                "array_index": int(task["array_index"]),
                "component_id": component_id,
                "n_vertices": int(task["n_vertices"]),
                "n_edges": int(task["n_edges"]),
                "hierarchy_nodes": int(metrics.get("hierarchy_node_count", 0)),
                "terminal_families": int(metrics.get("terminal_family_count", 0)),
                "runtime_seconds": float(metrics.get("runtime_seconds", 0.0)),
                "peak_rss_bytes": int(metrics.get("peak_rss_bytes", 0)),
                "result_path": component_root.relative_to(run_root).as_posix(),
                "manifest_sha256": sha256_file(manifest_path),
            }
        )
    merged.sort(key=lambda row: int(row["array_index"]))
    table = pa.Table.from_pylist(merged)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, output_path, compression="zstd", row_group_size=65_536)
    write_json(
        output_path.with_suffix(".manifest.json"),
        {
            "algorithm_version": "phase14-component-result-merge-v1",
            "task_manifest_sha256": sha256_file(task_manifest),
            "component_count": table.num_rows,
            "output_sha256": sha256_file(output_path),
        },
    )
    return table


@dataclass(frozen=True, slots=True)
class ParquetLayout:
    compression: str
    row_group_size: int


@dataclass(frozen=True, slots=True)
class ParquetProfile:
    compression: str
    row_group_size: int
    file_bytes: int
    row_groups: int
    rows: int
    write_seconds: float
    sequential_seconds: float
    random_seconds: float
    sequential_rows_per_second: float
    random_rows_per_second: float


def profile_parquet_layouts(
    table: pa.Table,
    output_directory: Path,
    layouts: tuple[ParquetLayout, ...],
    *,
    random_row_group_reads: int = 32,
    seed: int = 42,
) -> tuple[ParquetProfile, ...]:
    """Benchmark write, sequential mmap read, and randomized row-group read."""
    if not layouts:
        raise ValueError("At least one Parquet layout is required")
    if random_row_group_reads < 1:
        raise ValueError("random_row_group_reads must be positive")
    output_directory.mkdir(parents=True, exist_ok=True)
    profiles: list[ParquetProfile] = []
    for layout in layouts:
        path = output_directory / f"{layout.compression}-rg{layout.row_group_size}.parquet"
        started = time.perf_counter()
        pq.write_table(
            table,
            path,
            compression=layout.compression,
            row_group_size=layout.row_group_size,
        )
        write_seconds = time.perf_counter() - started

        started = time.perf_counter()
        sequential = pq.read_table(path, memory_map=True)
        sequential_seconds = time.perf_counter() - started

        parquet = pq.ParquetFile(path, memory_map=True)
        rng = random.Random(seed)
        row_group_indices = [
            rng.randrange(parquet.num_row_groups)
            for _ in range(min(random_row_group_reads, max(1, parquet.num_row_groups)))
        ]
        started = time.perf_counter()
        random_rows = sum(parquet.read_row_group(index).num_rows for index in row_group_indices)
        random_seconds = time.perf_counter() - started
        profiles.append(
            ParquetProfile(
                compression=layout.compression,
                row_group_size=layout.row_group_size,
                file_bytes=path.stat().st_size,
                row_groups=parquet.num_row_groups,
                rows=sequential.num_rows,
                write_seconds=write_seconds,
                sequential_seconds=sequential_seconds,
                random_seconds=random_seconds,
                sequential_rows_per_second=sequential.num_rows / max(sequential_seconds, 1e-12),
                random_rows_per_second=random_rows / max(random_seconds, 1e-12),
            )
        )
    return tuple(profiles)


def write_parquet_profiles(path: Path, profiles: tuple[ParquetProfile, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(
        pa.Table.from_pylist([asdict(profile) for profile in profiles]),
        path,
        compression="zstd",
    )
    write_json(
        path.with_suffix(".manifest.json"),
        {
            "algorithm_version": "phase14-parquet-profile-v1",
            "profile_count": len(profiles),
            "output_sha256": sha256_file(path),
        },
    )
