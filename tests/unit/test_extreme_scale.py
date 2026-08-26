from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from ogprofiler.benchmark.extreme_scale import (
    ParquetLayout,
    build_component_task_manifest,
    merge_component_task_results,
    profile_parquet_layouts,
)
from ogprofiler.hierarchy.loader import LocalArrayLookup
from ogprofiler.hierarchy.subtree import SubtreePath


def test_subtree_paths_are_deterministic_and_depth_first_sortable() -> None:
    root = SubtreePath(7)
    left = root.child(0)
    right = root.child(1)
    grandchild = left.child(2)
    assert root.task_id == "component-00000007/root"
    assert grandchild.task_id == "component-00000007/000000.000002"
    assert grandchild.parent == left
    assert sorted((right, grandchild, root, left)) == [root, left, grandchild, right]


def test_local_array_lookup_has_mapping_semantics_without_python_dict() -> None:
    lookup = LocalArrayLookup(
        np.asarray([10, 20, 30], dtype=np.int64),
        np.asarray([4, 5, 6], dtype=np.int32),
    )
    assert lookup[20] == 5
    assert list(lookup) == [10, 20, 30]
    assert dict(lookup) == {10: 4, 20: 5, 30: 6}


def test_component_task_manifest_is_largest_first_and_checksummed(tmp_path: Path) -> None:
    statistics = tmp_path / "statistics.parquet"
    pq.write_table(
        pa.table(
            {
                "component_id": [0, 1, 2],
                "n_vertices": [10, 20, 1],
                "n_edges": [30, 5, 0],
            }
        ),
        statistics,
    )
    output = tmp_path / "tasks.parquet"
    tasks = build_component_task_manifest(statistics, output)
    assert [task.component_id for task in tasks] == [0, 1]
    assert [task.array_index for task in tasks] == [0, 1]
    assert output.with_suffix(".manifest.json").is_file()


def test_parquet_profile_covers_sequential_and_random_access(tmp_path: Path) -> None:
    table = pa.table({"u": range(100), "v": range(1, 101), "weight": [0.5] * 100})
    profiles = profile_parquet_layouts(
        table,
        tmp_path,
        (ParquetLayout("zstd", 16),),
        random_row_group_reads=3,
    )
    assert len(profiles) == 1
    assert profiles[0].rows == 100
    assert profiles[0].row_groups == 7
    assert profiles[0].file_bytes > 0
    assert profiles[0].sequential_rows_per_second > 0
    assert profiles[0].random_rows_per_second > 0


def test_component_result_merge_is_a_checksum_completeness_barrier(tmp_path: Path) -> None:
    statistics = tmp_path / "statistics.parquet"
    pq.write_table(
        pa.table({"component_id": [0], "n_vertices": [3], "n_edges": [2]}), statistics
    )
    tasks = tmp_path / "tasks.parquet"
    build_component_task_manifest(statistics, tasks)
    component = tmp_path / "run/hierarchy/components/component=00000000"
    component.mkdir(parents=True)
    artifact = component / "nodes.parquet"
    artifact.write_bytes(b"nodes")
    from ogprofiler.core.manifest import sha256_file

    (component / "hierarchy-manifest.json").write_text(
        json.dumps(
            {
                "output_checksums": {"nodes.parquet": sha256_file(artifact)},
                "metrics": {
                    "hierarchy_node_count": 2,
                    "terminal_family_count": 1,
                    "runtime_seconds": 0.25,
                    "peak_rss_bytes": 1000,
                },
            }
        ),
        encoding="utf-8",
    )
    output = tmp_path / "merged.parquet"
    merged = merge_component_task_results(tasks, tmp_path / "run", output)
    assert merged.to_pylist()[0]["hierarchy_nodes"] == 2
    artifact.write_bytes(b"corrupt")
    try:
        merge_component_task_results(tasks, tmp_path / "run", output)
    except ValueError as error:
        assert "checksum mismatch" in str(error)
    else:
        raise AssertionError("corrupt component result was accepted")
