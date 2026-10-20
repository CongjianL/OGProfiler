import json

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from ogprofiler.evolution.stage import run_network_annotation_stage
from ogprofiler.exceptions import HierarchyError


def fixture(root):
    for d in ("input", "components", "hierarchy"):
        (root / d).mkdir()
    rows = [{"protein_id": 0, "component_id": 0}, {"protein_id": 1, "component_id": 1}]
    for name in ("index", "singleton_terminal_families"):
        pq.write_table(pa.Table.from_pylist(rows), root / f"components/{name}.parquet")
    pq.write_table(
        pa.Table.from_pylist([{"protein_id": 0}, {"protein_id": 1}]),
        root / "input/proteins.parquet",
    )
    (root / "hierarchy/scheduler-manifest.json").write_text(
        json.dumps(
            dict(
                hierarchy_status="RESOLVED",
                counts=dict(failed=0, unresolved=0, singleton_components=2),
            )
        )
    )


def test_network_annotation_accepts_verified_all_singleton_dataset(tmp_path):
    fixture(tmp_path)
    path, count, reused = run_network_annotation_stage(tmp_path, 0.5, [])
    report = json.loads(path.read_text())
    assert count == reused == 0
    assert report["component_ids"] == [] and report["event_counts"] == {}


@pytest.mark.parametrize(
    "corruption", ["missing_scheduler", "non_singleton", "missing_singleton", "missing_protein"]
)
def test_empty_hierarchy_does_not_hide_incomplete_inputs(tmp_path, corruption):
    fixture(tmp_path)
    if corruption == "missing_scheduler":
        (tmp_path / "hierarchy/scheduler-manifest.json").unlink()
    elif corruption == "non_singleton":
        pq.write_table(
            pa.Table.from_pylist(
                [{"protein_id": 0, "component_id": 0}, {"protein_id": 1, "component_id": 0}]
            ),
            tmp_path / "components/index.parquet",
        )
    elif corruption == "missing_singleton":
        (tmp_path / "components/singleton_terminal_families.parquet").unlink()
    else:
        pq.write_table(
            pa.Table.from_pylist([{"protein_id": 0}]), tmp_path / "input/proteins.parquet"
        )
    with pytest.raises(HierarchyError):
        run_network_annotation_stage(tmp_path, 0.5, [])
