from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from ogprofiler.exceptions import ExportError
from ogprofiler.output.results import assemble_export_tables


def _write_fixture(root: Path, membership_order: list[int] | None = None) -> None:
    (root / "input").mkdir(parents=True)
    (root / "components").mkdir()
    hierarchy = root / "hierarchy/components/component=00000000"
    evolution = root / "evolution/components/component=00000000"
    hierarchy.mkdir(parents=True)
    evolution.mkdir(parents=True)
    pq.write_table(
        pa.table(
            {
                "protein_id": [0, 1, 2],
                "species_id": pa.array([0, 1, 2], type=pa.int32()),
                "original_id": ["zeta", "alpha", "singleton"],
            }
        ),
        root / "input/proteins.parquet",
    )
    pq.write_table(
        pa.Table.from_pylist(
            [
                {
                    "cluster_id": 0,
                    "parent_id": None,
                    "component_id": 0,
                    "depth": 0,
                    "n_genes": 2,
                    "n_species": 2,
                    "resolution": None,
                    "quality": None,
                    "child_count": 0,
                    "terminal_reason": "NO_SPLIT",
                }
            ]
        ),
        hierarchy / "nodes.parquet",
    )
    order = membership_order or [0, 1]
    pq.write_table(
        pa.Table.from_pylist(
            [{"protein_id": protein_id, "terminal_cluster_id": 0} for protein_id in order]
        ),
        hierarchy / "members.parquet",
    )
    pq.write_table(
        pa.Table.from_pylist(
            [
                {
                    "cluster_id": 0,
                    "network_event": "AMBIGUOUS",
                    "overlap_score": 0.0,
                    "confidence": 0.0,
                }
            ]
        ),
        evolution / "events.parquet",
    )
    pq.write_table(
        pa.Table.from_pylist(
            [{"protein_id": 2, "component_id": 1, "terminal_reason": "SINGLETON"}]
        ),
        root / "components/singleton_terminal_families.parquet",
    )


def test_family_ids_ignore_membership_row_order_and_include_singletons(tmp_path: Path) -> None:
    _write_fixture(tmp_path, [1, 0])
    first = assemble_export_tables(tmp_path)
    assert [row["family_id"] for row in first.families] == ["OG000000000", "OG000000001"]
    assert [row["n_genes"] for row in first.families] == [2, 1]
    assert first.families[1]["network_event"] == "SPECIES_SPECIFIC"
    assert len(first.hierarchy) == len(first.events) == 2

    pq.write_table(
        pa.Table.from_pylist(
            [{"protein_id": protein_id, "terminal_cluster_id": 0} for protein_id in [0, 1]]
        ),
        tmp_path / "hierarchy/components/component=00000000/members.parquet",
    )
    second = assemble_export_tables(tmp_path)
    assert second.families == first.families
    assert second.members == first.members


def test_export_rejects_incomplete_terminal_partition(tmp_path: Path) -> None:
    _write_fixture(tmp_path)
    singleton_schema = pa.schema(
        [
            ("protein_id", pa.int64()),
            ("component_id", pa.int64()),
            ("terminal_reason", pa.string()),
        ]
    )
    pq.write_table(
        pa.Table.from_pylist([], schema=singleton_schema),
        tmp_path / "components/singleton_terminal_families.parquet",
    )
    with pytest.raises(ExportError, match="without terminal-family assignment"):
        assemble_export_tables(tmp_path)
