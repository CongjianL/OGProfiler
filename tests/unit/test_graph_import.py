from __future__ import annotations

from pathlib import Path

import igraph as ig
import pytest

from ogprofiler.exceptions import InputError
from ogprofiler.graph.components import extract_components
from ogprofiler.graph.edges import EdgeTable
from ogprofiler.graph.legacy import import_legacy_ssn


def _write_gml(path: Path, weight_attribute: str = "NBS") -> None:
    graph = ig.Graph(n=6, edges=[(0, 1), (1, 2), (3, 4)], directed=False)
    graph.vs["name"] = ["G1|g1", "G1|g0", "G0|g0", "G2|g0", "G2|g1", "singleton"]
    graph.es[weight_attribute] = [0.2, 0.8, 0.5]
    graph.write_gml(str(path))


def test_import_legacy_gml_assigns_stable_integer_ids(tmp_path: Path) -> None:
    path = tmp_path / "ssn.gml"
    _write_gml(path)
    table = import_legacy_ssn(path)
    assert table.vertices == (0, 1, 2, 3, 4, 5)
    assert dict(table.original_ids) == {
        0: "G0|g0",
        1: "G1|g0",
        2: "G1|g1",
        3: "G2|g0",
        4: "G2|g1",
        5: "singleton",
    }
    assert [(edge.source, edge.target, edge.weight) for edge in table.edges] == [
        (0, 1, 0.8),
        (1, 2, 0.2),
        (3, 4, 0.5),
    ]


def test_missing_weight_attribute_is_reported(tmp_path: Path) -> None:
    path = tmp_path / "ssn.gml"
    _write_gml(path, "score")
    with pytest.raises(InputError, match="NBS"):
        import_legacy_ssn(path)


def test_components_are_largest_first_and_include_singletons(tmp_path: Path) -> None:
    path = tmp_path / "ssn.gml"
    _write_gml(path)
    components = extract_components(import_legacy_ssn(path))
    assert [(item.component_id, item.vertices, item.n_edges) for item in components] == [
        (0, (0, 1, 2), 2),
        (1, (3, 4), 1),
        (2, (5,), 0),
    ]


def test_edge_table_parquet_round_trip_preserves_isolates(tmp_path: Path) -> None:
    table = EdgeTable.canonicalize(
        [10, 20, 30],
        [(20, 10, 0.2), (10, 20, 0.8)],
        {10: "a", 20: "b", 30: "c"},
    )
    path = tmp_path / "edges.parquet"
    table.write_parquet(path)
    assert EdgeTable.read_parquet(path) == EdgeTable.canonicalize(
        [10, 20, 30], [(10, 20, 0.8)], {10: "a", 20: "b", 30: "c"}
    )
