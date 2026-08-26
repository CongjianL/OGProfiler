from __future__ import annotations

import random
from pathlib import Path

import igraph as ig
import pyarrow as pa
import pyarrow.parquet as pq

from ogprofiler.graph.partition import (
    build_component_artifacts,
    load_component_edge_table,
    non_singleton_component_ids,
    read_component_edges,
)
from ogprofiler.graph.union_find import UnionFind
from ogprofiler.similarity.io import RETAINED_EDGE_SCHEMA


def _reference(vertices: list[int], edges: list[tuple[int, int]]) -> set[frozenset[int]]:
    adjacency = {vertex: set() for vertex in vertices}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    remaining = set(vertices)
    groups: set[frozenset[int]] = set()
    while remaining:
        stack, seen = [min(remaining)], set()
        while stack:
            vertex = stack.pop()
            if vertex not in seen:
                seen.add(vertex)
                stack.extend(adjacency[vertex] - seen)
        remaining -= seen
        groups.add(frozenset(seen))
    return groups


def test_union_find_random_graph_matches_reference_components() -> None:
    for seed in range(20):
        rng = random.Random(seed)
        vertices = list(range(40))
        edges = [
            (left, right)
            for left in vertices
            for right in range(left + 1, len(vertices))
            if rng.random() < 0.035
        ]
        union_find = UnionFind(vertices)
        for left, right in edges:
            union_find.union(left, right)
        observed: dict[int, set[int]] = {}
        for vertex in vertices:
            observed.setdefault(union_find.find(vertex), set()).add(vertex)
        assert {frozenset(group) for group in observed.values()} == _reference(vertices, edges)
        graph = ig.Graph(n=len(vertices), edges=edges, directed=False)
        expected_igraph = {
            frozenset(group) for group in graph.connected_components(mode="weak")
        }
        assert {frozenset(group) for group in observed.values()} == expected_igraph
        assert all(union_find.parent[vertex] == union_find.find(vertex) for vertex in vertices)


def _write_fixture(root: Path) -> tuple[Path, Path]:
    proteins = root / "proteins.parquet"
    edges = root / "retained_edges.parquet"
    pq.write_table(
        pa.table(
            {
                "protein_id": pa.array(range(7), type=pa.int64()),
                "species_id": pa.array([0, 1, 2, 0, 1, 2, 3], type=pa.int32()),
            }
        ),
        proteins,
    )
    rows = []
    for left, right, weight in [(3, 4, 0.7), (0, 1, 0.9), (1, 2, 0.8), (4, 5, 0.6)]:
        rows.append(
            {
                "u": left,
                "v": right,
                "u_species": left % 3,
                "v_species": right % 3,
                "score_uv": weight,
                "score_vu": weight,
                "weight": weight,
                "coverage": 100.0,
                "edge_type": "LRB",
            }
        )
    pq.write_table(pa.Table.from_pylist(rows, schema=RETAINED_EDGE_SCHEMA), edges)
    return proteins, edges


def test_streamed_partitions_are_largest_first_and_preserve_singletons(tmp_path: Path) -> None:
    proteins, edges = _write_fixture(tmp_path)
    output = tmp_path / "components"
    result = build_component_artifacts(proteins, edges, output, batch_size=1, max_open_files=1)
    assert (result.component_count, result.singleton_count, result.edge_count) == (3, 1, 4)
    assert pq.read_table(output / "index.parquet").to_pylist() == [
        {
            "protein_id": protein_id,
            "component_id": 0 if protein_id < 3 else 1 if protein_id < 6 else 2,
        }
        for protein_id in range(7)
    ]
    assert pq.read_table(output / "statistics.parquet").to_pylist() == [
        {"component_id": 0, "n_vertices": 3, "n_edges": 2, "n_species": 3},
        {"component_id": 1, "n_vertices": 3, "n_edges": 2, "n_species": 3},
        {"component_id": 2, "n_vertices": 1, "n_edges": 0, "n_species": 1},
    ]
    assert pq.read_table(output / "singleton_terminal_families.parquet").to_pylist() == [
        {"protein_id": 6, "component_id": 2, "terminal_reason": "SINGLETON"}
    ]
    assert read_component_edges(output, 0).num_rows == 2
    assert read_component_edges(output, 2).num_rows == 0
    assert load_component_edge_table(output, 1).vertices == (3, 4, 5)
    assert non_singleton_component_ids(output) == (0, 1)
