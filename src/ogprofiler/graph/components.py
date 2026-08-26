"""Deterministic connected-component extraction and partition export."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ogprofiler.graph.edges import EdgeTable, WeightedEdge
from ogprofiler.graph.union_find import UnionFind


@dataclass(frozen=True, slots=True)
class Component:
    component_id: int
    vertices: tuple[int, ...]
    edges: tuple[WeightedEdge, ...]
    original_ids: tuple[tuple[int, str], ...] = ()

    @property
    def n_vertices(self) -> int:
        return len(self.vertices)

    @property
    def n_edges(self) -> int:
        return len(self.edges)

    def as_edge_table(self) -> EdgeTable:
        return EdgeTable(self.vertices, self.edges, self.original_ids)

    def write(self, directory: Path) -> Path:
        path = directory / f"component_{self.component_id:08d}.parquet"
        self.as_edge_table().write_parquet(path)
        return path


def extract_components(edge_table: EdgeTable) -> tuple[Component, ...]:
    union_find = UnionFind(edge_table.vertices)
    for edge in edge_table.edges:
        union_find.union(edge.source, edge.target)

    vertices_by_root: dict[int, list[int]] = {}
    for vertex in edge_table.vertices:
        vertices_by_root.setdefault(union_find.find(vertex), []).append(vertex)
    edge_by_root: dict[int, list[WeightedEdge]] = {root: [] for root in vertices_by_root}
    for edge in edge_table.edges:
        edge_by_root[union_find.find(edge.source)].append(edge)

    labels = dict(edge_table.original_ids)
    unordered = [
        (
            tuple(sorted(vertices)),
            tuple(sorted(edge_by_root[root])),
        )
        for root, vertices in vertices_by_root.items()
    ]
    unordered.sort(key=lambda item: (-len(item[0]), item[0][0]))
    return tuple(
        Component(
            component_id=index,
            vertices=vertices,
            edges=edges,
            original_ids=tuple((vertex, labels[vertex]) for vertex in vertices if vertex in labels),
        )
        for index, (vertices, edges) in enumerate(unordered)
    )
