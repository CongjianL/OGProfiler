"""Compact canonical undirected edge-table representation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import igraph as ig
import pyarrow as pa
import pyarrow.parquet as pq

from ogprofiler.exceptions import EdgeConstructionError


@dataclass(frozen=True, order=True, slots=True)
class WeightedEdge:
    source: int
    target: int
    weight: float

    def __post_init__(self) -> None:
        if self.source >= self.target:
            raise EdgeConstructionError(
                f"Undirected edges must be canonical source < target: {self.source}, {self.target}"
            )


@dataclass(frozen=True, slots=True)
class EdgeTable:
    vertices: tuple[int, ...]
    edges: tuple[WeightedEdge, ...]
    original_ids: tuple[tuple[int, str], ...] = ()

    def __post_init__(self) -> None:
        if tuple(sorted(set(self.vertices))) != self.vertices:
            raise EdgeConstructionError("EdgeTable vertices must be sorted and unique")
        vertex_set = set(self.vertices)
        seen_pairs: set[tuple[int, int]] = set()
        for edge in self.edges:
            if edge.source not in vertex_set or edge.target not in vertex_set:
                raise EdgeConstructionError("Edge endpoint is absent from EdgeTable vertices")
            pair = (edge.source, edge.target)
            if pair in seen_pairs:
                raise EdgeConstructionError(f"Duplicate canonical edge: {pair}")
            seen_pairs.add(pair)
        if tuple(sorted(self.edges)) != self.edges:
            raise EdgeConstructionError("EdgeTable edges must be sorted")

    @classmethod
    def canonicalize(
        cls,
        vertices: list[int] | tuple[int, ...],
        edges: list[tuple[int, int, float]],
        original_ids: dict[int, str] | None = None,
    ) -> EdgeTable:
        merged: dict[tuple[int, int], float] = {}
        for left, right, weight in edges:
            if left == right:
                continue
            source, target = sorted((int(left), int(right)))
            pair = (source, target)
            merged[pair] = max(merged.get(pair, float("-inf")), float(weight))
        weighted = tuple(
            WeightedEdge(source, target, weight)
            for (source, target), weight in sorted(merged.items())
        )
        labels = tuple(sorted((original_ids or {}).items()))
        return cls(tuple(sorted(set(vertices))), weighted, labels)

    def to_igraph(self) -> tuple[ig.Graph, tuple[int, ...]]:
        global_ids = self.vertices
        local_by_global = {global_id: local for local, global_id in enumerate(global_ids)}
        graph = ig.Graph(
            n=len(global_ids),
            edges=[
                (local_by_global[edge.source], local_by_global[edge.target]) for edge in self.edges
            ],
            directed=False,
        )
        graph.vs["protein_id"] = list(global_ids)
        graph.es["weight"] = [edge.weight for edge in self.edges]
        return graph, global_ids

    def write_parquet(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        edge_table = pa.table(
            {
                "source": pa.array([edge.source for edge in self.edges], type=pa.int64()),
                "target": pa.array([edge.target for edge in self.edges], type=pa.int64()),
                "weight": pa.array([edge.weight for edge in self.edges], type=pa.float64()),
            }
        )
        pq.write_table(edge_table, path, compression="zstd")
        labels = dict(self.original_ids)
        vertex_table = pa.table(
            {
                "protein_id": pa.array(self.vertices, type=pa.int64()),
                "original_id": pa.array(
                    [labels.get(vertex) for vertex in self.vertices], type=pa.string()
                ),
            }
        )
        pq.write_table(vertex_table, path.with_suffix(".vertices.parquet"), compression="zstd")

    @classmethod
    def read_parquet(cls, path: Path) -> EdgeTable:
        edge_rows = pq.read_table(path).to_pylist()
        vertex_rows = pq.read_table(path.with_suffix(".vertices.parquet")).to_pylist()
        original_ids = {
            int(row["protein_id"]): str(row["original_id"])
            for row in vertex_rows
            if row["original_id"] is not None
        }
        return cls.canonicalize(
            [int(row["protein_id"]) for row in vertex_rows],
            [(int(row["source"]), int(row["target"]), float(row["weight"])) for row in edge_rows],
            original_ids,
        )
