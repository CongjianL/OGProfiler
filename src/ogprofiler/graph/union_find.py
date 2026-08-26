"""Streaming integer Union-Find with path compression and union by size."""

from __future__ import annotations

from collections.abc import Iterable

from ogprofiler.exceptions import ComponentError


class UnionFind:
    def __init__(self, vertices: Iterable[int]) -> None:
        ordered = tuple(vertices)
        if any(not isinstance(vertex, int) for vertex in ordered):
            raise ComponentError("Union-Find vertices must be integer protein IDs")
        if len(set(ordered)) != len(ordered):
            raise ComponentError("Union-Find vertices must be unique")
        self.parent = {vertex: vertex for vertex in ordered}
        self.size = {vertex: 1 for vertex in ordered}

    def find(self, vertex: int) -> int:
        parent = self.parent[vertex]
        while parent != self.parent[parent]:
            parent = self.parent[parent]
        while vertex != parent:
            next_vertex = self.parent[vertex]
            self.parent[vertex] = parent
            vertex = next_vertex
        return parent

    def union(self, left: int, right: int) -> int:
        if left not in self.parent or right not in self.parent:
            raise ComponentError(f"Edge endpoint is absent from protein index: {left}, {right}")
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return left_root
        if self.size[left_root] < self.size[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        self.size[left_root] += self.size[right_root]
        return left_root
