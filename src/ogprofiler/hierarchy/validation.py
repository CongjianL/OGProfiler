"""Structural invariant checks for hierarchy results."""

from __future__ import annotations

from collections import defaultdict

from ogprofiler.exceptions import HierarchyError
from ogprofiler.graph.components import Component
from ogprofiler.hierarchy.engine import HierarchyResult


def validate_hierarchy(component: Component, result: HierarchyResult) -> None:
    nodes = {node.cluster_id: node for node in result.nodes}
    if len(nodes) != len(result.nodes):
        raise HierarchyError("Hierarchy contains duplicate cluster IDs")
    roots = [node for node in result.nodes if node.parent_id is None]
    if len(roots) != 1:
        raise HierarchyError("Hierarchy must contain exactly one root")
    children: dict[int, list[int]] = defaultdict(list)
    for node in result.nodes:
        if node.component_id != component.component_id:
            raise HierarchyError(
                f"Node {node.cluster_id} has component {node.component_id}, "
                f"expected {component.component_id}"
            )
        if node.parent_id is not None:
            if node.parent_id not in nodes:
                raise HierarchyError(f"Missing parent {node.parent_id} for node {node.cluster_id}")
            if node.depth != nodes[node.parent_id].depth + 1:
                raise HierarchyError(f"Invalid depth for node {node.cluster_id}")
            children[node.parent_id].append(node.cluster_id)
        elif node.depth != 0:
            raise HierarchyError("Hierarchy root must have depth zero")

    for node in result.nodes:
        seen: set[int] = set()
        current = node
        while current.parent_id is not None:
            if current.cluster_id in seen:
                raise HierarchyError("Hierarchy contains a cycle")
            seen.add(current.cluster_id)
            current = nodes[current.parent_id]

    membership = dict(result.terminal_membership)
    if len(membership) != len(result.terminal_membership):
        raise HierarchyError("Duplicate protein membership")
    if set(membership) != set(component.vertices):
        raise HierarchyError("Every component protein must have exactly one terminal membership")
    terminal_ids = {node.cluster_id for node in result.nodes if node.terminal_reason is not None}
    if set(membership.values()) - terminal_ids:
        raise HierarchyError("Terminal membership points to a non-terminal node")

    descendants: dict[int, set[int]] = defaultdict(set)
    for protein_id, terminal_id in membership.items():
        current_id: int | None = terminal_id
        while current_id is not None:
            descendants[current_id].add(protein_id)
            current_id = nodes[current_id].parent_id
    for parent_id, child_ids in children.items():
        if len(child_ids) < 2:
            raise HierarchyError(f"Split node {parent_id} has fewer than two children")
        child_sets = [descendants[child_id] for child_id in child_ids]
        combined: set[int] = set()
        for child_set in child_sets:
            if combined & child_set:
                raise HierarchyError(f"Children of node {parent_id} overlap")
            combined.update(child_set)
        if combined != descendants[parent_id]:
            raise HierarchyError(f"Children of node {parent_id} do not cover the parent")
    for node in result.nodes:
        if node.child_count != len(children[node.cluster_id]):
            raise HierarchyError(f"Incorrect child_count for node {node.cluster_id}")
        if node.n_genes != len(descendants[node.cluster_id]):
            raise HierarchyError(f"Incorrect n_genes for node {node.cluster_id}")

    for candidate in result.resolution_candidates:
        if candidate.selected and candidate.selection_kind == "FALLBACK_KWAY":
            if (
                not candidate.kway_eligible
                or candidate.binary_eligible
                or candidate.original_violations
                or candidate.child_count <= 2
                or not candidate.phase.startswith("FALLBACK_KWAY/")
            ):
                raise HierarchyError("Fallback candidate lacks original-gated k-way evidence")
