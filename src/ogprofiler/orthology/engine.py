"""Ortholog candidates from cross-child genes at speciation-like nodes."""

from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from ogprofiler.exceptions import HierarchyError

ORTHOLOGY_SUPPORTING_EVENTS = frozenset({"SPECIATION_LIKE", "POLYTOMY"})


@dataclass(frozen=True, slots=True)
class OrthologCandidate:
    protein_a_id: int
    protein_b_id: int
    species_a_id: int
    species_b_id: int
    component_id: int
    supporting_cluster_id: int
    relationship: str = "CO_ORTHOLOG_CANDIDATE"


def _tree_indexes(
    nodes: Sequence[Mapping[str, Any]],
) -> tuple[dict[int, tuple[int, ...]], set[int]]:
    node_ids = {int(row["cluster_id"]) for row in nodes}
    if len(node_ids) != len(nodes):
        raise HierarchyError("Duplicate cluster IDs in orthology hierarchy")
    children: dict[int, list[int]] = {cluster_id: [] for cluster_id in node_ids}
    roots: set[int] = set()
    for row in nodes:
        cluster_id = int(row["cluster_id"])
        parent = row.get("parent_id")
        if parent is None:
            roots.add(cluster_id)
            continue
        parent_id = int(parent)
        if parent_id not in node_ids:
            raise HierarchyError(f"Unknown parent cluster {parent_id}")
        children[parent_id].append(cluster_id)
    if len(roots) != 1:
        raise HierarchyError("Component hierarchy must have exactly one root")
    result = {key: tuple(sorted(value)) for key, value in children.items()}
    visiting: set[int] = set()
    visited: set[int] = set()

    def visit(cluster_id: int) -> None:
        if cluster_id in visiting:
            raise HierarchyError("Cycle in component hierarchy")
        if cluster_id in visited:
            return
        visiting.add(cluster_id)
        for child_id in result[cluster_id]:
            visit(child_id)
        visiting.remove(cluster_id)
        visited.add(cluster_id)

    visit(next(iter(roots)))
    if visited != node_ids:
        raise HierarchyError("Disconnected nodes in component hierarchy")
    return result, node_ids


def generate_component_candidates(
    *,
    component_id: int,
    nodes: Sequence[Mapping[str, Any]],
    terminal_memberships: Sequence[Mapping[str, Any]],
    events: Sequence[Mapping[str, Any]],
    species_by_protein: Mapping[int, int],
) -> Iterator[OrthologCandidate]:
    children, node_ids = _tree_indexes(nodes)
    terminal_members: dict[int, list[int]] = {}
    seen_proteins: set[int] = set()
    for row in terminal_memberships:
        protein_id = int(row["protein_id"])
        cluster_id = int(row["terminal_cluster_id"])
        if cluster_id not in node_ids:
            raise HierarchyError(f"Membership references unknown cluster {cluster_id}")
        if protein_id in seen_proteins:
            raise HierarchyError(f"Duplicate terminal membership for protein {protein_id}")
        if protein_id not in species_by_protein:
            raise HierarchyError(f"Missing species metadata for protein {protein_id}")
        seen_proteins.add(protein_id)
        terminal_members.setdefault(cluster_id, []).append(protein_id)
    for members in terminal_members.values():
        members.sort()

    event_by_node = {int(row["cluster_id"]): str(row["network_event"]) for row in events}
    if set(event_by_node) != node_ids:
        raise HierarchyError("Hierarchy/event node mismatch during orthology traversal")

    def descendant_genes(cluster_id: int) -> list[int]:
        child_ids = children[cluster_id]
        if not child_ids:
            members = terminal_members.get(cluster_id)
            if not members:
                raise HierarchyError(f"Terminal cluster {cluster_id} has no members")
            return members
        genes: list[int] = []
        for child_id in child_ids:
            genes.extend(descendant_genes(child_id))
        return genes

    for cluster_id in sorted(node_ids):
        if event_by_node[cluster_id] not in ORTHOLOGY_SUPPORTING_EVENTS:
            continue
        child_ids = children[cluster_id]
        if len(child_ids) < 2:
            raise HierarchyError(
                f"Orthology-supporting node {cluster_id} has fewer than two children"
            )
        child_genes = [descendant_genes(child_id) for child_id in child_ids]
        for left_index, left_genes in enumerate(child_genes):
            for right_genes in child_genes[left_index + 1 :]:
                for left in left_genes:
                    for right in right_genes:
                        left_species = species_by_protein[left]
                        right_species = species_by_protein[right]
                        if left_species == right_species:
                            continue
                        if left < right:
                            yield OrthologCandidate(
                                left,
                                right,
                                left_species,
                                right_species,
                                component_id,
                                cluster_id,
                            )
                        else:
                            yield OrthologCandidate(
                                right,
                                left,
                                right_species,
                                left_species,
                                component_id,
                                cluster_id,
                            )


def generate_ortholog_candidates(
    components: Sequence[
        tuple[
            int,
            Sequence[Mapping[str, Any]],
            Sequence[Mapping[str, Any]],
            Sequence[Mapping[str, Any]],
        ]
    ],
    species_by_protein: Mapping[int, int],
) -> Iterator[OrthologCandidate]:
    for component_id, nodes, memberships, events in sorted(
        components, key=lambda value: value[0]
    ):
        yield from generate_component_candidates(
            component_id=component_id,
            nodes=nodes,
            terminal_memberships=memberships,
            events=events,
            species_by_protein=species_by_protein,
        )
