"""LCA-based gene/species-tree reconciliation evidence."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from ogprofiler.exceptions import PhylogenyError
from ogprofiler.phylogeny.newick import TreeNode, edge_rootings, leaf_names


@dataclass(frozen=True, slots=True)
class ReconciliationEvent:
    supporting_node: str
    phylo_event: str
    confidence: float
    species_node: str


@dataclass(frozen=True, slots=True)
class ReconciliationResult:
    events: tuple[ReconciliationEvent, ...]
    root_event: str
    confidence: float
    duplication_count: int


class ReconciliationBackend(Protocol):
    @property
    def name(self) -> str: ...

    def annotate(
        self,
        gene_tree: TreeNode,
        species_tree: TreeNode | None,
        leaf_species: dict[str, str],
    ) -> ReconciliationResult: ...


def _species_index(root: TreeNode) -> tuple[dict[str, int], dict[int, int | None], dict[int, int]]:
    leaf_to_id: dict[str, int] = {}
    parent: dict[int, int | None] = {}
    depth: dict[int, int] = {}
    next_id = 0

    def visit(node: TreeNode, parent_id: int | None, node_depth: int) -> int:
        nonlocal next_id
        node_id = next_id
        next_id += 1
        parent[node_id] = parent_id
        depth[node_id] = node_depth
        if node.is_leaf:
            if node.name is None:
                raise PhylogenyError("Species-tree leaf is unnamed")
            leaf_to_id[node.name] = node_id
        for child in node.children:
            visit(child, node_id, node_depth + 1)
        return node_id

    visit(root, None, 0)
    return leaf_to_id, parent, depth


def _lca(values: set[int], parent: dict[int, int | None], depth: dict[int, int]) -> int:
    if not values:
        raise PhylogenyError("LCA requires at least one species")
    ordered = sorted(values)
    current = ordered[0]
    for other in ordered[1:]:
        left, right = current, other
        while depth[left] > depth[right]:
            parent_id = parent[left]
            if parent_id is None:
                break
            left = parent_id
        while depth[right] > depth[left]:
            parent_id = parent[right]
            if parent_id is None:
                break
            right = parent_id
        while left != right:
            left_parent, right_parent = parent[left], parent[right]
            if left_parent is None or right_parent is None:
                raise PhylogenyError("Species-tree LCA traversal failed")
            left, right = left_parent, right_parent
        current = left
    return current


@dataclass(frozen=True, slots=True)
class LcaReconciliationBackend:
    name: str = "lca-reconciliation-v1"

    def annotate(
        self,
        gene_tree: TreeNode,
        species_tree: TreeNode | None,
        leaf_species: dict[str, str],
    ) -> ReconciliationResult:
        gene_leaves = set(leaf_names(gene_tree))
        if gene_leaves != set(leaf_species):
            raise PhylogenyError("Gene-tree leaves do not match the protein/species mapping")
        if species_tree is None:
            species_ids = {
                name: index
                for index, name in enumerate(sorted(set(leaf_species.values())))
            }
            parent: dict[int, int | None] = {value: None for value in species_ids.values()}
            depth = {value: 0 for value in species_ids.values()}
        else:
            species_ids, parent, depth = _species_index(species_tree)
            missing = sorted(set(leaf_species.values()) - set(species_ids))
            if missing:
                raise PhylogenyError(f"Species tree is missing selected species: {missing[0]}")

        events: list[ReconciliationEvent] = []
        internal_index = 0

        def visit(node: TreeNode) -> tuple[set[str], int]:
            nonlocal internal_index
            if node.is_leaf:
                if node.name is None:
                    raise PhylogenyError("Gene-tree leaf is unnamed")
                species = leaf_species[node.name]
                return {species}, species_ids[species]
            child_values = [visit(child) for child in node.children]
            species_sets = [value[0] for value in child_values]
            combined = set().union(*species_sets)
            if species_tree is None:
                mapped = min(species_ids[value] for value in combined)
            else:
                mapped = _lca({species_ids[value] for value in combined}, parent, depth)
            overlap = set()
            for left_index, left in enumerate(species_sets):
                for right in species_sets[left_index + 1 :]:
                    overlap.update(left & right)
            child_maps = [value[1] for value in child_values]
            duplication = bool(overlap) or (
                species_tree is not None and mapped in child_maps
            )
            event = "DUPLICATION" if duplication else "SPECIATION"
            confidence = 1.0 if not overlap else min(1.0, len(overlap) / len(combined))
            supporting = node.name or f"gene_node_{internal_index:06d}"
            internal_index += 1
            events.append(
                ReconciliationEvent(supporting, event, confidence, f"species_node_{mapped:06d}")
            )
            return combined, mapped

        visit(gene_tree)
        if not events:
            return ReconciliationResult((), "UNRESOLVED", 0.0, 0)
        root = events[-1]
        duplications = sum(event.phylo_event == "DUPLICATION" for event in events)
        return ReconciliationResult(tuple(events), root.phylo_event, root.confidence, duplications)


def species_tree_aware_root(
    gene_tree: TreeNode,
    species_tree: TreeNode,
    leaf_species: dict[str, str],
    backend: ReconciliationBackend,
) -> TreeNode:
    candidates = list(edge_rootings(gene_tree))
    if not candidates:
        return gene_tree
    scored = [
        (
            backend.annotate(candidate, species_tree, leaf_species).duplication_count,
            index,
            candidate,
        )
        for index, candidate in enumerate(candidates)
    ]
    return min(scored, key=lambda value: (value[0], value[1]))[2]
