"""Component-local V1 event compatibility on immutable hierarchy inputs."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from typing import Any

from ogprofiler.exceptions import HierarchyError
from ogprofiler.orthogroups.models import V1Event, V1EventAnnotation

V1_EVENT_ALGORITHM_VERSION = "v1-og-degree-integer-overlap-v1"


def validate_overlap_count(overlap_count: int) -> None:
    if isinstance(overlap_count, bool) or not isinstance(overlap_count, int) or overlap_count < 0:
        raise HierarchyError("V1 species overlap count must be a nonnegative integer")


def v1_event_for_node(
    parent_bitmap: int,
    smaller_neighbor_bitmaps: Sequence[int],
    n_genes: int,
    undirected_degree: int,
    overlap_count: int = 0,
) -> V1Event | None:
    """Exact V1 ordering; degrees 0/1 differ from parent/children counts.

    Eligible neighbors have STRICTLY smaller original gene counts. Invalid
    degree>1 shapes are reported, not silently reinterpreted as binary nodes.
    """
    validate_overlap_count(overlap_count)
    if parent_bitmap <= 0 or n_genes <= 0 or undirected_degree < 0:
        raise HierarchyError("V1 event view requires positive members/species and valid degree")
    if undirected_degree == 0:
        return "I" if n_genes == parent_bitmap.bit_count() else None
    if undirected_degree == 1:
        return None
    if len(smaller_neighbor_bitmaps) > 2:
        return "III-3"
    if len(smaller_neighbor_bitmaps) < 2:
        raise HierarchyError(
            "V1 GetEvolutionEvents raises IndexError for degree>1 with fewer than two "
            "strictly-smaller neighbors; unsupported legacy shape"
        )
    left, right = smaller_neighbor_bitmaps
    intersection = left & right
    if intersection == parent_bitmap:
        return "II"
    if intersection == 0:
        return "I"
    if (
        intersection.bit_count() <= overlap_count
        and intersection.bit_count() < parent_bitmap.bit_count()
    ):
        return "I"
    if left == parent_bitmap or right == parent_bitmap:
        return "III-1"
    return "III-2"


def annotate_v1_events(
    nodes: Sequence[Mapping[str, Any]],
    terminal_membership: Sequence[Mapping[str, Any]],
    species_by_protein: Mapping[int, int],
    overlap_count: int = 0,
) -> tuple[V1EventAnnotation, ...]:
    """Validate and annotate one component, using O(nodes+members) storage.

    Raw unset Event is Python None. selection_event is separately normalized
    to the string 'None', matching V1's unrefined pipeline. Input row order is
    retained as a traceable reference order for later stable candidate ties.
    """
    validate_overlap_count(overlap_count)
    by_id = {int(row["cluster_id"]): row for row in nodes}
    if not nodes or len(by_id) != len(nodes):
        raise HierarchyError("V1 event view requires nonempty unique hierarchy nodes")
    components = {int(row["component_id"]) for row in nodes}
    if len(components) != 1:
        raise HierarchyError("V1 event annotation is component-local")
    children: dict[int, list[int]] = defaultdict(list)
    roots = []
    for row in nodes:
        cluster = int(row["cluster_id"])
        parent = row["parent_id"]
        if parent is None:
            roots.append(cluster)
        elif int(parent) not in by_id or int(parent) == cluster:
            raise HierarchyError(f"Invalid parent for V1 event cluster {cluster}")
        else:
            children[int(parent)].append(cluster)
    if len(roots) != 1:
        raise HierarchyError("V1 event component must have exactly one root")
    order = []
    stack = [(roots[0], 0)]
    seen = set()
    while stack:
        cluster, depth = stack.pop()
        if cluster in seen:
            raise HierarchyError("Cycle in V1 event hierarchy")
        seen.add(cluster)
        if "depth" in by_id[cluster] and int(by_id[cluster]["depth"]) != depth:
            raise HierarchyError(f"Depth mismatch for V1 event cluster {cluster}")
        if "child_count" in by_id[cluster] and int(by_id[cluster]["child_count"]) != len(
            children[cluster]
        ):
            raise HierarchyError(f"Child count mismatch for V1 event cluster {cluster}")
        order.append(cluster)
        stack.extend((child, depth + 1) for child in reversed(children[cluster]))
    if len(seen) != len(nodes):
        raise HierarchyError("Disconnected or cyclic V1 event component")
    bitmaps: dict[int, int] = defaultdict(int)
    counts: dict[int, int] = defaultdict(int)
    proteins = set()
    for row in terminal_membership:
        protein, terminal = int(row["protein_id"]), int(row["terminal_cluster_id"])
        if terminal not in by_id or children[terminal] or protein in proteins:
            raise HierarchyError("Invalid or duplicate terminal membership in V1 event view")
        species = species_by_protein.get(protein)
        if species is None or species < 0:
            raise HierarchyError(f"Missing/invalid species for protein {protein}")
        proteins.add(protein)
        bitmaps[terminal] |= 1 << species
        counts[terminal] += 1
    for cluster in reversed(order):
        for child in children[cluster]:
            bitmaps[cluster] |= bitmaps[child]
            counts[cluster] += counts[child]
        row = by_id[cluster]
        if counts[cluster] != int(row["n_genes"]) or bitmaps[cluster].bit_count() != int(
            row["n_species"]
        ):
            raise HierarchyError(f"Member/species count mismatch for V1 event cluster {cluster}")
    annotations = []
    for index, row in enumerate(nodes):
        cluster = int(row["cluster_id"])
        degree = len(children[cluster]) + int(row["parent_id"] is not None)
        eligible = [bitmaps[c] for c in children[cluster] if counts[c] < counts[cluster]]
        event = v1_event_for_node(
            bitmaps[cluster], eligible, counts[cluster], degree, overlap_count
        )
        annotations.append(
            V1EventAnnotation(
                int(row["component_id"]), cluster, index, degree, len(eligible), event
            )
        )
    return tuple(annotations)
