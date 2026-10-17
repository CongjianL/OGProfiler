"""Pure V1-compatible OG selection with immutable component-local inputs."""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from collections.abc import Collection, Mapping, Sequence
from typing import Any, Literal

from ogprofiler.exceptions import HierarchyError
from ogprofiler.orthogroups.legacy_events import annotate_v1_events
from ogprofiler.orthogroups.models import (
    Orthogroup,
    OrthogroupResult,
    SelectionTrace,
    UnassignedProtein,
)

OG_EXTRACTION_ALGORITHM_VERSION = "v1-unrefined-active-view-v1"


class OrthogroupConflictError(HierarchyError):
    """Expose V1's overlapping diagnostic candidates; never publish them as OGs."""

    def __init__(self, result: OrthogroupResult, duplicate_members: tuple[int, ...]):
        self.result = result
        self.duplicate_members = duplicate_members
        super().__init__(
            f"V1-compatible OG selection overlaps in component {result.component_id}: "
            f"{len(duplicate_members)} duplicate proteins; inspect result.groups and result.trace"
        )


def extract_component_orthogroups(
    nodes: Sequence[Mapping[str, Any]],
    terminal_membership: Sequence[Mapping[str, Any]],
    species_by_protein: Mapping[int, int],
    original_by_protein: Mapping[int, str],
    *,
    total_species: int,
    overlap_count: int = 0,
    ssn_isolates: Collection[int] = (),
) -> OrthogroupResult:
    """V1 unrefined I/0/0 strategy; no graph mutation, I/O or process pools.

    total_species is the dataset count, not the component count. SSN isolates
    are explicit; a singleton hierarchy by itself does not imply SSN degree 0.
    Original node counts/events never change. Descendants are deactivated only
    at the end of a coverage level, exactly when V1 calls delete_vertices.
    """
    annotations = annotate_v1_events(nodes, terminal_membership, species_by_protein, overlap_count)
    if isinstance(total_species, bool) or not isinstance(total_species, int) or total_species < 1:
        raise HierarchyError("Dataset species count must be a positive integer")
    node_by_id = {int(row["cluster_id"]): row for row in nodes}
    component = annotations[0].component_id
    rank = {a.cluster_id: a.reference_order for a in annotations}
    events = {a.cluster_id: a for a in annotations}
    children: dict[int, list[int]] = defaultdict(list)
    neighbors: dict[int, list[int]] = defaultdict(list)
    root = -1
    for row in nodes:
        cluster = int(row["cluster_id"])
        if int(row["n_species"]) > total_species:
            raise HierarchyError("Component species count exceeds dataset species count")
        parent = row["parent_id"]
        if parent is None:
            root = cluster
        else:
            children[int(parent)].append(cluster)
            neighbors[int(parent)].append(cluster)
            neighbors[cluster].append(int(parent))
    for values in neighbors.values():
        values.sort(key=rank.__getitem__)
    terminal_proteins: dict[int, list[int]] = defaultdict(list)
    terminal_by_protein = {}
    identities = set()
    for row in terminal_membership:
        protein, terminal = int(row["protein_id"]), int(row["terminal_cluster_id"])
        original = original_by_protein.get(protein)
        if not isinstance(original, str) or not original:
            raise HierarchyError(f"Missing original ID for protein {protein}")
        identity = (species_by_protein[protein], original)
        if identity in identities:
            raise HierarchyError(f"Duplicate species/original protein identity: {identity}")
        identities.add(identity)
        terminal_by_protein[protein] = terminal
        terminal_proteins[terminal].append(protein)
    isolates = set(ssn_isolates)
    if not isolates <= terminal_by_protein.keys():
        raise HierarchyError("SSN isolate is outside this component's terminal membership")

    # A static DFS interval indexes original descendants; no per-node gene list.
    starts, ends = {}, {}
    all_proteins: list[int] = []
    stack = [(root, False)]
    while stack:
        cluster, exiting = stack.pop()
        if exiting:
            ends[cluster] = len(all_proteins)
            continue
        starts[cluster] = len(all_proteins)
        stack.append((cluster, True))
        if children[cluster]:
            stack.extend((child, False) for child in reversed(children[cluster]))
        else:
            all_proteins.extend(sorted(terminal_proteins[cluster]))

    active = set(node_by_id)
    groups: list[Orthogroup] = []
    trace: list[SelectionTrace] = []

    def append_group(cluster: int | None, proteins: list[int], level: int) -> None:
        members = tuple(sorted(proteins))
        member_keys = sorted((species_by_protein[p], original_by_protein[p]) for p in members)
        digest = hashlib.sha256(
            json.dumps(member_keys, ensure_ascii=False, separators=(",", ":")).encode()
        ).hexdigest()
        kind: Literal["EVENT_I", "RESIDUAL_NONE", "SSN_ISOLATE"] = (
            "SSN_ISOLATE" if level == 0 else ("RESIDUAL_NONE" if level == 1 else "EVENT_I")
        )
        groups.append(
            Orthogroup(
                component,
                len(groups),
                cluster,
                kind,
                events[cluster].v1_event if cluster is not None else None,
                level,
                members,
                len({species_by_protein[p] for p in members}),
                digest,
            )
        )
        trace.append(
            SelectionTrace(
                cluster,
                level,
                events[cluster].selection_event if cluster is not None else "None",
                "SELECTED",
                None,
            )
        )

    for level in range(total_species, 0, -1):
        candidates = [
            a.cluster_id
            for a in annotations
            if a.cluster_id in active
            and (
                a.selection_event == "None"
                if level == 1
                else a.selection_event == "I"
                and int(node_by_id[a.cluster_id]["n_species"]) == level
            )
        ]
        # Python's stable sort reproduces the original vertex-order tie.
        candidates.sort(key=lambda c: int(node_by_id[c]["n_genes"]), reverse=True)
        deleted: dict[int, int] = {}
        for candidate in candidates:
            if candidate in deleted:
                trace.append(
                    SelectionTrace(
                        candidate,
                        level,
                        events[candidate].selection_event,
                        "SKIPPED_CONSUMED",
                        deleted[candidate],
                    )
                )
                continue
            if level == 1:
                append_group(candidate, all_proteins[starts[candidate] : ends[candidate]], level)
                continue
            collected = []
            frontier = [candidate]
            while frontier:
                next_frontier = []
                for cluster in frontier:
                    adjacent = [other for other in neighbors[cluster] if other in active]
                    # All events have undergone None->'None' normalization.
                    if len(adjacent) > 1:
                        for other in adjacent:
                            if int(node_by_id[other]["n_genes"]) < int(
                                node_by_id[cluster]["n_genes"]
                            ):
                                next_frontier.append(other)
                                deleted[other] = candidate
                                trace.append(
                                    SelectionTrace(
                                        other,
                                        level,
                                        events[other].selection_event,
                                        "DESCENDANT_CONSUMED",
                                        candidate,
                                    )
                                )
                    else:
                        collected.extend(all_proteins[starts[cluster] : ends[cluster]])
                frontier = next_frontier
            append_group(candidate, collected, level)
        if level > 1:
            active.difference_update(deleted)
    for protein in sorted(isolates):
        append_group(None, [protein], 0)
    occurrences = Counter(p for group in groups for p in group.protein_ids)
    duplicates = tuple(sorted(p for p, count in occurrences.items() if count > 1))
    unassigned = tuple(
        UnassignedProtein(p, terminal_by_protein[p], "NO_V1_SELECTION")
        for p in sorted(terminal_by_protein)
        if p not in occurrences
    )
    result = OrthogroupResult(
        component,
        tuple(groups),
        tuple(trace),
        unassigned,
        tuple(a.cluster_id for a in annotations if a.cluster_id in active),
    )
    if duplicates:
        raise OrthogroupConflictError(result, duplicates)
    return result
