"""Species-overlap heuristics for binary and k-way hierarchy nodes."""

from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from typing import Any

from ogprofiler.exceptions import HierarchyError

NETWORK_EVENTS = {
    "SPECIATION_LIKE",
    "DUPLICATION_LIKE",
    "MIXED",
    "POLYTOMY",
    "AMBIGUOUS",
    "SPECIES_SPECIFIC",
}


@dataclass(frozen=True, slots=True)
class EventAnnotation:
    component_id: int
    cluster_id: int
    child_count: int
    network_event: str
    overlap_count: int
    overlap_score: float
    pairwise_overlap_summary: str
    confidence: float
    legacy_event: str | None


def _legacy_mapping(parent: int, children: list[int]) -> str:
    if len(children) > 2:
        return "III-3"
    intersection = children[0] & children[1]
    if intersection == parent:
        return "II"
    if intersection == 0:
        return "I"
    if children[0] == parent or children[1] == parent:
        return "III-1"
    return "III-2"


def _binary_event(score: float, threshold: float) -> tuple[str, float]:
    if score == 0:
        return "SPECIATION_LIKE", 1.0
    if threshold == 0 or score >= threshold:
        denominator = max(1.0 - threshold, 1e-12)
        return "DUPLICATION_LIKE", min(1.0, (score - threshold) / denominator)
    return "MIXED", max(0.0, 1.0 - score / threshold)


def annotate_network_events(
    nodes: list[dict[str, Any]],
    terminal_membership: list[dict[str, Any]],
    species_by_protein: dict[int, int],
    overlap_threshold: float,
) -> tuple[EventAnnotation, ...]:
    if not 0 <= overlap_threshold <= 1:
        raise HierarchyError("Network overlap threshold must be between 0 and 1")
    node_by_id = {int(row["cluster_id"]): row for row in nodes}
    children: dict[int, list[int]] = defaultdict(list)
    for row in nodes:
        parent = row["parent_id"]
        if parent is not None:
            children[int(parent)].append(int(row["cluster_id"]))
    for values in children.values():
        values.sort()

    bitmaps: dict[int, int] = defaultdict(int)
    for row in terminal_membership:
        protein_id = int(row["protein_id"])
        terminal_id = int(row["terminal_cluster_id"])
        bitmaps[terminal_id] |= 1 << species_by_protein[protein_id]
    for row in sorted(nodes, key=lambda item: int(item["depth"]), reverse=True):
        cluster_id = int(row["cluster_id"])
        for child_id in children[cluster_id]:
            bitmaps[cluster_id] |= bitmaps[child_id]
        if bitmaps[cluster_id].bit_count() != int(row["n_species"]):
            raise HierarchyError(f"Species bitmap mismatch for cluster {cluster_id}")

    annotations: list[EventAnnotation] = []
    for cluster_id in sorted(node_by_id):
        row = node_by_id[cluster_id]
        component_id = int(row["component_id"])
        child_ids = children[cluster_id]
        if "child_count" in row and int(row["child_count"]) != len(child_ids):
            raise HierarchyError(f"Child count mismatch for cluster {cluster_id}")
        child_bitmaps = [bitmaps[child_id] for child_id in child_ids]
        if not child_ids:
            event = "SPECIES_SPECIFIC" if bitmaps[cluster_id].bit_count() == 1 else "AMBIGUOUS"
            annotations.append(
                EventAnnotation(component_id, cluster_id, 0, event, 0, 0.0, "[]", 1.0, None)
            )
            continue

        pairs: list[dict[str, int | float]] = []
        scores: list[float] = []
        counts: list[int] = []
        for left_index, right_index in combinations(range(len(child_ids)), 2):
            left, right = child_bitmaps[left_index], child_bitmaps[right_index]
            overlap = (left & right).bit_count()
            denominator = min(left.bit_count(), right.bit_count())
            score = overlap / denominator if denominator else 0.0
            counts.append(overlap)
            scores.append(score)
            pairs.append(
                {
                    "left_cluster_id": child_ids[left_index],
                    "right_cluster_id": child_ids[right_index],
                    "overlap_count": overlap,
                    "normalized_overlap": score,
                }
            )
        maximum_count = max(counts)
        maximum_score = max(scores)
        if len(child_ids) == 2:
            event, confidence = _binary_event(scores[0], overlap_threshold)
        elif all(score == 0 for score in scores):
            event, confidence = "POLYTOMY", 1.0
        elif all(
            score > 0 and (overlap_threshold == 0 or score >= overlap_threshold)
            for score in scores
        ):
            event, confidence = "DUPLICATION_LIKE", min(scores)
        else:
            event, confidence = "MIXED", 1.0 - (max(scores) - min(scores))
        summary = json.dumps(pairs, sort_keys=True, separators=(",", ":"))
        annotations.append(
            EventAnnotation(
                component_id,
                cluster_id,
                len(child_ids),
                event,
                maximum_count,
                maximum_score,
                summary,
                max(0.0, min(1.0, confidence)),
                _legacy_mapping(bitmaps[cluster_id], child_bitmaps),
            )
        )
    return tuple(annotations)
