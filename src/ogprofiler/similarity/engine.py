"""Coverage, RBH/LRB retention, canonicalization, and symmetrization."""

from __future__ import annotations

import math
from dataclasses import dataclass

from ogprofiler.exceptions import EdgeConstructionError
from ogprofiler.similarity.models import NormalizedHit, RetainedEdge


@dataclass(frozen=True, slots=True)
class EdgeBuildConfig:
    method: str = "lrb"
    normalization: str = "legacy_nbs"
    nbs_fallback: str = "v1_zero"
    apply_coverage_filter: bool = False
    min_query_coverage: float = 0.0
    min_target_coverage: float = 0.0
    min_bidirectional_coverage: float = 0.0
    best_hit_tolerance: float = 1e-3
    symmetrization: str = "forward"


def filter_coverage(hits: list[NormalizedHit], config: EdgeBuildConfig) -> list[NormalizedHit]:
    return [
        item
        for item in hits
        if item.hit.query_coverage >= config.min_query_coverage
        and item.hit.target_coverage >= config.min_target_coverage
        and min(item.hit.query_coverage, item.hit.target_coverage)
        >= config.min_bidirectional_coverage
    ]


def _deduplicate(hits: list[NormalizedHit]) -> dict[tuple[int, int], NormalizedHit]:
    result: dict[tuple[int, int], NormalizedHit] = {}
    for item in hits:
        key = (item.hit.query_id, item.hit.target_id)
        previous = result.get(key)
        if previous is None or (item.normalized_score, -item.hit.evalue) > (
            previous.normalized_score,
            -previous.hit.evalue,
        ):
            result[key] = item
    return result


def best_hit_keys(
    hits: dict[tuple[int, int], NormalizedHit], tolerance: float
) -> set[tuple[int, int]]:
    per_species_max: dict[tuple[int, int], float] = {}
    external_max: dict[int, float] = {}
    for item in hits.values():
        hit = item.hit
        if hit.query_species != hit.target_species:
            group = (hit.query_id, hit.target_species)
            per_species_max[group] = max(
                per_species_max.get(group, -math.inf), item.normalized_score
            )
            external_max[hit.query_id] = max(
                external_max.get(hit.query_id, -math.inf), item.normalized_score
            )
    selected: set[tuple[int, int]] = set()
    for key, item in hits.items():
        hit = item.hit
        if hit.query_species == hit.target_species:
            threshold = external_max.get(hit.query_id, -1.0) - tolerance
        else:
            threshold = per_species_max[(hit.query_id, hit.target_species)] - tolerance
        if item.normalized_score > threshold:
            selected.add(key)
    return selected


def _select_directional(
    hits: dict[tuple[int, int], NormalizedHit], config: EdgeBuildConfig
) -> list[NormalizedHit]:
    best = best_hit_keys(hits, config.best_hit_tolerance)
    rbh = {
        key
        for key in best
        if (key[1], key[0]) in best and hits[key].hit.query_species != hits[key].hit.target_species
    }
    if config.method in {"rbh", "arb"}:
        selected = rbh
    elif config.method == "ar":
        selected = {key for key in hits if (key[1], key[0]) in hits}
    elif config.method == "lrb":
        rbh_scores: dict[int, list[float]] = {}
        external_best: dict[int, float] = {}
        for key, item in hits.items():
            if item.hit.query_species != item.hit.target_species:
                external_best[key[0]] = max(external_best.get(key[0], 0.0), item.normalized_score)
            if key in rbh:
                rbh_scores.setdefault(key[0], []).append(item.normalized_score)
        thresholds = {query: min(scores) for query, scores in rbh_scores.items()}
        selected = set()
        for key, item in hits.items():
            threshold = thresholds.get(key[0], external_best.get(key[0], 0.0) + 1e-6)
            if item.normalized_score >= threshold:
                selected.add(key)
    else:
        raise EdgeConstructionError(f"Unsupported edge method: {config.method}")
    return [hits[key] for key in sorted(selected)]


def _symmetrize(forward: float, reverse: float, method: str) -> float:
    if method == "forward":
        return forward if forward > 0 else reverse
    if method == "max":
        return max(forward, reverse)
    if method == "min":
        return min(forward, reverse)
    if method == "mean":
        return (forward + reverse) / 2
    if method == "geometric_mean":
        return math.sqrt(forward * reverse)
    raise EdgeConstructionError(f"Unsupported symmetrization: {method}")


def build_retained_edges(
    normalized_hits: list[NormalizedHit], config: EdgeBuildConfig
) -> tuple[list[NormalizedHit], list[RetainedEdge]]:
    filtered = filter_coverage(normalized_hits, config) if config.apply_coverage_filter else normalized_hits
    selected = _select_directional(_deduplicate(filtered), config)
    grouped: dict[tuple[int, int], list[NormalizedHit]] = {}
    for item in selected:
        u, v = sorted((item.hit.query_id, item.hit.target_id))
        grouped.setdefault((u, v), []).append(item)
    edges: list[RetainedEdge] = []
    for (u, v), directions in sorted(grouped.items()):
        score_uv = 0.0
        score_vu = 0.0
        species: dict[int, int] = {}
        coverages: list[float] = []
        for item in directions:
            hit = item.hit
            species[hit.query_id] = hit.query_species
            species[hit.target_id] = hit.target_species
            coverages.append(min(hit.query_coverage, hit.target_coverage))
            if hit.query_id == u:
                score_uv = max(score_uv, item.normalized_score)
            else:
                score_vu = max(score_vu, item.normalized_score)
        edge_type = config.method.upper()
        if species[u] == species[v]:
            edge_type = f"PARALOG_{edge_type}"
        edges.append(
            RetainedEdge(
                u=u,
                v=v,
                u_species=species[u],
                v_species=species[v],
                score_uv=score_uv,
                score_vu=score_vu,
                weight=_symmetrize(score_uv, score_vu, config.symmetrization),
                coverage=min(coverages),
                edge_type=edge_type,
            )
        )
    return selected, edges
