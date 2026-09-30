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
    symmetrization: str = "mean"


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


def lrb_cutoffs(
    hits: dict[tuple[int, int], NormalizedHit], tolerance: float = 1e-3
) -> dict[int, float]:
    """Reproduce OF GetMostDistant_s, including repeated-row assignment.

    In each species-pair RBH matrix, CSR nonzero traversal is ordered by
    query then target. Repeated query indexing keeps the last target's score,
    evaluated against the cutoff from previous species, not an intra-pair min.
    """
    best = best_hit_keys(hits, tolerance)
    external: dict[int, float] = {}
    last_rbh: dict[tuple[int, int], float] = {}
    for key in sorted(hits):
        item = hits[key]
        h = item.hit
        if h.query_species == h.target_species:
            continue
        external[h.query_id] = max(external.get(h.query_id, 0.0), item.normalized_score)
        if key in best and (key[1], key[0]) in best:
            last_rbh[h.query_id, h.target_species] = item.normalized_score
    cutoffs: dict[int, float] = {}
    for (query, _), score in last_rbh.items():
        cutoffs[query] = min(cutoffs.get(query, 1e9), score)
    queries = {item.hit.query_id for item in hits.values()}
    return {
        query: (external.get(query, 0.0) + 1e-6)
        if cutoffs.get(query, 1e9) > 1e8
        else cutoffs[query]
        for query in queries
    }


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
        thresholds = lrb_cutoffs(hits, config.best_hit_tolerance)
        selected = set()
        for key, item in hits.items():
            threshold = thresholds[key[0]]
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
    filtered = (
        filter_coverage(normalized_hits, config)
        if config.apply_coverage_filter
        else normalized_hits
    )
    complete = _deduplicate(filtered)
    selected = _select_directional(complete, config)
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
        if config.method == "lrb":
            # OF W=(C+C.T)*B uses COMPLETE B, not only selected directions.
            # One direction passed: factor 1; both passed: factor 2.
            factor = len(directions)
            forward = complete.get((u, v))
            reverse = complete.get((v, u))
            score_uv = factor * forward.normalized_score if forward is not None else 0.0
            score_vu = factor * reverse.normalized_score if reverse is not None else 0.0
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
