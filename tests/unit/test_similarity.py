from __future__ import annotations

import math

import pytest

from ogprofiler.similarity.engine import (
    EdgeBuildConfig,
    best_hit_keys,
    build_retained_edges,
    filter_coverage,
)
from ogprofiler.similarity.models import DirectionalHit, NormalizedHit
from ogprofiler.similarity.normalization import legacy_nbs, retain_top_data


def hit(
    query: int,
    target: int,
    query_species: int,
    target_species: int,
    score: float,
    coverage: float = 100.0,
) -> NormalizedHit:
    return NormalizedHit(
        DirectionalHit(
            query,
            target,
            query_species,
            target_species,
            score,
            50.0,
            coverage,
            coverage,
            1e-10,
        ),
        score,
    )


def test_legacy_nbs_reproduces_log_length_fit_and_small_sample_fallback() -> None:
    raw = [
        DirectionalHit(0, 1, 0, 1, 40.0, 90.0, 100.0, 100.0, 1e-20),
        DirectionalHit(0, 2, 0, 1, 80.0, 90.0, 100.0, 100.0, 1e-20),
        DirectionalHit(3, 0, 2, 0, 7.0, 90.0, 100.0, 100.0, 1e-10),
    ]
    # The first group follows bit_score = 2 * sqrt(length_product).
    lengths = {0: 100, 1: 4, 2: 16, 3: 5}
    normalized = legacy_nbs(raw, lengths)
    assert [item.normalized_score for item in normalized[:2]] == pytest.approx([1.0, 1.0])
    assert normalized[2].normalized_score == pytest.approx(1.0)


def test_v1_top_bin_logic_keeps_overlapping_95th_percentile_members() -> None:
    lengths = [float(value) for value in range(1, 101)]
    scores = [float(value) for value in range(1, 101)]
    top_lengths, top_scores = retain_top_data(lengths, scores)
    assert top_lengths == [20.0, 21.0, 40.0, 41.0, 60.0, 61.0, 80.0, 81.0, 100.0]
    assert top_scores == top_lengths


def test_coverage_filter_applies_directional_and_minimum_thresholds() -> None:
    values = [hit(0, 1, 0, 1, 1.0, 80.0), hit(0, 2, 0, 1, 1.0, 49.0)]
    config = EdgeBuildConfig(
        min_query_coverage=50,
        min_target_coverage=50,
        min_bidirectional_coverage=75,
    )
    assert filter_coverage(values, config) == [values[0]]


def test_best_hits_rbh_and_same_species_paralog_handling() -> None:
    values = [
        hit(0, 2, 0, 1, 10.0),
        hit(2, 0, 1, 0, 9.0),
        hit(0, 3, 0, 1, 9.998),
        hit(3, 0, 1, 0, 8.0),
        hit(0, 1, 0, 0, 10.0),
    ]
    keyed = {(item.hit.query_id, item.hit.target_id): item for item in values}
    best = best_hit_keys(keyed, tolerance=1e-3)
    assert (0, 2) in best and (2, 0) in best
    assert (0, 3) not in best
    assert (0, 1) in best

    selected, edges = build_retained_edges(values, EdgeBuildConfig(method="rbh"))
    assert {(item.hit.query_id, item.hit.target_id) for item in selected} == {(0, 2), (2, 0)}
    assert len(edges) == 1
    assert edges[0].edge_type == "RBH"


def test_lrb_threshold_paralogs_reverse_only_and_no_rbh_fallback() -> None:
    values = [
        hit(0, 2, 0, 1, 8.0),
        hit(2, 0, 1, 0, 8.0),
        hit(0, 3, 0, 2, 9.0),  # retained above the most-distant RBH threshold
        hit(0, 1, 0, 0, 8.5),  # same-species paralog retained by the same threshold
        hit(4, 0, 2, 0, 5.0),  # no RBH: legacy fallback is best + epsilon
        hit(0, 4, 0, 2, 4.0),
        hit(0, 5, 0, 2, 10.0),
        hit(5, 0, 2, 0, 10.0),
    ]
    selected, edges = build_retained_edges(values, EdgeBuildConfig(method="lrb"))
    keys = {(item.hit.query_id, item.hit.target_id) for item in selected}
    assert (0, 3) in keys
    assert (0, 1) in keys
    assert (4, 0) not in keys
    edge_by_pair = {(edge.u, edge.v): edge for edge in edges}
    assert edge_by_pair[(0, 3)].score_uv == 9.0
    assert edge_by_pair[(0, 3)].score_vu == 0.0
    assert edge_by_pair[(0, 1)].edge_type == "PARALOG_LRB"


@pytest.mark.parametrize(
    ("method", "expected"),
    [("max", 9.0), ("min", 4.0), ("mean", 6.5), ("geometric_mean", 6.0)],
)
def test_symmetrization_is_order_invariant(method: str, expected: float) -> None:
    values = [hit(0, 1, 0, 1, 9.0), hit(1, 0, 1, 0, 4.0)]
    config = EdgeBuildConfig(method="ar", symmetrization=method)
    first = build_retained_edges(values, config)[1]
    second = build_retained_edges(list(reversed(values)), config)[1]
    assert first == second
    assert first[0].weight == pytest.approx(expected)
    assert math.isfinite(first[0].weight)
