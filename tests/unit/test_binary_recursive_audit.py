"""Observable finite experimental search contract, independent candidate evaluator."""

from dataclasses import replace

from benchmarks.og_extraction.binary_recursive_audit import finite_binary_search
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig, SplitCandidate


def candidate(gamma, count, valid=True, violations=()):
    return SplitCandidate(
        gamma=gamma,
        membership=(0, 0, 1),
        child_count=count,
        quality=1,
        min_child_size=1,
        max_child_fraction=2 / 3,
        tiny_fragment_fraction=1,
        stability=1,
        adjusted_rand_index=1,
        normalized_mutual_info=1,
        inter_edge_fraction=0,
        intra_edge_fraction=1,
        valid=valid,
        rejection_reason=violations[0] if violations else None,
        violations=violations,
        policy_valid=valid,
    )


def test_shared_schedule_finds_binary_between_coarse_points_and_refines():
    def evaluate(gamma):
        return candidate(gamma, 1 if gamma < 0.79 else 2 if gamma <= 0.81 else 3)

    result = finite_binary_search(ResolutionSearchConfig(), evaluate)
    assert result.search_status == "ACCEPTED"
    assert result.selected.gamma == 0.8
    assert [c.gamma for c in result.candidates if c.phase == "TARGET_PROBE"] == [0.96, 0.8]
    assert len(result.candidates) == 16
    assert len({c.gamma for c in result.candidates}) == 16


def test_binary_keeps_original_rejection_and_complete_schedule_is_not_exhaustion():
    result = finite_binary_search(
        ResolutionSearchConfig(), lambda g: candidate(g, 2, False, ("UNSTABLE",))
    )
    assert result.selected is None
    assert result.search_status == result.terminal_reason == "REJECTED_ALL_TESTED"
    assert len(result.candidates) <= 21
    assert all(c.violations == ("UNSTABLE",) for c in result.candidates)


def test_budget_blocks_requested_point_without_kway_fallback():
    cfg = replace(ResolutionSearchConfig(), max_candidate_evaluations=2)
    result = finite_binary_search(cfg, lambda g: candidate(g, 3))
    assert len(result.candidates) == 2
    assert result.search_status == "EVALUATION_BUDGET_EXHAUSTED"
    assert result.selected is None
    assert all(c.violations == ("TARGET_CHILD_COUNT",) for c in result.candidates)


def test_target_policy_boundaries_are_independent_of_stop_size():
    from benchmarks.og_extraction.binary_recursive_audit import uses_binary_target

    assert not uses_binary_target(2)
    assert uses_binary_target(3)
    assert uses_binary_target(9999)
    assert not uses_binary_target(10000)


def test_failed_full_recursion_retains_members_as_unresolved_not_terminal():
    from unittest.mock import patch

    import igraph as ig

    from benchmarks.og_extraction.binary_recursive_audit import experimental_search
    from ogprofiler.graph.components import Component
    from ogprofiler.hierarchy.engine import HierarchyConfig, infer_component_hierarchy
    from ogprofiler.hierarchy.validation import validate_hierarchy

    component = Component(0, (0, 1, 2), ())
    graph = ig.Graph.Full(3)
    graph.es["weight"] = [1.0, 1.0, 1.0]
    with patch("ogprofiler.hierarchy.engine.search_resolution", experimental_search):
        result = infer_component_hierarchy(component, HierarchyConfig(), {0: 0, 1: 1, 2: 2}, graph)
    validate_hierarchy(component, result)
    assert result.nodes[0].split_status == "UNRESOLVED"
    assert result.nodes[0].terminal_reason == "REJECTED_ALL_TESTED"
    assert sorted(result.terminal_membership) == [(0, 0), (1, 0), (2, 0)]
    assert result.metrics.leiden_calls <= 72


def test_all_phases_share_exactly_24_points_and_full_cap_can_be_accepted():
    import math

    calls = []

    def evaluate(gamma):
        calls.append(gamma)
        count = (
            2 if math.isclose(gamma, math.sqrt(10), rel_tol=1e-12) else (1 if gamma < 0.7 else 3)
        )
        return candidate(gamma, count)

    result = finite_binary_search(ResolutionSearchConfig(), evaluate)
    assert result.search_status == "ACCEPTED"
    assert len(calls) == len(set(calls)) == 24
    assert math.isclose(result.selected.gamma, math.sqrt(10), rel_tol=1e-12)
