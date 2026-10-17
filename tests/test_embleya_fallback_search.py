from dataclasses import replace

import pytest

from benchmarks.og_extraction.embleya_fallback_search import search_reserved_fallback
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig, SplitCandidate
from ogprofiler.hierarchy.soft_search import search_soft_binary_24


def candidate(gamma, count, violation=()):
    return SplitCandidate(
        gamma=gamma,
        membership=(0, 0, 1),
        child_count=count,
        quality=1,
        min_child_size=1,
        max_child_fraction=0.5,
        tiny_fragment_fraction=0,
        stability=1,
        adjusted_rand_index=1,
        normalized_mutual_info=1,
        inter_edge_fraction=0,
        intra_edge_fraction=1,
        valid=not violation,
        rejection_reason=violation[0] if violation else None,
        violations=violation,
        policy_valid=not violation,
        original_violations=violation,
        binary_eligible=not violation and count == 2,
        kway_eligible=not violation,
    )


def evaluator(gamma):
    if gamma < 0.03:
        return candidate(gamma, 1, ("NO_SPLIT",))
    if gamma < 0.08:
        return candidate(gamma, 2, ("MAX_CHILD_FRACTION",))
    if gamma < 1.5:
        return candidate(gamma, 4, ("UNSTABLE",))
    return candidate(gamma, 6)


def test_reserved_points_refine_fallback_without_exceeding_shared_budget():
    config = ResolutionSearchConfig(topology_policy="soft_binary_24_v2")
    calls = []

    def observe(gamma):
        calls.append(gamma)
        return evaluator(gamma)

    previous = search_soft_binary_24(config, evaluator)
    result = search_reserved_fallback(config, observe)
    assert previous.selected.gamma == 2.56
    assert result.selected.gamma == 1.6
    assert result.selected.phase == "FALLBACK_KWAY/FALLBACK_REFINE"
    assert result.selected.original_violations == ()
    assert len(calls) == len(set(calls)) == 24
    assert sum(c.phase.endswith("FALLBACK_REFINE") for c in result.candidates) == 3


def test_early_binary_unchanged():
    config = ResolutionSearchConfig(topology_policy="soft_binary_24_v2")

    def evaluate(gamma):
        return candidate(gamma, 2) if gamma >= 0.04 else candidate(gamma, 1, ("NO_SPLIT",))

    assert search_reserved_fallback(config, evaluate) == search_soft_binary_24(config, evaluate)


@pytest.mark.parametrize("cap", [5, 21, 24])
def test_all_rejected_preserves_failure_and_budget(cap):
    config = ResolutionSearchConfig(
        topology_policy="soft_binary_24_v2", max_candidate_evaluations=cap
    )

    def evaluate(gamma):
        return candidate(gamma, 2, ("UNSTABLE",))

    old = search_soft_binary_24(config, evaluate)
    new = search_reserved_fallback(config, evaluate)
    assert new == old and new.selected is None
    assert len(new.candidates) <= cap


def test_smaller_cap_does_not_enable_early_fallback():
    config = replace(
        ResolutionSearchConfig(topology_policy="soft_binary_24_v2"), max_candidate_evaluations=21
    )
    assert search_reserved_fallback(config, evaluator) == search_soft_binary_24(config, evaluator)


def test_experimental_adapter_obeys_real_robust_call_budget():
    import igraph as ig

    from benchmarks.og_extraction.embleya_subtree_repair import experimental_search
    from ogprofiler.hierarchy.leiden import LeidenCallCounter

    graph = ig.Graph.Full(6)
    graph.es["weight"] = [1.0] * graph.ecount()
    config = ResolutionSearchConfig(topology_policy="soft_binary_24_v2")
    counter = LeidenCallCounter(n_iterations=10)
    result = experimental_search(
        graph,
        config,
        method="rber",
        weights="weight",
        seed=42,
        counter=counter,
        stability_mode="robust",
    )
    assert counter.count == 3 * len(result.candidates)
    assert len(result.candidates) <= 24
    assert len({c.gamma for c in result.candidates}) == len(result.candidates)
    if result.selected:
        assert not result.selected.original_violations
