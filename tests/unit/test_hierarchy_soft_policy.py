"""Public default topology configuration and hierarchy search contract."""

from dataclasses import replace

import igraph as ig
import pytest

from ogprofiler.config import hierarchy_config, load_config
from ogprofiler.exceptions import InputError
from ogprofiler.hierarchy.leiden import LeidenCallCounter
from ogprofiler.hierarchy.resolution import search_resolution


def test_soft_policy_is_default_and_incompatible_with_legacy():
    default = hierarchy_config(load_config()["hierarchy"])
    assert default.resolution.topology_policy == "soft_binary_24_v2"
    assert default.max_depth == 42
    soft = hierarchy_config(
        load_config(overrides=["hierarchy.topology_policy=soft_binary_24_v2"])["hierarchy"]
    )
    assert soft.resolution.topology_policy == "soft_binary_24_v2"
    assert (soft.recursion_stop_size, soft.leiden_iterations) == (1, 10)
    assert soft.resolution.stability_threshold == 0.9
    with pytest.raises(InputError):
        load_config(
            overrides=[
                "hierarchy.topology_policy=soft_binary_24_v2",
                "hierarchy.admission_policy=legacy_strict",
                "hierarchy.resolution_strategy=adaptive",
            ]
        )


def test_soft_search_publishes_explicit_original_gated_fallback_in_shared_budget():
    cfg = hierarchy_config(
        load_config(overrides=["hierarchy.topology_policy=soft_binary_24_v2"])["hierarchy"]
    )
    graph = ig.Graph.Full(3)
    graph.es["weight"] = [1.0, 1.0, 1.0]
    counter = LeidenCallCounter(n_iterations=10)
    result = search_resolution(
        graph,
        cfg.resolution,
        method="rber",
        weights="weight",
        seed=42,
        counter=counter,
        stability_mode="robust",
    )
    assert result.selected is not None
    assert result.selected.selection_kind == "FALLBACK_KWAY"
    assert result.selected.kway_eligible and not result.selected.binary_eligible
    assert not result.selected.original_violations
    assert result.selected.phase.startswith("FALLBACK_KWAY/")
    assert counter.count == len(result.candidates) * 3 <= 72
    assert len({c.gamma for c in result.candidates}) == len(result.candidates)
    shortened = replace(cfg.resolution, max_candidate_evaluations=1)
    failed = search_resolution(
        graph,
        shortened,
        method="rber",
        weights="weight",
        seed=42,
        counter=LeidenCallCounter(n_iterations=10),
        stability_mode="robust",
    )
    assert failed.selected is None
    assert failed.search_status == "EVALUATION_BUDGET_EXHAUSTED"


def _candidate(gamma, count, violations=()):
    from ogprofiler.hierarchy.resolution import SplitCandidate

    valid = not violations
    return SplitCandidate(
        gamma=gamma,
        membership=tuple(range(count)),
        child_count=count,
        quality=1.0,
        min_child_size=1,
        max_child_fraction=1 / count,
        tiny_fragment_fraction=0.0,
        stability=1.0,
        adjusted_rand_index=1.0,
        normalized_mutual_info=1.0,
        inter_edge_fraction=0.0,
        intra_edge_fraction=1.0,
        valid=valid,
        rejection_reason=violations[0] if violations else None,
        violations=violations,
        structural_valid=True,
        policy_valid=valid,
        binary_eligible=valid and count == 2,
        kway_eligible=valid,
        original_violations=violations,
    )


@pytest.mark.parametrize("violation", ["UNSTABLE", "MAX_CHILD_FRACTION"])
def test_soft_full_budget_rejection_preserves_original_gate(violation):
    from ogprofiler.hierarchy.soft_search import search_soft_binary_24

    config = hierarchy_config(
        load_config(overrides=["hierarchy.topology_policy=soft_binary_24_v2"])["hierarchy"]
    ).resolution
    result = search_soft_binary_24(config, lambda g: _candidate(g, 3, (violation,)))
    assert len(result.candidates) == 24
    assert result.search_status == "REJECTED_ALL_TESTED"
    assert result.selected is None
    assert all(c.original_violations == (violation,) for c in result.candidates)
    assert not any(c.kway_eligible for c in result.candidates)


def test_soft_binary_priority_and_explicit_refinement_truncation():
    from ogprofiler.hierarchy.soft_search import search_soft_binary_24

    config = hierarchy_config(
        load_config(overrides=["hierarchy.topology_policy=soft_binary_24_v2"])["hierarchy"]
    ).resolution
    result = search_soft_binary_24(config, lambda g: _candidate(g, 3 if g < 0.06 else 2))
    assert result.selected.selection_kind == "BINARY"
    assert result.selected.gamma >= 0.06
    assert len(result.candidates) <= 24
    truncated = search_soft_binary_24(
        replace(config, max_candidate_evaluations=1), lambda g: _candidate(g, 2)
    )
    assert truncated.search_status == "ACCEPTED"
    assert truncated.selected.refinement_truncated


def test_soft_rejects_oversized_budget_and_keeps_two_protein_policy():
    with pytest.raises(InputError):
        load_config(
            overrides=[
                "hierarchy.topology_policy=soft_binary_24_v2",
                "hierarchy.max_candidate_evaluations=25",
            ]
        )
    cfg = hierarchy_config(
        load_config(overrides=["hierarchy.topology_policy=soft_binary_24_v2"])["hierarchy"]
    )
    graph = ig.Graph.Full(2)
    graph.es["weight"] = [1.0]
    result = search_resolution(
        graph,
        cfg.resolution,
        method="rber",
        weights="weight",
        seed=42,
        counter=LeidenCallCounter(n_iterations=10),
        stability_mode="robust",
    )
    assert result.selected.selection_kind == "KWAY"


def test_h4_migration_changes_only_explicit_topology_opt_in():
    from benchmarks.og_extraction.hierarchy_regression import migrate_config

    previous = load_config()
    default = migrate_config(previous)
    soft = migrate_config(previous, topology_policy="soft_binary_24_v2")
    assert default["hierarchy"]["topology_policy"] == "kway_v1"
    assert soft["hierarchy"]["topology_policy"] == "soft_binary_24_v2"
    soft["hierarchy"]["topology_policy"] = "kway_v1"
    assert soft == default
