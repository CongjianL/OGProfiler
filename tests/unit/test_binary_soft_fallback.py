"""Explicit topology-policy experiment, not implicit terminalization."""

from dataclasses import replace

from benchmarks.og_extraction.binary_coverage_schedule import (
    SOFT,
    UPPER_GUARD,
    coverage_search,
)
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig
from tests.unit.test_binary_recursive_audit import candidate


def test_guard_supplies_upper_side_after_invalid_binary_at_lower_endpoint():
    result = coverage_search(
        ResolutionSearchConfig(),
        lambda g: candidate(g, 2 if g < 0.04 else 3, g >= 0.04, ("UNSTABLE",) if g < 0.04 else ()),
        protocol=UPPER_GUARD,
    )
    assert [c.gamma for c in result.candidates if c.phase == "UPPER_GUARD"] == [0.02, 0.04]
    assert result.selected is None
    assert len(result.candidates) <= 24


def test_soft_selects_only_measured_original_gated_kway_without_extra_calls():
    calls = []

    def evaluate(gamma):
        calls.append(gamma)
        return candidate(gamma, 3)

    result = coverage_search(ResolutionSearchConfig(), evaluate, protocol=SOFT)
    assert len(calls) == len(result.candidates) == 24
    assert result.selected.gamma == 0.01
    assert result.selected.child_count == 3
    assert result.selected.phase.startswith("FALLBACK_KWAY/")
    assert result.selected.violations == ()
    assert result.search_status == "ACCEPTED"


def test_soft_does_not_relax_original_gates_or_fallback_on_truncated_search():
    result = coverage_search(
        ResolutionSearchConfig(), lambda g: candidate(g, 3, False, ("UNSTABLE",)), protocol=SOFT
    )
    assert result.selected is None
    assert result.search_status == "REJECTED_ALL_TESTED"
    shortened = coverage_search(
        replace(ResolutionSearchConfig(), max_candidate_evaluations=2),
        lambda g: candidate(g, 3),
        protocol=SOFT,
    )
    assert shortened.selected is None
    assert shortened.search_status == "EVALUATION_BUDGET_EXHAUSTED"


def test_soft_keeps_binary_preference_even_when_lower_gamma_kway_passes():
    result = coverage_search(
        ResolutionSearchConfig(), lambda g: candidate(g, 3 if g < 0.06 else 2), protocol=SOFT
    )
    assert result.selected.child_count == 2
    assert not result.selected.phase.startswith("FALLBACK_KWAY/")


def test_expanded_guard_covers_rejected_nonbinary_split_with_same_24_point_cap():
    from benchmarks.og_extraction.binary_coverage_schedule import SOFT_V2

    calls = []

    def evaluate(gamma):
        calls.append(gamma)
        return candidate(gamma, 3, gamma >= 0.4, ("MAX_CHILD_FRACTION",) if gamma < 0.4 else ())

    result = coverage_search(ResolutionSearchConfig(), evaluate, protocol=SOFT_V2)
    assert result.selected.child_count == 3
    assert result.selected.phase.startswith("FALLBACK_KWAY/")
    assert len(calls) == len(result.candidates) == 24
    assert sum("UPPER_GUARD" in c.phase for c in result.candidates) > 3
    assert result.selected.violations == ()
