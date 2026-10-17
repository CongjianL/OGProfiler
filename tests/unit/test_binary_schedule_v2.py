"""Shared-budget experimental schedule behavior at its evaluator seam."""

from dataclasses import replace

from benchmarks.og_extraction.binary_schedule_v2 import finite_binary_search_v2
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig
from tests.unit.test_binary_recursive_audit import candidate


def test_rejected_binary_probes_both_boundaries_without_relaxing_gates():
    def evaluate(gamma):
        count = 1 if gamma < 0.6 else 2 if gamma <= 0.65 else 3
        return candidate(gamma, count, False, ("UNSTABLE",))

    result = finite_binary_search_v2(ResolutionSearchConfig(), evaluate)
    probes = [c.gamma for c in result.candidates if c.phase == "TARGET_PROBE"]
    assert probes[:2] == [0.48, 0.96]
    assert result.selected is None
    assert len(result.candidates) == 24
    assert result.search_status == "REJECTED_ALL_TESTED"
    assert all("UNSTABLE" in c.violations for c in result.candidates)


def test_failure_borrows_refinement_slots_but_never_exceeds_24():
    result = finite_binary_search_v2(ResolutionSearchConfig(), lambda g: candidate(g, 3))
    assert len(result.candidates) == 24
    assert sum(c.phase == "BORROWED_REFINE" for c in result.candidates) == 3
    assert result.selected is None
    assert result.search_status == "REJECTED_ALL_TESTED"


def test_lower_budget_truncation_is_exhausted_not_complete_rejection():
    result = finite_binary_search_v2(
        replace(ResolutionSearchConfig(), max_candidate_evaluations=12), lambda g: candidate(g, 3)
    )
    assert len(result.candidates) == 12
    assert result.search_status == "EVALUATION_BUDGET_EXHAUSTED"


def test_qualifying_binary_still_requires_original_fraction_gate():
    result = finite_binary_search_v2(
        ResolutionSearchConfig(), lambda g: candidate(g, 2, False, ("MAX_CHILD_FRACTION",))
    )
    assert result.selected is None
    assert all(c.violations == ("MAX_CHILD_FRACTION",) for c in result.candidates)


def test_deeper_target_search_recovers_valid_binary_missed_by_five_probe_schedule():
    def evaluate(gamma):
        count = 1 if gamma < 0.969 else 2 if gamma <= 0.977 else 3
        return candidate(gamma, count)

    result = finite_binary_search_v2(ResolutionSearchConfig(), evaluate)
    assert result.search_status == "ACCEPTED"
    assert result.selected.gamma == 0.97
    assert len(result.candidates) <= 24
    assert sum(c.phase == "REFINE" for c in result.candidates) == 3
