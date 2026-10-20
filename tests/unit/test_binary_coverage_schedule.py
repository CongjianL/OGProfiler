"""Predeclared public evaluator seam for adaptive coverage arms."""

from dataclasses import replace

import pytest

from benchmarks.og_extraction.binary_coverage_schedule import FAIR, STABILITY, coverage_search
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig
from tests.unit.test_binary_recursive_audit import candidate


@pytest.mark.parametrize("protocol", [FAIR, STABILITY])
def test_count_bracket_releases_endpoint_and_rescue_quota_without_new_budget(protocol):
    result = coverage_search(
        ResolutionSearchConfig(),
        lambda g: candidate(g, 1 if g < 0.79 else 2 if g <= 0.81 else 3),
        protocol=protocol,
    )
    assert result.search_status == "ACCEPTED"
    assert result.selected.gamma == 0.8
    assert len(result.candidates) <= 24
    assert not any(c.phase in ("ENDPOINT", "RESCUE") for c in result.candidates)


@pytest.mark.parametrize("protocol", [FAIR, STABILITY])
def test_unstable_binary_is_never_published_and_original_gates_remain(protocol):
    result = coverage_search(
        ResolutionSearchConfig(),
        lambda g: candidate(g, 1 if g < 0.6 else 2 if g <= 0.65 else 3, False, ("UNSTABLE",)),
        protocol=protocol,
    )
    assert result.selected is None
    assert result.search_status == "REJECTED_ALL_TESTED"
    assert len(result.candidates) == 24
    assert all("UNSTABLE" in c.violations for c in result.candidates)
    assert sum(c.phase.startswith("BORROWED_") for c in result.candidates) == 3


def test_stability_arm_actually_probes_internal_binary_region():
    def evaluate(gamma):
        return candidate(gamma, 1 if gamma < 0.6 else 2 if gamma < 1 else 3, False, ("UNSTABLE",))

    a = coverage_search(ResolutionSearchConfig(), evaluate, protocol=FAIR)
    b = coverage_search(ResolutionSearchConfig(), evaluate, protocol=STABILITY)
    assert any("STABILITY_BAND" in c.phase for c in b.candidates)
    assert not any("STABILITY_BAND" in c.phase for c in a.candidates)
    assert len(b.candidates) == len(a.candidates) == 24


@pytest.mark.parametrize("protocol", [FAIR, STABILITY])
def test_explicit_smaller_cap_blocks_search_and_unknown_protocol_has_no_implicit_default(protocol):
    cfg = replace(ResolutionSearchConfig(), max_candidate_evaluations=12)
    result = coverage_search(cfg, lambda g: candidate(g, 3), protocol=protocol)
    assert result.search_status == "EVALUATION_BUDGET_EXHAUSTED"
    assert len(result.candidates) == 12
    with pytest.raises(ValueError, match="Unknown"):
        coverage_search(cfg, lambda g: candidate(g, 3), protocol="unknown")
