"""Lower-bound experiments retain the same public finite-search contract."""

from dataclasses import replace

from benchmarks.og_extraction.binary_coverage_schedule import FAIR, coverage_search
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig
from tests.unit.test_binary_recursive_audit import candidate


def test_lower_bound_admits_low_gamma_only_in_explicit_experimental_configuration():
    frozen = ResolutionSearchConfig()

    def evaluate(gamma):
        return candidate(gamma, 1 if gamma < 0.0004 else 2 if gamma <= 0.002 else 3)

    control = coverage_search(frozen, evaluate, protocol=FAIR)
    assert control.selected is None
    for lower in (0.001, 0.0001):
        cfg = replace(frozen, gamma_min=lower)
        result = coverage_search(cfg, evaluate, protocol=FAIR)
        assert result.search_status == "ACCEPTED"
        assert result.selected.gamma < frozen.gamma_min
        assert all(lower <= c.gamma <= frozen.gamma_max for c in result.candidates)
        assert len(result.candidates) <= 24
    assert frozen.gamma_min == 0.01


def test_lower_gamma_binary_still_requires_original_stability_and_fraction_gates():
    cfg = replace(ResolutionSearchConfig(), gamma_min=0.001)
    result = coverage_search(
        cfg, lambda g: candidate(g, 2, False, ("UNSTABLE", "MAX_CHILD_FRACTION")), protocol=FAIR
    )
    assert result.selected is None
    assert result.search_status == "REJECTED_ALL_TESTED"
    assert len(result.candidates) <= 24
    assert all(c.violations == ("UNSTABLE", "MAX_CHILD_FRACTION") for c in result.candidates)
