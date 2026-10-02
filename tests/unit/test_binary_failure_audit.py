"""Finite diagnostic evidence is distinct from admissible search points."""

import math

from benchmarks.og_extraction.binary_failure_audit import diagnostic_points
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig


def test_out_of_range_probe_is_explicit_and_only_for_lower_boundary_already_multigroup():
    config = ResolutionSearchConfig()
    points = diagnostic_points(
        config, [dict(gamma=0.01, child_count=3), dict(gamma=10.0, child_count=132)]
    )
    below = [(g, scope) for g, scope in points if g < 0.01]
    assert len(below) == 32
    assert all(scope == "BELOW_MIN_DIAGNOSTIC_ONLY" for _, scope in below)
    assert all(0.01 <= g <= 10 for g, scope in points if scope != "BELOW_MIN_DIAGNOSTIC_ONLY")
    assert len(points) == 161


def test_transition_diagnostic_is_finite_deduplicated_and_within_original_bounds():
    points = diagnostic_points(
        ResolutionSearchConfig(), [dict(gamma=0.94, child_count=1), dict(gamma=0.96, child_count=3)]
    )
    assert len(points) <= 194
    assert all(0.01 <= g <= 10 for g, _ in points)
    assert any(math.isclose(g, 0.9584375, rel_tol=1e-12) for g, _ in points)
