"""Finite bank and unchanged V1 event semantics for target diagnostics."""

from dataclasses import replace
from types import SimpleNamespace

import pytest

from benchmarks.og_extraction.binary_target_audit import candidate_bank, event
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig


def test_common_bank_preserves_original_points_and_deduplicates_within_budget():
    bank = candidate_bank([0.08, 0.08, 0.625], ResolutionSearchConfig())
    assert len(bank) == 20
    assert bank.count(0.08) == 1
    assert 0.625 in bank
    assert 10.0 in bank
    assert bank == sorted(bank)
    with pytest.raises(ValueError, match="exceeds candidate budget"):
        candidate_bank([], replace(ResolutionSearchConfig(), max_candidate_evaluations=18))


def test_binary_unique_species_event_is_i_but_overlapping_binary_stays_iii1():
    unique = SimpleNamespace(membership=(0, 0, 1), child_count=2)
    assert event(unique, [0, 1, 2], True) == "I"
    overlap = SimpleNamespace(membership=(0, 0, 0, 1), child_count=2)
    assert event(overlap, [0, 1, 2, 1], False) == "III-1"


def test_single_copy_k_way_does_not_gain_binary_event_qualification():
    candidate = SimpleNamespace(membership=tuple(range(11)), child_count=11)
    assert event(candidate, list(range(11)), True) == "III-3"
