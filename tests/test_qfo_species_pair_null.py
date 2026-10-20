import pytest

from benchmarks.qfo.cut_species_null import decompose
from benchmarks.qfo.species_pair_null import conditioned_null


def test_species_segregation_deficit_disappears():
    pred = {0: 0, 1: 0, 2: 1, 3: 1}
    species = dict(pred)
    edges = [(0, 1, 1.0), (2, 3, 1.0)]
    assert decompose(pred, species, edges)["gain"] == pytest.approx(0.5)
    assert conditioned_null(pred, species, edges)["gain"] == pytest.approx(0)


def test_single_species_reduces_to_unconditional_null():
    pred = {0: 0, 1: 0, 2: 1, 3: 1}
    species = dict.fromkeys(pred, 0)
    edges = [(0, 1, 2.0), (0, 2, 1.0), (2, 3, 3.0)]
    assert conditioned_null(pred, species, edges)["gain"] == pytest.approx(
        decompose(pred, species, edges)["gain"]
    )


def test_bipartite_expectation_retains_mixed_species_cut_signal():
    pred = {0: 0, 1: 0, 2: 1, 3: 1}
    species = {0: 0, 1: 1, 2: 0, 3: 1}
    result = conditioned_null(pred, species, [(0, 1, 1.0), (2, 3, 1.0)])
    assert result["gain"] == pytest.approx(0.5)
    assert result["cross_species_gain"] == pytest.approx(0.5)
    assert result["block_conservation_verified"]


def test_cross_block_cut_is_neutral_when_groups_are_species():
    result = conditioned_null({0: 0, 1: 1}, {0: 0, 1: 1}, [(0, 1, 1.0)])
    assert result["gain"] == 0
    assert result["species_blocks"][0]["expected_between_groups_weight"] == 1


def test_scaling_zero_blocks_and_fixed_prediction():
    pred = {0: 0, 1: 0, 2: 1, 3: 1}
    original = dict(pred)
    species = {0: 0, 1: 1, 2: 0, 3: 1}
    edges = [(0, 1, 1.0), (0, 2, 0.0), (2, 3, 2.0)]
    a = conditioned_null(pred, species, edges)
    b = conditioned_null(pred, species, [(u, v, 7 * w) for u, v, w in edges])
    assert a["gain"] == pytest.approx(b["gain"])
    assert pred == original
    assert conditioned_null(pred, species, [])["gain"] == 0
    assert conditioned_null(dict.fromkeys(pred, 0), species, edges)["gain"] == 0


def test_negative_signed_gain_is_not_clamped():
    pred = {0: 0, 1: 1, 2: 0, 3: 1}
    species = dict.fromkeys(pred, 0)
    result = conditioned_null(pred, species, [(0, 1, 1.0), (2, 3, 1.0)])
    assert result["gain"] == pytest.approx(-0.5)


def test_multiple_species_blocks_match_explicit_pair_expectations():
    pred = {0: 0, 1: 0, 2: 1, 3: 1, 4: 2, 5: 2}
    species = {0: 0, 1: 1, 2: 0, 3: 1, 4: 0, 5: 2}
    edges = [
        (0, 1, 2.0),
        (0, 2, 3.0),
        (1, 2, 1.0),
        (2, 3, 4.0),
        (3, 4, 2.0),
        (4, 5, 5.0),
        (0, 5, 1.0),
    ]
    # Species 0/1: W=9; endpoint strengths [2,5,2] and [3,6,0].
    # Between expectation = (9*9 - (2*3 + 5*6 + 2*0))/9 = 5.
    # Species 0/0: only group 0--1 edge, W=3: expectation 3*3/(2*3)=1.5.
    # Species 0/2: W=6; strengths [1,0,5] and [0,0,6]: expectation 1.
    # Observed between groups: 3+1+2+1=7; total W=18.
    result = conditioned_null(pred, species, edges)
    assert result["gain"] == pytest.approx((5 + 1.5 + 1 - 7) / 18)
    assert sum(r["gain_contribution"] for r in result["species_blocks"]) == pytest.approx(
        result["gain"]
    )
