import pytest

from benchmarks.qfo.cut_species_null import decompose


def test_split_by_species_has_cross_species_null_deficit():
    result = decompose(
        {0: 0, 1: 0, 2: 1, 3: 1}, {0: 0, 1: 0, 2: 1, 3: 1}, [(0, 1, 1.0), (2, 3, 1.0)]
    )
    assert result["gain"] == pytest.approx(0.5)
    assert result["species_terms"]["same_species"]["gain_contribution"] == 0
    assert result["species_terms"]["cross_species"]["gain_contribution"] == pytest.approx(0.5)


def test_mixed_species_subgroups_decompose_additively():
    result = decompose(
        {0: 0, 1: 0, 2: 1, 3: 1}, {0: 0, 1: 1, 2: 0, 3: 1}, [(0, 1, 1.0), (2, 3, 1.0)]
    )
    assert result["gain"] == pytest.approx(0.5)
    assert result["species_terms"]["same_species"]["gain_contribution"] == pytest.approx(0.25)
    assert result["species_terms"]["cross_species"]["gain_contribution"] == pytest.approx(0.25)
    assert result["group_species_pair_overlap"][0]["shared_species"] == 2


def test_single_group_and_zero_weight_explicit():
    result = decompose({0: 0, 1: 0}, {0: 0, 1: 1}, [(0, 1, 1.0)])
    assert result["gain"] == 0
    assert result["group_species_pair_overlap"] == []
    zero = decompose({0: 0, 1: 1}, {0: 0, 1: 1}, [])
    assert zero["gain"] == 0
    assert zero["total_weight"] == 0


def test_negative_species_contributions_are_retained():
    result = decompose({0: 0, 1: 1}, {0: 0, 1: 1}, [(0, 1, 1.0)])
    assert result["species_terms"]["cross_species"]["gain_contribution"] == pytest.approx(-0.5)
    assert result["gain"] == pytest.approx(-0.5)


def test_weight_rescaling_preserves_gain():
    pred = {0: 0, 1: 0, 2: 1, 3: 1}
    species = {0: 0, 1: 1, 2: 0, 3: 1}
    a = decompose(pred, species, [(0, 1, 1.0), (2, 3, 1.0)])
    b = decompose(pred, species, [(0, 1, 7.0), (2, 3, 7.0)])
    assert a["gain"] == pytest.approx(b["gain"])
