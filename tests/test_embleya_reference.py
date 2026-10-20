from itertools import combinations

import pytest

from benchmarks.og_extraction.embleya_reference import compare


@pytest.mark.parametrize(
    "prediction",
    [
        {0: "a", 1: "a", 2: "b", 3: "b"},
        {0: "a", 1: "a", 2: "a", 3: "a"},
        {0: "a", 1: "b", 2: "c", 3: "d"},
        {0: "a", 2: "a"},
    ],
)
def test_contingency_matches_enumerated_pairs(prediction):
    reference = {0: "r", 1: "r", 2: "s", 3: "s"}
    result = compare(reference, prediction, dict.fromkeys(reference, 0))
    expected = {(a, b) for a, b in combinations(reference, 2) if reference[a] == reference[b]}
    observed = {(a, b) for a, b in combinations(prediction, 2) if prediction[a] == prediction[b]}
    assert result["true_positive_pairs"] == len(expected & observed)
    assert result["discordant_pairs"] == len(observed - expected)
    assert result["lost_pairs"] == len(expected - observed)


def test_component_ceiling_and_extra_predictions():
    result = compare({0: "r", 1: "r", 2: "r"}, {0: "p", 1: "p", 2: "q", 3: "p"}, {0: 0, 1: 0, 2: 1})
    assert result["reference_pairs_across_components"] == 2
    assert result["component_pair_recall_ceiling"] == 1 / 3
    assert result["excluded_predictions"] == 1
    assert result["missing_predictions"] == 0
    assert result["bcubed_precision"] == 1
    assert result["bcubed_recall"] == pytest.approx(5 / 9)
