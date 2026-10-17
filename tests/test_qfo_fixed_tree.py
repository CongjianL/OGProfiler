import itertools
import random
from collections import Counter

import pytest

from benchmarks.qfo.fixed_tree import boundary_losses


def test_pure_lca_is_available_selection_opportunity():
    nodes = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
    ]
    result = boundary_losses(nodes, [(10, 1), (11, 2)], {10: "A", 11: "A"}, {10: "P1", 11: "P2"})
    assert result["pure_assigned_lca_tree_available"] == 1
    assert result["within_component_lost_pairs"] == 1


def test_mixed_lca_requires_contamination_but_does_not_prove_cause():
    nodes = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
    ]
    result = boundary_losses(
        nodes,
        [(10, 1), (11, 2), (12, 2)],
        {10: "A", 11: "A", 12: "B"},
        {10: "P1", 11: "P2", 12: "P2"},
    )
    assert result["mixed_assigned_lca_requires_contamination"] == 1
    assert result["pure_assigned_lca_tree_available"] == 0


def test_loss_inside_indivisible_structural_leaf_is_separate():
    result = boundary_losses(
        [dict(cluster_id=0, parent_id=None)],
        [(10, 0), (11, 0)],
        {10: "A", 11: "A"},
        {10: "P1", 11: "P2"},
    )
    assert result["within_structural_leaf"] == 1


def test_lca_contingency_counts_match_explicit_pairs():
    nodes = [dict(cluster_id=i, parent_id=None if i == 0 else (i - 1) // 2) for i in range(7)]
    membership = [(p, 3 + p // 2) for p in range(8)]
    leaf = dict(membership)

    def path(c):
        result = [c]
        while c:
            c = (c - 1) // 2
            result.append(c)
        return result

    rng = random.Random(42)
    for _ in range(40):
        ref = {p: rng.choice("ABC") for p in range(8) if rng.random() > 0.2}
        actual = {p: rng.choice("XYZ") for p in range(8)}
        brute = Counter()
        for a, b in itertools.combinations(ref, 2):
            if ref[a] != ref[b] or actual[a] == actual[b]:
                continue
            ancestors = set(path(leaf[b]))
            lca = next(c for c in path(leaf[a]) if c in ancestors)
            labels = {ref[p] for p in ref if lca in path(leaf[p])}
            key = (
                "within_structural_leaf"
                if lca >= 3
                else "pure_assigned_lca_tree_available"
                if len(labels) == 1
                else "mixed_assigned_lca_requires_contamination"
            )
            brute[key] += 1
        result = boundary_losses(nodes, membership, ref, actual)
        for key in (
            "within_structural_leaf",
            "pure_assigned_lca_tree_available",
            "mixed_assigned_lca_requires_contamination",
        ):
            assert result.get(key, 0) == brute[key]
        assert result["within_component_lost_pairs"] == sum(brute.values())


def test_duplicate_membership_fails():
    with pytest.raises(ValueError, match="Invalid terminal"):
        boundary_losses(
            [dict(cluster_id=0, parent_id=None)], [(10, 0), (10, 0)], {10: "A"}, {10: "P"}
        )
