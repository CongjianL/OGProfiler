import pytest

from benchmarks.qfo.internal_boundary import diagnostic_pairs, inspect


def inputs():
    return [
        dict(cluster_id=0, parent_id=None, n_genes=4),
        dict(cluster_id=1, parent_id=0, n_genes=2),
        dict(cluster_id=2, parent_id=0, n_genes=2),
    ], [(0, 1), (1, 1), (2, 2), (3, 2)]


def test_disconnected_dense_children_split_under_fixed_objective():
    nodes, members = inputs()
    result, pred = inspect(nodes, members, [(0, 1, 1.0), (2, 3, 1.0)])
    assert result["split"]
    assert result["merge_score"] == pytest.approx(0)
    assert result["cut_gain"] == pytest.approx(0.5)
    assert set(pred) == {0, 1, 2, 3}
    assert diagnostic_pairs({0: "A", 1: "A", 2: "B", 3: "B"}, pred) == dict(tp=2, fp=0, lost=0)


def test_homogeneous_clique_keeps_complete_merge():
    nodes, members = inputs()
    edges = [(u, v, 1.0) for u in range(4) for v in range(u + 1, 4)]
    result, pred = inspect(nodes, members, edges)
    assert not result["split"]
    assert result["cut_gain"] == pytest.approx(0)
    assert len(set(pred.values())) == 1


def test_zero_weight_ties_keep_parent_and_scaling_invariance():
    nodes, members = inputs()
    zero, _ = inspect(nodes, members, [])
    assert not zero["split"]
    a, _ = inspect(nodes, members, [(0, 1, 1.0), (2, 3, 1.0)])
    b, _ = inspect(nodes, members, [(0, 1, 7.0), (2, 3, 7.0)])
    assert a["cut_gain"] == pytest.approx(b["cut_gain"])
