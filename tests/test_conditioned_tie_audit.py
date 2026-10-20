import pytest

from benchmarks.og_extraction.conditioned_tie_audit import audit_tree


def fixture():
    nodes = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
        dict(cluster_id=3, parent_id=1),
        dict(cluster_id=4, parent_id=1),
    ]
    members = [(0, 3), (1, 3), (2, 4), (3, 4), (4, 2), (5, 2), (6, 2), (7, 2)]
    species = {0: 0, 1: 1, 2: 2, 3: 3, 4: 0, 5: 1, 6: 2, 7: 3}
    return nodes, members, species


def test_nonzero_exact_local_tie_keeps_parent():
    n, m, s = fixture()
    result = audit_tree(n, m, [(0, 1, 0.1), (2, 3, 0.3), (4, 5, 0.2), (6, 7, 0.4)], s)
    row = next(r for r in result["local_rows"] if r["cluster_id"] == 1)
    assert row["exact_keep_score"] > 0
    assert row["exact_gap_sign"] == 0
    assert row["exact_tie"] and row["exact_keep"]
    assert 1 in result["exact_selected"]
    assert result["exact_objective_identity_verified"]
    for r in result["local_rows"]:
        if r["exact_split_minus_keep"] is not None:
            assert sum(v["gain"] for v in r["exact_child_boundary"].values()) == pytest.approx(
                r["exact_split_minus_keep"]
            )
    assert set(result["exact_prediction"]) == set(range(8))


def test_genuine_positive_split_and_zero_root():
    n = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
    ]
    m = [(0, 1), (1, 1), (2, 2), (3, 2)]
    result = audit_tree(n, m, [(0, 1, 1.0), (2, 3, 1.0)], {0: 0, 1: 1, 2: 0, 3: 1})
    root = next(r for r in result["local_rows"] if r["cluster_id"] == 0)
    assert root["exact_gap_sign"] == 1
    assert result["exact_selected"] == [1, 2]
    assert result["exact_score"] == pytest.approx(0.5)
    zero = audit_tree(n, m, [], {0: 0, 1: 1, 2: 0, 3: 1})
    assert zero["exact_selected"] == [0]
    assert zero["exact_score"] == 0


def test_edge_order_does_not_change_exact_cut_or_signs():
    n, m, s = fixture()
    edges = [(0, 1, 0.1), (2, 3, 0.3), (4, 5, 0.2), (6, 7, 0.4)]
    a, b = audit_tree(n, m, edges, s), audit_tree(n, m, list(reversed(edges)), s)
    assert a["exact_selected"] == b["exact_selected"]
    assert [r["exact_gap_sign"] for r in a["local_rows"]] == [
        r["exact_gap_sign"] for r in b["local_rows"]
    ]
