import math

import pytest

from benchmarks.og_extraction.conditioned_tie_audit import audit_tree
from benchmarks.og_extraction.exact_species_pair_graph import exact_scores
from benchmarks.og_extraction.species_pair_graph import species_pair_cut
from benchmarks.og_extraction.tree_cut import optimal_cut


def fixture():
    nodes = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
    ]
    return nodes, [(0, 1), (1, 1), (2, 2), (3, 2)], {0: 0, 1: 1, 2: 0, 3: 1}


def test_positive_gain_below_previous_epsilon_is_preserved():
    n, m, s = fixture()
    x = math.nextafter(1.0, 0.0)
    edges = [(0, 1, 1.0), (2, 3, 1.0), (0, 3, x), (1, 2, x)]
    raw = species_pair_cut(n, m, edges, s, exact=False)
    repaired = species_pair_cut(n, m, edges, s)
    assert raw["selected"] == [0]
    assert repaired["selected"] == [1, 2]
    assert 0 < repaired["score"] < 1e-12
    assert repaired["objective_identity_verified"]
    scores, _ = exact_scores(n, m, edges, s)
    assert scores[1] + scores[2] > scores[0]
    assert optimal_cut(n, scores).selected == (1, 2)


def test_nonzero_exact_tie_and_repaired_cut_match_auditor():
    n = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
        dict(cluster_id=3, parent_id=1),
        dict(cluster_id=4, parent_id=1),
    ]
    m = [(0, 3), (1, 3), (2, 4), (3, 4), (4, 2), (5, 2), (6, 2), (7, 2)]
    s = {0: 0, 1: 1, 2: 2, 3: 3, 4: 0, 5: 1, 6: 2, 7: 3}
    e = [(0, 1, 0.1), (2, 3, 0.3), (4, 5, 0.2), (6, 7, 0.4)]
    cut = species_pair_cut(n, m, e, s)
    audit = audit_tree(n, m, e, s)
    assert cut["selected"] == audit["exact_selected"]
    assert 1 in cut["selected"]
    assert cut["score"] == pytest.approx(audit["exact_score"])


def test_exact_cut_is_order_invariant_and_zero_keeps_root():
    n, m, s = fixture()
    e = [(0, 1, 1.0), (2, 3, 1.0), (0, 2, 0.2)]
    a, b = species_pair_cut(n, m, e, s), species_pair_cut(n, m, list(reversed(e)), s)
    assert a["selected"] == b["selected"]
    assert a["score"] == b["score"]
    assert species_pair_cut(n, m, [], s)["selected"] == [0]


def test_21396_block_detail_matches_boundary_gain():
    n, m, s = fixture()
    n[0]["cluster_id"] = 21396
    n[1]["parent_id"] = n[2]["parent_id"] = 21396
    audit = audit_tree(n, m, [(0, 1, 1.0), (2, 3, 1.0)], s)
    row = next(r for r in audit["local_rows"] if r["cluster_id"] == 21396)
    blocks = row["exact_child_boundary"]["species_blocks"]
    assert sum(b["gain"] for b in blocks) == pytest.approx(row["exact_split_minus_keep"])
    assert blocks[0]["source_block_weight"] == 2
    assert len(blocks[0]["endpoint_strengths"]) == 2
