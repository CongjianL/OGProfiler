import itertools

import pytest

from benchmarks.og_extraction.species_pair_graph import species_pair_cut
from benchmarks.og_extraction.tree_cut import topology
from benchmarks.qfo.internal_boundary import inspect
from benchmarks.qfo.species_pair_null import conditioned_null


def fixture():
    nodes = [
        dict(cluster_id=0, parent_id=None, n_genes=4),
        dict(cluster_id=1, parent_id=0, n_genes=2),
        dict(cluster_id=2, parent_id=0, n_genes=2),
    ]
    return nodes, [(0, 1), (1, 1), (2, 2), (3, 2)]


def test_species_segregation_changes_optimal_cut_not_posthoc_filter():
    n, m = fixture()
    e = [(0, 1, 1.0), (2, 3, 1.0)]
    assert inspect(n, m, e)[0]["split"]
    cut = species_pair_cut(n, m, e, {0: 0, 1: 0, 2: 1, 3: 1})
    assert cut["selected"] == [0]
    assert cut["score"] == 0


def test_mixed_species_boundary_still_splits_and_scales():
    n, m = fixture()
    s = {0: 0, 1: 1, 2: 0, 3: 1}
    e = [(0, 1, 1.0), (2, 3, 1.0)]
    a = species_pair_cut(n, m, e, s)
    b = species_pair_cut(n, m, [(u, v, 7 * w) for u, v, w in e], s)
    assert a["selected"] == b["selected"] == [1, 2]
    assert a["score"] == pytest.approx(0.5)
    assert b["score"] == pytest.approx(a["score"])
    assert a["objective_identity_verified"]


def test_single_species_matches_original_and_zero_edges_keep_parent():
    n, m = fixture()
    e = [(0, 1, 2.0), (2, 3, 3.0), (0, 2, 0.2)]
    s = dict.fromkeys(range(4), 0)
    a = species_pair_cut(n, m, e, s)
    old, _ = inspect(n, m, e)
    assert a["selected"] == old["selected"]
    assert a["score"] == pytest.approx(old["best_cut_score"])
    assert species_pair_cut(n, m, [], s)["selected"] == [0]


def test_dp_matches_exhaustive_complete_cuts_and_independent_gain():
    nodes = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
        dict(cluster_id=3, parent_id=1),
        dict(cluster_id=4, parent_id=1),
        dict(cluster_id=5, parent_id=2),
        dict(cluster_id=6, parent_id=2),
    ]
    membership = [(0, 3), (1, 4), (2, 5), (3, 6)]
    species = {0: 0, 1: 1, 2: 0, 3: 1}
    edges = [(0, 1, 4.0), (0, 2, 0.5), (1, 3, 0.3), (2, 3, 2.0)]
    result = species_pair_cut(nodes, membership, edges, species)
    _, children, _ = topology(nodes)

    def cuts(c):
        yield (c,)
        if children[c]:
            for combo in itertools.product(*(list(cuts(ch)) for ch in children[c])):
                yield sum(combo, ())

    scores = {r["cluster_id"]: r["keep_score"] for r in result["node_scores"]}
    enumerated = [sum(scores[c] for c in cut) for cut in cuts(0)]
    assert result["score"] == pytest.approx(max(enumerated))
    prediction = {r["protein_id"]: r["cluster_id"] for r in result["members"]}
    assert set(prediction) == set(range(4))
    assert result["score"] == pytest.approx(conditioned_null(prediction, species, edges)["gain"])


def test_duplicate_or_outside_edges_rejected_and_leaves_indivisible():
    n, m = fixture()
    s = dict.fromkeys(range(4), 0)
    with pytest.raises(ValueError):
        species_pair_cut(n, m, [(0, 1, 1.0), (0, 1, 1.0)], s)
    with pytest.raises(ValueError):
        species_pair_cut(n, m, [(0, 9, 1.0)], s)
    result = species_pair_cut(n, m, [(0, 2, 1.0)], s)
    assert set(result["selected"]).issubset({0, 1, 2})
