from fractions import Fraction

import pytest

from benchmarks.og_extraction.complete_objective_audit import annotate_losses, summarize
from benchmarks.og_extraction.species_block_degeneracy import annotate_reference, profile
from benchmarks.og_extraction.species_pair_graph import species_pair_cut


def audit(nodes, members, edges, species, reference):
    cut = species_pair_cut(nodes, members, edges, species)
    result = profile(
        nodes,
        members,
        edges,
        species,
        set(cut["selected"]),
        context_flow=True,
        objective_decomposition=True,
    )
    rows = result["legal_internal_nodes"]
    annotate_reference(rows, nodes, members, reference)
    annotate_losses(rows, nodes, members, reference)
    return rows, cut


def part(row, stage, kind):
    return next(
        p
        for p in row["objective_decomposition"]["parts"]
        if p["stage"] == stage and p["species_channel"] == kind
    )


def test_same_and_cross_actual_null_channels_with_hand_computed_complete_cut():
    nodes = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
    ]
    members = [(0, 1), (1, 1), (2, 2), (3, 2)]
    # Each child has a unit same-species edge, plus one unit cross-species bridge.
    rows, cut = audit(
        nodes,
        members,
        [(0, 1, 1.0), (2, 3, 1.0), (0, 2, 1.0)],
        {0: 0, 1: 0, 2: 1, 3: 1},
        dict.fromkeys(range(4), "A"),
    )
    r = rows[0]
    same = part(r, "direct", "same_species")
    cross = part(r, "direct", "cross_species")
    assert same["actual"] == same["expected"] == 0
    assert cross["actual"] == cross["expected"] == 1
    assert cut["selected"] == [0]
    assert part(r, "optimal", "cross_species")["actual"] == 0
    assert r["objective_decomposition"]["local_complete_tp"] == 0
    assert r["objective_decomposition"]["actual_direct_lca_tp"] == 0


def test_same_species_gain_from_disconnected_children_and_no_cross_channel():
    nodes = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
    ]
    members = [(0, 1), (1, 1), (2, 2), (3, 2)]
    rows, _ = audit(
        nodes,
        members,
        [(0, 1, 1.0), (2, 3, 1.0)],
        dict.fromkeys(range(4), 0),
        dict.fromkeys(range(4), "A"),
    )
    r = rows[0]
    same = part(r, "optimal", "same_species")
    assert same["actual"] == 0 and same["expected"] == 1 and same["gain"] == 0.5
    assert part(r, "optimal", "cross_species")["gain_exact"] == "0"
    assert r["objective_decomposition"]["category"] == "internal_positive"
    assert r["objective_decomposition"]["actual_direct_lca_tp"] == 4


def test_direct_and_descendant_terms_are_distinct_and_telescope_without_nested_losses():
    nodes = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
        dict(cluster_id=3, parent_id=1),
        dict(cluster_id=4, parent_id=1),
    ]
    members = [(0, 3), (1, 3), (2, 4), (3, 4), (4, 2), (5, 2)]
    rows, _ = audit(
        nodes,
        members,
        [(0, 1, 1.0), (2, 3, 1.0), (4, 5, 1.0)],
        dict.fromkeys(range(6), 0),
        dict.fromkeys(range(6), "A"),
    )
    root = next(r for r in rows if r["cluster_id"] == 0)
    assert part(root, "direct", "same_species")["gain"] == pytest.approx(4 / 9)
    assert part(root, "descendant", "same_species")["gain"] == pytest.approx(2 / 9)
    assert part(root, "optimal", "same_species")["gain"] == pytest.approx(2 / 3)
    assert root["objective_decomposition"]["actual_direct_lca_tp"] == 8
    assert root["objective_decomposition"]["local_complete_tp"] == 12
    summary = summarize(
        [
            dict(
                component_id=0,
                source_cluster_id=0,
                reference_cohort="clean",
                objective=dict(block_degeneracy=dict(legal_internal_nodes=rows)),
            )
        ]
    )
    assert sum(s["actual_direct_lca_tp"] for s in summary["strata"]) == 12
    assert summary["sources"][0]["actual_lca_channels"][0]["expected"] == 2
    assert not summary["matching_performed"] and not summary["objective_changed"]


def test_random_channel_exact_identities_and_zero_graph():
    import random

    rng = random.Random(42)
    nodes = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
        dict(cluster_id=3, parent_id=1),
        dict(cluster_id=4, parent_id=1),
    ]
    members = [(0, 3), (1, 3), (2, 4), (3, 4), (4, 2), (5, 2)]
    for _ in range(40):
        species = {p: rng.randrange(3) for p, _ in members}
        edges = [
            (u, v, rng.randrange(1, 7) / 4)
            for u in range(6)
            for v in range(u + 1, 6)
            if rng.random() < 0.4
        ]
        rows, cut = audit(nodes, members, edges, species, {p: rng.randrange(2) for p, _ in members})
        for row in rows:
            d = row["objective_decomposition"]
            for kind in ("same_species", "cross_species"):
                a = part(row, "direct", kind)
                b = part(row, "descendant", kind)
                c = part(row, "forced_split", kind)
                assert Fraction(a["gain_exact"]) + Fraction(b["gain_exact"]) == Fraction(
                    c["gain_exact"]
                )
            assert sum(
                Fraction(part(row, "direct", k)["gain_exact"])
                for k in ("same_species", "cross_species")
            ) == Fraction(d["direct_gain_exact"])
    rows, _ = audit(nodes, members, [], dict.fromkeys(range(6), 0), dict.fromkeys(range(6), "A"))
    assert all(r["objective_decomposition"]["local_keep"] for r in rows)
    assert all(p["gain_exact"] == "0" for r in rows for p in r["objective_decomposition"]["parts"])


def test_direct_nonpositive_descendant_split_fixture():
    nodes = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
        dict(cluster_id=3, parent_id=1),
        dict(cluster_id=4, parent_id=1),
    ]
    members = [(0, 3), (1, 3), (2, 4), (3, 4), (4, 2), (5, 2)]
    # This fixture is filled from a graph with independently inspected direct/optimal gain.
    rows, _ = audit(
        nodes,
        members,
        [(0, 1, 3.0), (0, 4, 1.0), (2, 5, 2.0)],
        dict.fromkeys(range(6), 0),
        dict.fromkeys(range(6), "A"),
    )
    root = rows[0]
    d = root["objective_decomposition"]
    assert d["direct_gain_exact"] == "-1/8"
    assert d["forced_split_gain_exact"] == d["optimal_gain_exact"] == "5/72"
    assert part(root, "descendant", "same_species")["gain_exact"] == "7/36"
    assert d["category"] == "nonpositive_direct_with_descendant_optimum"
    assert root["dp_state"] == "active_split"
    assert d["actual_direct_lca_tp"] == 8
