from benchmarks.og_extraction.species_block_degeneracy import annotate_reference, profile
from benchmarks.og_extraction.species_pair_graph import species_pair_cut


def fixture():
    n = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
        dict(cluster_id=3, parent_id=1),
        dict(cluster_id=4, parent_id=1),
    ]
    m = [(0, 3), (1, 3), (2, 4), (3, 2)]
    s = {0: 0, 1: 1, 2: 1, 3: 0}
    return n, m, s


def test_connected_single_side_equal_not_missing_edge():
    n, m, s = fixture()
    r = profile(n, m, [(0, 1, 2.0), (0, 2, 1.0)], s, {0})
    row = next(r for r in r["legal_internal_nodes"] if r["cluster_id"] == 1)
    assert row["single_side_connected_equal_blocks"] == 1
    assert row["one_side_opposite_distributed_blocks"] == 1
    assert row["connected_cross_weight"] == 1
    assert row["single_side_equal_weight_fraction"] == 1
    assert row["direct_gain_sign"] == 0
    assert row["dp_state"] == "inactive"
    assert row["blocks"][0]["expected"] == row["blocks"][0]["observed"] == 1


def test_external_node_context_breaks_single_side_carrier():
    n, m, s = fixture()
    r = profile(n, m, [(0, 1, 2.0), (0, 2, 1.0), (2, 3, 0.5)], s, {1, 2})
    row = next(r for r in r["legal_internal_nodes"] if r["cluster_id"] == 1)
    assert row["single_side_connected_equal_blocks"] == 0
    assert row["decisive_blocks"] == 1
    assert row["blocks"][0]["nonzero_deficit"]


def test_species_segregation_classified_as_both_sides_concentrated():
    n = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
    ]
    r = profile(n, [(0, 1), (1, 2)], [(0, 1, 1.0)], {0: 0, 1: 1}, {0})
    row = r["legal_internal_nodes"][0]
    assert row["both_sides_concentrated_blocks"] == 1
    assert row["one_side_opposite_distributed_blocks"] == 0
    assert row["single_side_connected_equal_blocks"] == 1
    zero = profile(n, [(0, 1), (1, 2)], [], {0: 0, 1: 1}, {0})
    assert zero["legal_internal_nodes"][0]["cross_connected_blocks"] == 0
    assert zero["legal_internal_nodes"][0]["single_side_equal_weight_fraction"] is None


def test_lca_pairs_attribute_selected_complete_cut_losses():
    n, m, s = fixture()
    e = [(0, 1, 2.0), (0, 2, 1.0), (2, 3, 0.5)]
    cut = species_pair_cut(n, m, e, s)
    r = profile(n, m, e, s, set(cut["selected"]))
    annotate_reference(r["legal_internal_nodes"], n, m, {0: "A", 1: "A", 2: "A", 3: "B"})
    lost = sum(
        row["direct_removed_tp"]
        for row in r["legal_internal_nodes"]
        if row["dp_state"] == "active_split"
    )
    joint = {}
    ref = {0: "A", 1: "A", 2: "A", 3: "B"}
    for row in cut["members"]:
        key = ref[row["protein_id"]], row["cluster_id"]
        joint[key] = joint.get(key, 0) + 1
    tp = sum(v * (v - 1) // 2 for v in joint.values())
    assert lost == 3 - tp
    assert all(row["objective_identity_verified"] for row in r["legal_internal_nodes"])


def block_at_one(n, m, e, s):
    rows = profile(n, m, e, s, {1, 2}, context_flow=True)["legal_internal_nodes"]
    row = next(r for r in rows if r["cluster_id"] == 1)
    return row, row["blocks"][0]


def test_context_cancelled_block_is_internal_only_without_boundary():
    n, m, s = fixture()
    row, b = block_at_one(n, m, [(0, 1, 2.0), (0, 2, 1.0)], s)
    assert b["expected_internal_internal"] == b["observed"] == 1
    assert b["expected_internal_external"] == b["expected_external_external"] == 0
    assert b["node_boundary_weight"] == b["outside_node_weight"] == 0
    assert row["context_groups"][0]["blocks"] == 1
    assert not b["source_context_changes_sign"]


def test_boundary_context_turns_internal_negative_into_positive():
    n, m, s = fixture()
    e = [(0, 1, 2.0), (0, 2, 1.0), (2, 3, 0.5)]
    row, b = block_at_one(n, m, e, s)
    assert b["internal_only_sign"] == -1
    assert b["total_deficit_sign"] == 1
    assert b["source_context_changes_sign"]
    assert b["node_internal_weight"] == 3
    assert b["node_boundary_weight"] == 0.5
    assert b["outside_node_weight"] == 0
    assert row["context_groups"][1]["internal_nonpositive_total_positive"] == 1
    assert (
        profile(n, m, e, s, {1, 2})["legal_internal_nodes"][1]["direct_gain"] == row["direct_gain"]
    )


def test_context_cross_external_and_outside_edge_conservation():
    import pytest

    n, m, s = fixture()
    m = m + [(4, 2)]
    s = dict(s)
    s[4] = 1
    e = [(0, 1, 2.0), (0, 2, 1.0), (2, 3, 0.5), (0, 4, 0.25), (3, 4, 0.75)]
    _, b = block_at_one(n, m, e, s)
    assert b["expected_external_external"] == pytest.approx(0.125 / 4.5)
    assert b["node_internal_weight"] == 3
    assert b["node_boundary_weight"] == 0.75
    assert b["outside_node_weight"] == 0.75
    assert b["source_block_weight"] == 4.5
    assert b["context_identity_verified"]


def test_same_species_context_endpoint_factors_and_empty_graph():
    import pytest

    n, m, _ = fixture()
    s = {p: 0 for p, _ in m}
    e = [(0, 1, 2.0), (0, 2, 1.0), (2, 3, 0.5)]
    _, b = block_at_one(n, m, e, s)
    assert b["expected_internal_internal"] == pytest.approx(10 / 14)
    assert b["expected_internal_external"] == pytest.approx(5 / 14)
    assert b["expected_external_external"] == 0
    assert b["node_internal_weight"] == 3
    assert b["node_boundary_weight"] == 0.5
    assert b["outside_node_weight"] == 0
    empty = profile(n, m, [], s, {0}, context_flow=True)
    assert all(
        all(g["blocks"] == 0 for g in r["context_groups"]) for r in empty["legal_internal_nodes"]
    )


def test_context_exact_identities_across_random_species_and_edges():
    import random

    rng = random.Random(42)
    n, m, _ = fixture()
    m += [(4, 2), (5, 4)]
    for _ in range(30):
        s = {p: rng.randrange(3) for p, _ in m}
        e = [
            (u, v, rng.randrange(1, 9) / 4)
            for u in range(6)
            for v in range(u + 1, 6)
            if rng.random() < 0.5
        ]
        cut = species_pair_cut(n, m, e, s)
        original = profile(n, m, e, s, set(cut["selected"]))
        context = profile(n, m, e, s, set(cut["selected"]), context_flow=True)
        for before, after in zip(
            original["legal_internal_nodes"], context["legal_internal_nodes"], strict=True
        ):
            assert before["direct_gain"] == after["direct_gain"]
            assert before["dp_state"] == after["dp_state"]
            assert after["context_identity_verified"]
