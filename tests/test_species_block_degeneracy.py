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
