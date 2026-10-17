from copy import deepcopy

import pytest

from benchmarks.og_extraction.positive_boundary import annotate_positive, compare_positive
from benchmarks.og_extraction.species_block_degeneracy import annotate_reference, profile
from benchmarks.og_extraction.species_pair_graph import species_pair_cut


def run_case(nodes, members, edges, reference):
    species = {p: 0 for p, _ in members}
    cut = species_pair_cut(nodes, members, edges, species)
    result = profile(
        nodes,
        members,
        edges,
        species,
        set(cut["selected"]),
        context_flow=True,
        positive_boundary=True,
    )
    rows = result["legal_internal_nodes"]
    annotate_reference(rows, nodes, members, reference)
    annotate_positive(rows, nodes, members, reference)
    return dict(
        component_id=0,
        source_cluster_id=0,
        reference_cohort="polluted",
        objective=dict(block_degeneracy=result),
    )


def test_same_fixed_boundary_objective_distinguishes_pure_and_mixed_pair_losses():
    nodes = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
    ]
    members = [(0, 1), (1, 1), (2, 2), (3, 2)]
    edges = [(0, 1, 1.0), (2, 3, 1.0)]
    pure = run_case(nodes, members, edges, dict.fromkeys(range(4), "A"))
    mixed = run_case(nodes, members, edges, {0: "A", 1: "A", 2: "B", 3: "B"})
    a = pure["objective"]["block_degeneracy"]["legal_internal_nodes"][0]["positive_boundary"]
    b = mixed["objective"]["block_degeneracy"]["legal_internal_nodes"][0]["positive_boundary"]
    assert a["complete_gain"] == b["complete_gain"] == 0.5
    assert a["complete_cross_actual"] == 0
    assert a["complete_cross_expected"] == 1
    assert a["actual_direct_lca_tp"] == a["local_complete_removed_tp"] == 4
    assert a["actual_direct_lca_fp"] == 0
    assert b["actual_direct_lca_tp"] == 0
    assert b["actual_direct_lca_fp"] == b["local_complete_removed_fp"] == 4
    comparison = compare_positive([pure, mixed])
    assert len(comparison["matched_cells"]) == 1
    assert comparison["matched_cells"][0]["pure"]["mean_complete_gain"] == 0.5
    assert comparison["matched_cells"][0]["mixed"]["mean_complete_gain"] == 0.5
    assert len(comparison["targets"]) == 2


def test_complete_frontier_gain_and_nested_losses_are_separate_from_direct_boundary():
    nodes = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
        dict(cluster_id=3, parent_id=1),
        dict(cluster_id=4, parent_id=1),
    ]
    members = [(0, 3), (1, 3), (2, 4), (3, 4), (4, 2), (5, 2)]
    case = run_case(
        nodes, members, [(0, 1, 1.0), (2, 3, 1.0), (4, 5, 1.0)], dict.fromkeys(range(6), "A")
    )
    rows = case["objective"]["block_degeneracy"]["legal_internal_nodes"]
    root = next(r for r in rows if r["cluster_id"] == 0)
    b = root["positive_boundary"]
    assert b["complete_gain"] == pytest.approx(2 / 3)
    assert root["direct_gain"] == pytest.approx(4 / 9)
    assert b["descendant_gain"] == pytest.approx(2 / 9)
    assert b["frontier_count"] == 3 and b["deeper_frontier"]
    assert b["actual_direct_lca_tp"] == 8
    assert b["local_complete_removed_tp"] == b["actual_subtree_removed_tp"] == 12
    comparison = compare_positive([case])
    assert comparison["strata"][0]["actual_direct_lca_tp"] == 12
    assert comparison["matched_cells"] == []
    copy = deepcopy(case)
    copy["reference_cohort"] = "clean_tp_gain"
    assert compare_positive([case, copy])["matched_cells"] == []


def test_random_fixed_cuts_keep_exact_boundary_and_loss_coverage():
    import random
    from collections import Counter

    rng = random.Random(42)
    nodes = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
        dict(cluster_id=3, parent_id=1),
        dict(cluster_id=4, parent_id=1),
    ]
    members = [(0, 3), (1, 3), (2, 4), (3, 4), (4, 2), (5, 2)]
    for _ in range(20):
        species = {p: rng.randrange(3) for p, _ in members}
        edges = [
            (u, v, rng.randrange(1, 5) / 4)
            for u in range(6)
            for v in range(u + 1, 6)
            if rng.random() < 0.5
        ]
        reference = {p: rng.randrange(2) for p in range(6)}
        cut = species_pair_cut(nodes, members, edges, species)
        result = profile(
            nodes,
            members,
            edges,
            species,
            set(cut["selected"]),
            context_flow=True,
            positive_boundary=True,
        )
        rows = result["legal_internal_nodes"]
        annotate_reference(rows, nodes, members, reference)
        annotate_positive(rows, nodes, members, reference)
        comparison = compare_positive(
            [
                dict(
                    component_id=0,
                    source_cluster_id=0,
                    reference_cohort="polluted",
                    objective=dict(block_degeneracy=result),
                )
            ]
        )
        ref_counts = Counter(reference.values())
        joint = Counter((reference[r["protein_id"]], r["cluster_id"]) for r in cut["members"])
        loss = sum(v * (v - 1) // 2 for v in ref_counts.values()) - sum(
            v * (v - 1) // 2 for v in joint.values()
        )
        assert sum(r["actual_direct_lca_tp"] for r in comparison["coverage"]) == loss
        for row in rows:
            if "positive_boundary" in row:
                b = row["positive_boundary"]
                assert b["complete_gain"] == pytest.approx(
                    row["direct_gain"] + b["descendant_gain"]
                )
                assert b["complete_boundary_identity_verified"]
