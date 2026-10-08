from test_frozen_graph_batch import frozen  # noqa: F401

from benchmarks.og_extraction.frozen_graph_batch import run_batch
from benchmarks.og_extraction.graph_merge_diagnosis import boundary_evidence, diagnose, transitions


def test_exact_transitions_not_net_difference():
    r = transitions(
        [0, 1, 2, 3, 4], {0: "a", 1: "a", 2: "b", 3: "b"}, {0: "x", 1: "y", 2: "x", 3: "z", 4: "z"}
    )
    assert r == dict(assigned_size=4, tp=2, fp=4, new_tp=2, new_fp=3, retained_fp=1)


def test_single_bridge_vs_broad_connections():
    branches = {0: 0, 1: 0, 2: 1, 3: 1}
    a = boundary_evidence(branches, [(0, 2, 4)])
    b = boundary_evidence(branches, [(0, 2, 1), (0, 3, 1), (1, 2, 1), (1, 3, 1)])
    assert a["boundary_weight"] == b["boundary_weight"] == 4
    assert a["endpoint_coverage"] == 0.5 and b["endpoint_coverage"] == 1
    assert a["largest_edge_share"] == 1 and b["largest_edge_share"] == 0.25
    assert a["endpoint_weight_hhi"] == 0.5 and b["endpoint_weight_hhi"] == 0.25
    assert boundary_evidence(branches, [])["largest_edge_share"] is None


def test_batch_diagnosis_end_to_end(frozen, tmp_path):  # noqa: F811
    run, ref, _ = frozen
    (run / "results").mkdir()
    (run / "results/members.tsv").write_text(
        "family_id\tprotein_id\tspecies_id\toriginal_id\na\t0\t0\tp0\nb\t1\t1\tp1\nc\t2\t2\tp2\n"
    )
    batch = tmp_path / "batch"
    run_batch(run, ref, batch, [0.1])
    r = diagnose(run, ref, batch, tmp_path / "diagnosis")
    assert r["totals"]["new_tp"] == 1 and r["totals"]["new_fp"] == 0
    assert r["clean_gain_controls"] == 1
    assert next(x for x in r["rows"] if x["new_tp"])["boundary"]["endpoint_coverage"] == 1


def test_relative_density_and_child_strength_hand_calculation():
    from benchmarks.og_extraction.graph_merge_diagnosis import relative_strength

    r = relative_strength({0: 0, 1: 0, 2: 1, 3: 1}, [(0, 1, 4.0), (2, 3, 4.0), (1, 2, 0.5)])
    assert r["within_possible_pairs"] == 2 and r["cross_possible_pairs"] == 4
    assert r["within_density"] == 4 and r["cross_density"] == 0.125
    assert r["cross_to_within_density"] == 0.03125
    assert r["child_strengths"][0]["outgoing_strength_fraction"] == 0.5 / 8.5
    assert r["child_pair_density_cv"] == 0


def test_relative_undefined_and_missing_pairs_are_not_smoothed():
    from benchmarks.og_extraction.graph_merge_diagnosis import relative_strength

    r = relative_strength({0: 0, 1: 1, 2: 2}, [(0, 1, 3.0)])
    assert r["within_status"] == "no_possible_pairs"
    assert r["cross_to_within_density"] is None
    assert r["cross_density"] == 1 and r["child_pair_density_mean"] == 1
    assert r["zero_weight_child_pairs"] == 2
    assert r["child_strengths"][2]["outgoing_strength_fraction"] is None
    r = relative_strength({0: 0, 1: 0, 2: 1}, [(1, 2, 1.0)])
    assert r["within_status"] == "zero_weight" and r["cross_to_within_density"] is None
