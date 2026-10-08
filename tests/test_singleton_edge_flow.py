import pytest

from benchmarks.og_extraction.singleton_edge_flow import trace_flow


def fixture():
    n = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
        dict(cluster_id=3, parent_id=1),
        dict(cluster_id=4, parent_id=1),
        dict(cluster_id=5, parent_id=1),
    ]
    m = [(0, 3), (1, 3), (2, 4), (3, 5), (4, 2)]
    s = {0: 0, 1: 1, 2: 1, 3: 1, 4: 0, 9: 0}
    edges = [(0, 1, 2.0), (0, 2, 0.5), (3, 4, 0.2)]
    incident = [(0, 2, 0.5), (3, 4, 0.2), (2, 9, 0.3)]
    return n, m, s, edges, incident


def test_raw_pair_weights_external_scope_and_tree_constraint():
    n, m, s, e, incident = fixture()
    r = trace_flow(n, m, e, incident, s, 1)
    assert r["direct_children"] == [3, 4, 5]
    assert r["singletons"] == {2: 4, 3: 5}
    assert r["incident_source_identity_verified"]
    assert {e["scope"] for e in r["incident_edges"]} == {
        "node_child",
        "outside_node_inside_source",
        "outside_source",
    }
    pair = next(p for p in r["child_pair_species_blocks"] if p["a"] == 3 and p["b"] == 4)
    assert pair["observed"] == 0.5
    assert any(p["observed"] == 0 and p["expected"] > 0 for p in r["child_pair_species_blocks"])
    assert [p["tree_representable"] for p in r["candidate_partitions"]] == [
        True,
        True,
        False,
        False,
        False,
    ]
    split = r["candidate_partitions"][1]
    assert split["gain_vs_keep"] == pytest.approx(
        sum(p["gain"] for p in r["child_pair_species_blocks"])
    )


def test_missing_source_incident_rejected():
    n, m, s, e, incident = fixture()
    with pytest.raises(ValueError, match="disagree"):
        trace_flow(n, m, e, incident[1:], s, 1)


def test_duplicate_raw_incident_rejected():
    n, m, s, e, incident = fixture()
    with pytest.raises(ValueError, match="incident"):
        trace_flow(n, m, e, incident + incident[:1], s, 1)
