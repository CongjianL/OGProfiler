"""Graph objectives are checked against hand-worked cuts, not OF labels."""

import pytest

from benchmarks.og_extraction.fixed_tree_graph import graph_cut

NODES = [dict(cluster_id=0, parent_id=None, n_genes=4)] + [
    dict(cluster_id=i, parent_id=0, n_genes=2) for i in (1, 2)
]
MEMBERS = [(0, 1), (1, 1), (2, 2), (3, 2)]
EDGES = [(0, 1, 4.0), (2, 3, 4.0), (1, 2, 0.5)]


def test_weak_bridge_splits_but_strong_bridge_keeps_root():
    r = graph_cut(NODES, MEMBERS, EDGES, pair_penalty=1)
    assert r["selected"] == [1, 2]
    assert r["score"] == 6  # split: 8 - 2; root: 8.5 - 6 = 2.5
    assert len(r["members"]) == 4
    strong = graph_cut(NODES, MEMBERS, EDGES[:2] + [(1, 2, 5.0)], pair_penalty=1)
    assert strong["selected"] == [0] and strong["score"] == 7


def test_penalty_extremes_ties_and_terminal_indivisibility():
    assert graph_cut(NODES, MEMBERS, EDGES, pair_penalty=0)["selected"] == [0]
    large = graph_cut(NODES, MEMBERS, EDGES, pair_penalty=100)
    assert large["selected"] == [1, 2] and large["score"] == -192
    assert large["diagnostics"]["all_structural_leaves_selected"]
    tied = graph_cut(NODES, MEMBERS, EDGES[:2] + [(1, 2, 4.0)], pair_penalty=1)
    assert tied["selected"] == [0] and tied["score"] == 6


@pytest.mark.parametrize("penalty", [-1, float("nan"), float("inf"), True])
def test_invalid_penalty(penalty):
    with pytest.raises(ValueError, match="pair_penalty"):
        graph_cut(NODES, MEMBERS, EDGES, pair_penalty=penalty)


@pytest.mark.parametrize(
    "edges",
    [
        EDGES + [EDGES[0]],
        [(1, 0, 1)],
        [(0, 0, 1)],
        [(0, 9, 1)],
        [(0, 1, -1)],
        [(0, 1, float("nan"))],
    ],
)
def test_invalid_edges(edges):
    with pytest.raises(ValueError):
        graph_cut(NODES, MEMBERS, edges, pair_penalty=1)


def test_membership_and_failure_gates():
    with pytest.raises(ValueError, match="UNRESOLVED"):
        graph_cut(
            [dict(NODES[0], split_status="UNRESOLVED"), *NODES[1:]], MEMBERS, EDGES, pair_penalty=1
        )
    for members in [MEMBERS + [MEMBERS[0]], MEMBERS[:-1], [(0, 0), *MEMBERS[1:]]]:
        with pytest.raises(ValueError):
            graph_cut(NODES, members, EDGES, pair_penalty=1)


def test_order_invariance_and_no_event_qualification():
    r = graph_cut(NODES, MEMBERS, EDGES, pair_penalty=1)
    shuffled = graph_cut(
        [dict(n, event="III-3") for n in reversed(NODES)],
        list(reversed(MEMBERS)),
        list(reversed(EDGES)),
        pair_penalty=1,
    )
    assert r == shuffled
    singleton = graph_cut([dict(cluster_id=7, parent_id=None)], [(42, 7)], [], pair_penalty=1)
    assert singleton["members"] == [dict(protein_id=42, cluster_id=7)]


def test_multilevel_lca_and_polytomy():
    nodes = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
        dict(cluster_id=3, parent_id=1),
        dict(cluster_id=4, parent_id=1),
        dict(cluster_id=5, parent_id=1),
    ]
    r = graph_cut(
        nodes,
        [(0, 3), (1, 4), (2, 5), (3, 2)],
        [(0, 1, 3.0), (1, 2, 3.0), (0, 2, 3.0), (2, 3, 0.5)],
        pair_penalty=1,
    )
    assert r["selected"] == [1, 2]
    assert r["score"] == 6
    assert {x["cluster_id"]: x["internal_weight"] for x in r["node_scores"]} == {
        0: 9.5,
        1: 9.0,
        2: 0.0,
        3: 0.0,
        4: 0.0,
        5: 0.0,
    }


def test_cli_provenance_and_output_protection(tmp_path):
    import json
    import subprocess
    import sys

    source, out = tmp_path / "input.json", tmp_path / "out.json"
    source.write_text(json.dumps(dict(nodes=NODES, membership=MEMBERS, edges=EDGES)))
    command = [
        sys.executable,
        "-m",
        "benchmarks.og_extraction.fixed_tree_graph",
        "--input",
        str(source),
        "--out",
        str(out),
        "--pair-penalty",
        "1",
    ]
    subprocess.run(command, check=True, capture_output=True)
    report = json.loads(out.read_text())
    assert report["selected"] == [1, 2]
    assert not report["reference_labels_used"]
    assert len(report["input_sha256"]) == 64
    before = out.read_bytes()
    assert subprocess.run(command, capture_output=True).returncode != 0
    assert out.read_bytes() == before
    source.write_text(
        json.dumps(dict(nodes=NODES, membership=MEMBERS, edges=EDGES, reference={"p0": "OG0"}))
    )
    bad_command = list(command)
    bad_command[bad_command.index("--out") + 1] = str(tmp_path / "bad.json")
    failure = subprocess.run(bad_command, capture_output=True)
    assert failure.returncode != 0
    assert b"no reference labels" in failure.stderr
    assert not (tmp_path / "bad.json").exists()
