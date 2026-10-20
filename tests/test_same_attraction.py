from collections import Counter
from fractions import Fraction
from itertools import product

import pytest

from benchmarks.og_extraction.same_attraction import intervene
from benchmarks.og_extraction.species_pair_graph import species_pair_cut
from benchmarks.og_extraction.tree_cut import topology
from benchmarks.qfo.intervention_metrics import annotate_intervention, metrics, summarize


def fixture(deep=False):
    nodes = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
    ]
    if deep:
        nodes += [dict(cluster_id=3, parent_id=1), dict(cluster_id=4, parent_id=1)]
        members = [(0, 3), (1, 3), (2, 4), (3, 4), (4, 2), (5, 2)]
    else:
        members = [(0, 1), (1, 1), (2, 2), (3, 2)]
    return nodes, members


def run(nodes, members, edges, species):
    baseline = species_pair_cut(nodes, members, edges, species)
    return intervene(nodes, members, edges, species, baseline["selected"]), baseline


def test_same_only_keeps_root_with_nonzero_baseline_not_per_group_shift():
    nodes, m = fixture()
    result, base = run(nodes, m, [(0, 1, 1.0), (2, 3, 1.0)], dict.fromkeys(range(4), 0))
    assert base["selected"] == [1, 2]
    assert result["candidate_selected"] == [0]
    assert result["root_baseline_exact"] == "1"
    assert result["candidate_score"] == 1
    assert result["candidate_gain_over_merge"] == 0
    assert result["baseline_gain_candidate"] == 0
    annotate_intervention(result, nodes, m, dict.fromkeys(range(4), "A"))
    assert result["pair_transitions"]["recovered_tp"] == 4
    assert result["pair_transitions"]["new_merge_fp"] == 0


def test_cross_only_exactly_reproduces_original_score_cut_and_zero_root():
    nodes, m = fixture()
    result, base = run(nodes, m, [(0, 1, 1.0), (2, 3, 1.0)], {0: 0, 1: 1, 2: 0, 3: 1})
    assert result["candidate_selected"] == base["selected"]
    assert result["candidate_gain_over_merge"] == base["score"] == 0.5
    assert result["root_baseline_exact"] == "0"
    assert all(not r["local_decision_changed"] for r in result["node_transitions"])


def test_zero_graph_and_invalid_frozen_baseline():
    nodes, m = fixture()
    result, _ = run(nodes, m, [], dict.fromkeys(range(4), 0))
    assert result["candidate_selected"] == [0]
    assert result["root_baseline_exact"] == "0"
    assert result["candidate_gain_over_merge"] == 0
    with pytest.raises(ValueError, match="frozen exact cut"):
        intervene(nodes, m, [], dict.fromkeys(range(4), 0), [1, 2])


def enumerate_cuts(c, children):
    cuts = [(c,)]
    if children[c]:
        cuts += [
            sum(parts, ())
            for parts in product(*(enumerate_cuts(ch, children) for ch in children[c]))
        ]
    return cuts


def brute_score(nodes, members, edges, species, cut):
    by_id, _, order = topology(nodes)
    owner = {}
    for c in order:
        owner[c] = c if c in cut else owner.get(by_id[c]["parent_id"])
    pred = {p: owner[leaf] for p, leaf in members}
    weights, k = Counter(), Counter()
    within = Fraction()
    for u, v, w in edges:
        w = Fraction.from_float(w)
        s, t = species[u], species[v]
        weights[min(s, t), max(s, t)] += w
        k[pred[u], s, t] += w
        k[pred[v], t, s] += w
        if pred[u] == pred[v]:
            within += w
    expected = sum(
        (k[g, s, t] * k[g, t, s] / w for (s, t), w in weights.items() if s != t and w for g in cut),
        Fraction(),
    )
    total = sum(weights.values(), Fraction())
    return (within - expected) / total if total else Fraction()


def test_random_candidate_matches_exhaustive_legal_cut_and_unique_lca_transitions():
    import random

    rng = random.Random(42)
    nodes, m = fixture(True)
    _, children, order = topology(nodes)
    cuts = enumerate_cuts(order[0], children)
    for _ in range(30):
        species = {p: rng.randrange(3) for p, _ in m}
        edges = [
            (u, v, rng.randrange(1, 5) / 4)
            for u in range(6)
            for v in range(u + 1, 6)
            if rng.random() < 0.5
        ]
        r, _ = run(nodes, m, edges, species)
        best = max(brute_score(nodes, m, edges, species, c) for c in cuts)
        assert r["candidate_score"] == float(best)
        assert brute_score(nodes, m, edges, species, r["candidate_selected"]) == best
        ref = {p: rng.randrange(3) for p, _ in m}
        annotate_intervention(r, nodes, m, ref)
        a = metrics(ref, r["baseline_prediction"])
        b = metrics(ref, r["candidate_prediction"])
        t = Counter(r["pair_transitions"])
        assert t["new_split_tp"] == t["removed_fp"] == 0
        assert r["coarsening_verified"]
        assert b["tp"] - a["tp"] == t["recovered_tp"] - t["new_split_tp"]
        assert b["fp"] - a["fp"] == t["new_merge_fp"] - t["removed_fp"]


def test_full_hybrid_keeps_explicit_outside_namespace_and_reference_universe():
    nodes, m = fixture()
    r, _ = run(nodes, m, [(0, 1, 1.0), (2, 3, 1.0)], dict.fromkeys(range(4), 0))
    reference = {0: "A", 1: "A", 2: "B", 3: "B", 4: "A", 5: "B"}
    annotate_intervention(r, nodes, m, reference)
    summary = summarize(
        [
            dict(
                component_id=0,
                source_cluster_id=0,
                reference_cohort="polluted",
                objective=dict(same_attraction_intervention=r),
            )
        ],
        reference,
        dict.fromkeys(range(7), "X"),
    )
    a = summary["full_reference_universe"]["baseline"]
    b = summary["full_reference_universe"]["candidate"]
    assert summary["total_proteins"] == 7 and summary["target_proteins"] == 4
    assert a["assigned_proteins"] == b["assigned_proteins"] == 6
    assert a["reference_pairs"] == b["reference_pairs"] == 6
    assert a["tp"] == b["tp"] == 2 and a["fp"] == 1 and b["fp"] == 5
    assert summary["pair_transitions"]["new_merge_fp"] == 4
    assert summary["full_reference_universe"]["frozen_production"]["predicted_pairs"] == 15
    assert summary["target_reference_universe"]["baseline"]["reference_pairs"] == 2
    assert metrics({}, {})["pair_f1"] == 0


def test_ancestor_consumption_change_is_not_a_new_local_split():
    nodes, m = fixture(True)
    r, _ = run(nodes, m, [(0, 1, 1.0), (2, 3, 1.0), (4, 5, 1.0)], dict.fromkeys(range(6), 0))
    root = next(x for x in r["node_transitions"] if x["cluster_id"] == 0)
    child = next(x for x in r["node_transitions"] if x["cluster_id"] == 1)
    assert root["baseline_state"] == "active_split" and root["candidate_state"] == "active_keep"
    assert child["baseline_state"] == "active_split" and child["candidate_state"] == "inactive"
    assert child["candidate_selected_ancestor"] == 0
    assert child["root_path"] == [0, 1]
    assert all(
        not x["baseline_local_keep"] or x["candidate_local_keep"] for x in r["node_transitions"]
    )


def test_frozen_batch_adapter_runs_candidate_and_scores_full_reference_universe(tmp_path):
    import json

    import pyarrow as pa
    import pyarrow.parquet as pq

    from benchmarks.qfo.internal_boundary import diagnose
    from ogprofiler.core.manifest import sha256_file

    run_dir = tmp_path / "run"
    trace = tmp_path / "trace"
    replay = tmp_path / "replay"
    ref = tmp_path / "reference"
    for path in (trace, replay, ref, run_dir / "results"):
        path.mkdir(parents=True)
    nodes, m = fixture()

    def parquet(path, rows):
        path.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.Table.from_pylist(rows), path)

    parquet(
        run_dir / "input/proteins.parquet",
        [dict(protein_id=p, species_id=0, original_id=f"p{p}") for p in range(4)],
    )
    parquet(run_dir / "input/species.parquet", [dict(species_id=0, species_name="S")])
    folder = run_dir / "hierarchy/components/component=00000000"
    parquet(folder / "nodes.parquet", nodes)
    parquet(folder / "members.parquet", [dict(protein_id=p, terminal_cluster_id=c) for p, c in m])
    parquet(
        run_dir / "components/edges/component=00000000/part.parquet",
        [dict(u=0, v=1, weight=1.0), dict(u=2, v=3, weight=1.0)],
    )
    (run_dir / "results/members.tsv").write_text("protein_id\tfamily_id\n0\tX\n1\tX\n2\tY\n3\tY\n")
    (ref / "Orthogroups.tsv").write_text("Orthogroup\tS\nA\tp0, p1, p2, p3\n")
    (trace / "summary.json").write_text(
        json.dumps(dict(replay_verified=True, immutable_inputs_verified=True))
    )
    (trace / "input-hashes.json").write_text(
        json.dumps({str(ref / "Orthogroups.tsv"): sha256_file(ref / "Orthogroups.tsv")})
    )
    (trace / "groups.json").write_text(
        json.dumps(
            [
                dict(
                    scope="complete_source",
                    msa_admitted=True,
                    source_full_size=4,
                    component_id=0,
                    source_cluster_id=0,
                    original_event="III-1",
                    output_transitions=dict(tp=6, fp=0, new_tp=4, new_fp=0),
                )
            ]
        )
    )
    (replay / "summary.json").write_text(
        json.dumps(
            dict(
                context_flow=True,
                exact_conditioned=True,
                coverage_verified=True,
                immutable_inputs_verified=True,
            )
        )
    )
    (replay / "input-hashes.json").write_text("{}")
    (replay / "cuts-DIAGNOSTIC-ONLY.json").write_text(
        json.dumps([dict(component_id=0, source_cluster_id=0, prediction={0: 1, 1: 1, 2: 2, 3: 2})])
    )
    out = tmp_path / "output"
    diagnose(
        run_dir,
        ref,
        trace,
        out,
        condition_species_pairs=True,
        exact_conditioned=True,
        block_degeneracy=True,
        context_flow=True,
        endpoint_replay_dir=replay,
        same_attraction_intervention=True,
    )
    summary = json.loads((out / "summary.json").read_text())
    assert summary["endpoint_replay_verified"] and not summary["production_changed"]
    result = summary["intervention_summary"]
    assert result["full_reference_universe"]["baseline"]["pair_f1"] == 0.5
    assert result["full_reference_universe"]["candidate"]["pair_f1"] == 1
    assert result["pair_transitions"] == dict(
        recovered_tp=4, new_split_tp=0, new_merge_fp=0, removed_fp=0
    )
    detail = json.loads((out / "cases.json").read_text())[0]["objective"][
        "same_attraction_intervention"
    ]
    assert detail["root_baseline_exact"] == "1" and detail["coarsening_verified"]
