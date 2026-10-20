from itertools import product

import pytest

from benchmarks.og_extraction.embleya_representation import node_profiles, strata
from benchmarks.og_extraction.tree_cut import optimal_cut, pair_f1_oracle

NODES = [
    dict(cluster_id=0, parent_id=None),
    dict(cluster_id=1, parent_id=0),
    dict(cluster_id=2, parent_id=0),
]


def test_keep_split_and_deterministic_tie():
    assert optimal_cut(NODES, {0: 3, 1: 1, 2: 1}).selected == (0,)
    assert optimal_cut(NODES, {0: 1, 1: 1, 2: 1}).selected == (1, 2)
    assert optimal_cut(NODES, {0: 2, 1: 1, 2: 1}).selected == (0,)
    assert optimal_cut(NODES, {0: 3, 1: 1, 2: 1}, eligible={1, 2}).selected == (1, 2)
    assert optimal_cut(NODES, {0: 3, 1: 1, 2: 1}, split_scores={0: 2}).selected == (1, 2)


def test_coverage_and_bad_inputs():
    with pytest.raises(ValueError, match="complete"):
        optimal_cut(NODES, {0: 0, 1: 0, 2: 0}, eligible={1})
    with pytest.raises(ValueError, match="UNRESOLVED"):
        optimal_cut([*NODES[:2], dict(NODES[2], split_status="UNRESOLVED")], {0: 100, 1: 0, 2: 0})
    with pytest.raises(ValueError, match="Non-finite"):
        optimal_cut(NODES, {0: float("nan"), 1: 0, 2: 0})
    with pytest.raises(ValueError, match="Disconnected"):
        optimal_cut(
            [
                dict(cluster_id=0, parent_id=None),
                dict(cluster_id=1, parent_id=2),
                dict(cluster_id=2, parent_id=1),
            ],
            {0: 0, 1: 0, 2: 0},
        )
    with pytest.raises(ValueError, match="exactly"):
        optimal_cut(NODES, {0: 0})


def test_global_f1_oracle_matches_exhaustive_forest_cuts():
    # Each component has four proteins. A: two clean pairs; B: one family of four.
    # A root merges reference groups; B root should be kept, independently of events.
    forest = [
        dict(nodes=NODES, tp={0: 2, 1: 1, 2: 1}, pp={0: 6, 1: 1, 2: 1}),
        dict(nodes=NODES, tp={0: 6, 1: 1, 2: 1}, pp={0: 6, 1: 1, 2: 1}),
    ]
    values = []
    for cuts in product([(0,), (1, 2)], repeat=2):
        tp = sum(sum(t["tp"][c] for c in cut) for t, cut in zip(forest, cuts, strict=True))
        pp = sum(sum(t["pp"][c] for c in cut) for t, cut in zip(forest, cuts, strict=True))
        values.append(2 * tp / (pp + 8))
    r = pair_f1_oracle(forest, 8)
    assert r["pair_f1"] == max(values) == 1
    assert r["cuts"] == [(1, 2), (0,)]


def test_oracle_cross_component_pairs_and_inseparable_terminal():
    tree = dict(nodes=[dict(cluster_id=0, parent_id=None)], tp={0: 1}, pp={0: 3})
    r = pair_f1_oracle([tree], 3)
    assert r["pair_f1"] == pytest.approx(1 / 3)
    with pytest.raises(ValueError, match="converge"):
        pair_f1_oracle([tree], 3, max_iterations=1)


def test_polytomy_candidate_excluded_by_legacy_but_represented():
    nodes = [
        dict(
            cluster_id=0,
            parent_id=None,
            component_id=0,
            depth=0,
            n_genes=3,
            n_species=3,
            child_count=3,
        )
    ]
    nodes += [
        dict(
            cluster_id=i,
            parent_id=0,
            component_id=0,
            depth=1,
            n_genes=1,
            n_species=1,
            child_count=0,
        )
        for i in (1, 2, 3)
    ]
    counts, full, eligible, events = node_profiles(
        nodes, [(i, i) for i in (1, 2, 3)], {i: "OF1" for i in (1, 2, 3)}, {1: 0, 2: 1, 3: 2}
    )
    assert counts[0] == {"OF1": 3} and full[0] == 3
    assert events[0] == "III-3" and 0 not in eligible
    assert eligible == {1, 2, 3}
    with pytest.raises(ValueError, match="Duplicate"):
        node_profiles(nodes, [(1, 1), (1, 2), (3, 3)], {1: "OF1"}, {1: 0, 3: 2})


def test_pair_strata_not_ortholog_accuracy():
    ref = {0: "a", 1: "a", 2: "a"}
    pred = {0: "x", 1: "x", 2: "y"}
    species = {0: 0, 1: 0, 2: 1}
    r = strata(ref, pred, species)
    assert r["same_species"]["recall"] == 1
    assert r["cross_species"]["recall"] == 0
    assert r["cross_species"]["reference_pairs"] == 2


def test_deep_tree_iterative_and_complete():
    nodes = []
    scores = {}
    parent = None
    for i in range(1100):
        c = 2 * i
        nodes.append(dict(cluster_id=c, parent_id=parent))
        scores[c] = -1
        nodes.append(dict(cluster_id=c + 1, parent_id=c))
        scores[c + 1] = 1
        parent = c
    nodes.append(dict(cluster_id=2200, parent_id=parent))
    scores[2200] = 1
    cut = optimal_cut(nodes, scores)
    assert len(cut.selected) == 1101


def test_end_to_end_fixed_tree_oracle_and_manifest_gate(tmp_path):
    import csv
    import json

    import pyarrow as pa
    import pyarrow.parquet as pq

    from benchmarks.og_extraction.embleya_representation import run_audit
    from ogprofiler.core.manifest import sha256_file

    run = tmp_path / "run"
    ref = tmp_path / "ref"
    ref.mkdir()
    for name in ("input", "components", "results", "hierarchy/components/component=00000000"):
        (run / name).mkdir(parents=True, exist_ok=True)

    def write(path, rows):
        pq.write_table(pa.Table.from_pylist(rows), path)

    write(
        run / "input/proteins.parquet",
        [dict(protein_id=i, species_id=i, original_id=f"p{i}", length=10) for i in range(4)],
    )
    write(
        run / "input/species.parquet",
        [dict(species_id=i, species_name=f"s{i}", source_file=f"s{i}.faa") for i in range(4)],
    )
    write(
        run / "components/index.parquet",
        [dict(protein_id=i, component_id=0 if i < 3 else 1) for i in range(4)],
    )
    (run / "run.yaml").write_text("{}\n")
    (ref / "Orthogroups.tsv").write_text("Orthogroup\ts0\ts1\ts2\ts3\nOG0\tp0\tp1\tp2\t\n")
    (ref / "Orthogroups_UnassignedGenes.tsv").write_text(
        "Orthogroup\ts0\ts1\ts2\ts3\nU0\t\t\t\tp3\n"
    )
    with (run / "results/members.tsv").open("w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["family_id", "protein_id", "species_id", "original_id"])
        for i in range(4):
            w.writerow([f"g{i}", i, i, f"p{i}"])
    folder = run / "hierarchy/components/component=00000000"
    nodes = [
        dict(
            cluster_id=0,
            parent_id=None,
            component_id=0,
            depth=0,
            n_genes=3,
            n_species=3,
            child_count=3,
            split_status="SPLIT",
        )
    ]
    nodes += [
        dict(
            cluster_id=i + 1,
            parent_id=0,
            component_id=0,
            depth=1,
            n_genes=1,
            n_species=1,
            child_count=0,
            split_status="TERMINAL",
        )
        for i in range(3)
    ]
    write(folder / "nodes.parquet", nodes)
    write(
        folder / "members.parquet",
        [dict(protein_id=i, terminal_cluster_id=i + 1) for i in range(3)],
    )
    manifest = dict(
        hierarchy_status="RESOLVED",
        structural_validation_passed=True,
        output_checksums={n: sha256_file(folder / n) for n in ["nodes.parquet", "members.parquet"]},
        input_checksums={"input/proteins.parquet": sha256_file(run / "input/proteins.parquet")},
    )
    (folder / "hierarchy-manifest.json").write_text(json.dumps(manifest))
    out = tmp_path / "audit"
    run_audit(run, ref, out)
    report = json.loads((out / "report.json").read_text())
    assert report["eligible_oracle_cut"]["pair_f1"] == 0
    assert report["eligible_oracle_strata"]["cross_species"]["recall"] == 0
    assert report["cut_family_comparison"]["unrestricted_minus_eligible_pair_f1"] == 1
    assert len(json.loads((out / "eligible-oracle-cut-DIAGNOSTIC-ONLY.json").read_text())) == 4
    assert report["oracle_cut"]["pair_f1"] == 1
    assert report["actual"]["pair_recall"] == 0
    assert report["positive_eligibility_gap"] == 1
    assert report["whole_hierarchy_macro_best_f1"] == 1
    assert report["eligible_macro_best_f1"] == 0.5
    assert len(json.loads((out / "oracle-cut-DIAGNOSTIC-ONLY.json").read_text())) == 2
    (folder / "members.parquet").write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="Artifact mismatch"):
        run_audit(run, ref, tmp_path / "corrupt-audit")


def test_eligible_oracle_matches_exhaustive_feasible_cuts():
    forest = [
        dict(nodes=NODES, tp={0: 2, 1: 1, 2: 1}, pp={0: 6, 1: 1, 2: 1}, eligible={0, 1, 2}),
        dict(nodes=NODES, tp={0: 6, 1: 1, 2: 1}, pp={0: 6, 1: 1, 2: 1}, eligible={1, 2}),
    ]
    # Only the second root is disallowed; enumerate the two feasible forest cuts.
    values = []
    for cut in [(0,), (1, 2)]:
        tp = sum(forest[0]["tp"][c] for c in cut) + 2
        pp = sum(forest[0]["pp"][c] for c in cut) + 2
        values.append(2 * tp / (pp + 8))
    result = pair_f1_oracle(forest, 8, restrict_eligible=True)
    assert result["cuts"] == [(1, 2), (1, 2)]
    assert result["pair_f1"] == max(values) == pytest.approx(2 / 3)
    assert pair_f1_oracle(forest, 8)["pair_f1"] == 1


def test_eligible_oracle_rejects_missing_or_infeasible_eligibility():
    tree = dict(nodes=NODES, tp={0: 2, 1: 1, 2: 1}, pp={0: 6, 1: 1, 2: 1})
    with pytest.raises(ValueError, match="explicit eligible"):
        pair_f1_oracle([tree], 2, restrict_eligible=True)
    for allowed in [set(), {1}]:
        with pytest.raises(ValueError, match="complete"):
            pair_f1_oracle([dict(tree, eligible=allowed)], 2, restrict_eligible=True)
    with pytest.raises(ValueError, match="Unknown eligible"):
        pair_f1_oracle([dict(tree, eligible={99})], 2, restrict_eligible=True)
