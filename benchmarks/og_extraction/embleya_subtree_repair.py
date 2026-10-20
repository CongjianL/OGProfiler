"""Frozen component-0 diagnosis and bounded six-subtree repair experiment.

OF labels select three diagnosed cases and score outputs, never select gamma.
Control selection uses only node topology/size. This is not whole-dataset scoring.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from dataclasses import asdict, replace
from itertools import combinations
from pathlib import Path
from unittest.mock import patch

import pyarrow.parquet as pq

from benchmarks.og_extraction.embleya_fallback_search import PROTOCOL, search_reserved_fallback
from benchmarks.og_extraction.embleya_reference import compare, read_of
from benchmarks.og_extraction.refog_split_audit import load_component
from ogprofiler.config import hierarchy_config, load_config
from ogprofiler.core.manifest import sha256_file
from ogprofiler.graph.components import Component
from ogprofiler.graph.partition import load_component_edge_table
from ogprofiler.hierarchy.engine import infer_component_hierarchy
from ogprofiler.hierarchy.resolution import search_resolution
from ogprofiler.hierarchy.validation import validate_hierarchy
from ogprofiler.orthogroups.engine import extract_component_orthogroups


def experimental_search(graph, config, **kwargs):
    if not 3 <= graph.vcount() < 10000:
        return search_resolution(graph, config, **kwargs)

    def evaluate(gamma):
        fixed = replace(config, gamma_min=gamma, gamma_max=gamma, topology_policy="kway_v1")
        result = search_resolution(graph, fixed, **kwargs)
        assert len(result.candidates) == 1
        return result.candidates[0]

    return search_reserved_fallback(config, evaluate)


def clades(nodes, members):
    by_id = {n["cluster_id"]: n for n in nodes}
    leaves = defaultdict(set)
    for protein, leaf in members:
        while leaf is not None:
            leaves[leaf].add(protein)
            leaf = by_id[leaf]["parent_id"]
    return leaves


def diagnose(run, reference_dir):
    proteins = pq.read_table(run / "input/proteins.parquet").to_pylist()
    species = {
        r["species_name"]: r["species_id"]
        for r in pq.read_table(run / "input/species.parquet").to_pylist()
    }
    keys = {(r["species_id"], r["original_id"]): r["protein_id"] for r in proteins}
    reference = read_of(reference_dir / "Orthogroups.tsv", species, keys)
    labels = {keys[k]: v for k, v in reference.items()}
    targets = {
        og: sorted(p for p, r in labels.items() if r == og)
        for og in ("OG0000000", "OG0000001", "OG0000002")
    }
    view = load_component(run, 0, sum(targets.values(), []))
    groups = {
        int(r["protein_id"]): r["family_id"]
        for r in csv.DictReader((run / "results/members.tsv").open(), delimiter="\t")
    }
    report = {}
    for og, ids in targets.items():
        losses, reasons = Counter(), Counter()
        for a, b in combinations(ids, 2):
            if groups[a] == groups[b]:
                continue
            common = []
            for x, y in zip(view["paths"][a], view["paths"][b], strict=False):
                if x != y:
                    break
                common.append(x)
            losses[common[-1]] += 1
            eligible = [
                c
                for c in common
                if view["events"][c]["v1_event"] is None
                or (view["events"][c]["v1_event"] == "I" and view["nodes"][c]["n_species"] > 1)
            ]
            reasons["COMMON_ELIGIBLE" if eligible else "NO_COMMON_ELIGIBLE"] += 1
        report[og] = dict(
            lost_pairs=sum(losses.values()),
            reasons=dict(reasons),
            nodes=[
                dict(
                    view["nodes"][c],
                    lost_pairs=n,
                    event=view["events"][c]["v1_event"],
                    candidate=view["selected"].get(c),
                )
                for c, n in losses.most_common()
            ],
        )
    return report, proteins, labels


def run_experiment(run, reference_dir, out):
    out.mkdir(parents=True, exist_ok=False)
    config = hierarchy_config(load_config(str(run / "run.yaml"))["hierarchy"])
    assert config.resolution.topology_policy == "soft_binary_24_v2"
    assert config.max_depth == 42 and config.leiden_iterations == 10
    diagnosis, proteins, labels = diagnose(run, reference_dir)
    (out / "diagnosis.json").write_text(json.dumps(diagnosis, indent=2))
    folder = run / "hierarchy/components/component=00000000"
    manifest = json.loads((folder / "hierarchy-manifest.json").read_text())
    import igraph
    import leidenalg

    assert manifest["parameters"]["environment"] == {
        "igraph": igraph.__version__,
        "leidenalg": leidenalg.__version__,
    }
    checks = {run / k: v for k, v in manifest["input_checksums"].items()}
    checks.update({folder / k: v for k, v in manifest["output_checksums"].items()})
    for path, digest in checks.items():
        assert sha256_file(path) == digest
    nodes = pq.ParquetFile(folder / "nodes.parquet").read().to_pylist()
    members = pq.ParquetFile(folder / "members.parquet").read().to_pylist()
    by_id = {n["cluster_id"]: n for n in nodes}
    # Component is small enough for ID-only ancestor membership indexing; no graphs per node.
    sets = clades(nodes, [(r["protein_id"], r["terminal_cluster_id"]) for r in members])
    panel = {d["nodes"][0]["cluster_id"]: "diagnosed" for d in diagnosis.values()}
    occupied = set().union(*(sets[c] for c in panel))
    for node in nodes:
        c = node["cluster_id"]
        if (
            100 <= node["n_genes"] <= 200
            and node["selection_kind"] == "BINARY"
            and not sets[c] & occupied
        ):
            panel[c] = "control"
            occupied.update(sets[c])
            if len(panel) == 6:
                break
    assert len(panel) == 6
    table = load_component_edge_table(run / "components", 0)
    graph, global_ids = table.to_igraph()
    indices = {p: i for i, p in enumerate(global_ids)}
    species = {p["protein_id"]: p["species_id"] for p in proteins}
    original = {p["protein_id"]: p["original_id"] for p in proteins}
    results = []
    for cluster, role in panel.items():
        ids = tuple(sorted(sets[cluster]))
        component = Component(0, ids, ())
        small = graph.induced_subgraph([indices[p] for p in ids])
        assert tuple(small.vs["protein_id"]) == ids
        local_config = replace(config, max_depth=config.max_depth - by_id[cluster]["depth"])
        result_row = dict(
            source_cluster=cluster,
            role=role,
            proteins=len(ids),
            remaining_depth=local_config.max_depth,
            variants={},
        )
        for name, search in [("baseline", search_resolution), ("refined", experimental_search)]:
            with patch("ogprofiler.hierarchy.engine.search_resolution", search):
                result = infer_component_hierarchy(component, local_config, species, small)
                repeated = infer_component_hierarchy(component, local_config, species, small)
            validate_hierarchy(component, result)
            assert result.nodes == repeated.nodes
            assert result.terminal_membership == repeated.terminal_membership
            assert result.resolution_candidates == repeated.resolution_candidates
            counts = Counter(c.cluster_id for c in result.resolution_candidates)
            assert max(counts.values(), default=0) <= 24
            assert result.metrics.leiden_calls == 3 * len(result.resolution_candidates)
            assert all(
                c.valid and not c.original_violations
                for c in result.resolution_candidates
                if c.selected
            )
            new_nodes = [asdict(n) for n in result.nodes]
            if name == "baseline":
                actual = clades(new_nodes, result.terminal_membership)
                old = {frozenset(s): by_id[c] for c, s in sets.items() if s <= sets[cluster]}
                assert len(actual) == len(old)
                for n in new_nodes:
                    prev = old[frozenset(actual[n["cluster_id"]])]
                    for field in (
                        "n_genes",
                        "n_species",
                        "child_count",
                        "resolution",
                        "split_status",
                        "terminal_reason",
                        "selection_phase",
                        "selection_kind",
                    ):
                        assert n[field] == prev[field], (cluster, field, n[field], prev[field])
                result_row["frozen_baseline_clades_match"] = True
            unresolved = [n for n in new_nodes if n["split_status"] == "UNRESOLVED"]
            detail = dict(
                nodes=new_nodes,
                candidates=[asdict(c) for c in result.resolution_candidates],
                membership=result.terminal_membership,
                metrics=asdict(result.metrics),
                repeat_identical=True,
                unresolved=len(unresolved),
            )
            if not unresolved:
                og = extract_component_orthogroups(
                    new_nodes,
                    [
                        dict(protein_id=p, terminal_cluster_id=c)
                        for p, c in result.terminal_membership
                    ],
                    species,
                    original,
                    total_species=len(set(species.values())),
                )
                prediction = {p: g.local_group_id for g in og.groups for p in g.protein_ids}
                assert set(prediction) == set(ids) and not og.unassigned
                ref = {p: labels[p] for p in ids if p in labels}
                detail["local_reference_score"] = compare(ref, prediction, dict.fromkeys(ids, 0))
                detail["groups"] = [list(g.protein_ids) for g in og.groups]
            (out / f"{cluster}-{name}.json").write_text(json.dumps(detail, indent=2))
            result_row["variants"][name] = {
                k: v
                for k, v in detail.items()
                if k not in ("nodes", "candidates", "membership", "groups")
            }
        results.append(result_row)
        (out / "progress.json").write_text(json.dumps(results, indent=2))
    assert all(sha256_file(p) == h for p, h in checks.items())
    resolved = all(v["unresolved"] == 0 for r in results for v in r["variants"].values())
    non_regression, strict_gain = False, False
    if resolved:
        scores = [
            (
                r["variants"]["baseline"]["local_reference_score"],
                r["variants"]["refined"]["local_reference_score"],
            )
            for r in results
        ]
        non_regression = all(
            b["true_positive_pairs"] >= a["true_positive_pairs"]
            and b["discordant_pairs"] <= a["discordant_pairs"]
            for a, b in scores
        )
        strict_gain = any(
            b["true_positive_pairs"] > a["true_positive_pairs"]
            or b["discordant_pairs"] < a["discordant_pairs"]
            for a, b in scores
        )
    report = dict(
        protocol=PROTOCOL,
        completed=True,
        immutable_inputs_verified=True,
        production_changed=False,
        whole_dataset_score=False,
        engineering_checks_passed=True,
        all_subtrees_resolved=resolved,
        panel_non_regression=non_regression,
        panel_strict_gain=strict_gain,
        repair_supported_on_panel=resolved and non_regression and strict_gain,
        results=results,
        source_run=str(run),
        input_hashes={str(p): h for p, h in checks.items()},
        reference_hash=sha256_file(reference_dir / "Orthogroups.tsv"),
    )
    (out / "report.json").write_text(json.dumps(report, indent=2))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--reference-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    run_experiment(args.run, args.reference_dir, args.out)


if __name__ == "__main__":
    main()
