"""H4/H5 frozen-SSN diagnostics; no changes to OG decisions or official scorer."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import pyarrow.parquet as pq
import yaml

from benchmarks.og_extraction.orthobench import read_groups, score
from ogprofiler.config import DEFAULT_CONFIG, deep_merge, validate_config
from ogprofiler.core.manifest import sha256_file


def validate_h4(run: Path, parallel: Path):
    """Gate H5 on complete component0 coverage, policy and deterministic replay."""
    from ogprofiler.config import hierarchy_config, load_config
    from ogprofiler.hierarchy.stage import hierarchy_component_is_verified

    config = hierarchy_config(load_config(str(run / "run.yaml"))["hierarchy"])
    summary = summarize_hierarchy(run, (0,))
    folder = run / "hierarchy/components/component=00000000"
    other = parallel / "hierarchy/components/component=00000000"
    tables = {
        name: pq.ParquetFile(folder / f"{name}.parquet").read()
        for name in ("nodes", "members", "candidates")
    }
    nodes = tables["nodes"].to_pylist()
    members = tables["members"].to_pylist()
    candidates = tables["candidates"].to_pylist()
    expected = pq.read_table(
        run / "components/index.parquet",
        columns=["protein_id"],
        filters=[("component_id", "=", 0)],
    )["protein_id"].to_pylist()
    ids = [row["protein_id"] for row in members]
    leaf_sizes = Counter(row["terminal_cluster_id"] for row in members)
    leaves = {n["cluster_id"]: n["n_genes"] for n in nodes if n["child_count"] == 0}
    parallel_config = hierarchy_config(
        load_config(str(parallel / "run.yaml"), ["hierarchy.subtree_workers=2"])["hierarchy"]
    )
    checks = dict(
        resolved=summary["resolved"],
        complete_unique_members=len(ids) == len(set(ids)) and set(ids) == set(expected),
        leaf_sizes_match=dict(leaf_sizes) == leaves,
        candidate_budget=summary["candidate_budget_passed"],
        robust_call_count=summary["leiden_calls"] == 3 * len(candidates),
        selected_policy=all(
            c["valid"]
            and not c["violations"]
            and c["stability_evaluated"]
            and c["stability"] >= config.resolution.stability_threshold
            for c in candidates
            if c["selected"]
        ),
        serial_manifest=hierarchy_component_is_verified(run, 0, config),
        parallel_manifest=hierarchy_component_is_verified(parallel, 0, parallel_config),
        serial_parallel_equal=all(
            table.equals(pq.ParquetFile(other / f"{name}.parquet").read())
            for name, table in tables.items()
        ),
    )
    return dict(passed=all(checks.values()), checks=checks, hierarchy=summary)


def migrate_config(previous):
    config = deep_merge(DEFAULT_CONFIG, previous)
    config["hierarchy"].update(
        admission_policy="nonempty_children_v1",
        resolution_strategy="bounded_adaptive_v2",
        min_family_size=None,
        recursion_stop_size=1,
        max_candidate_evaluations=24,
        max_coarse_candidates=10,
        rescue_grid_points=8,
        component_leiden_call_budget=None,
        leiden_iterations=10,
    )
    validate_config(config)
    return config


def summarize_hierarchy(run: Path, components=None):
    folders = sorted((run / "hierarchy/components").glob("component=*"))
    if components is not None:
        folders = [p for p in folders if int(p.name.split("=")[1]) in components]
    counts, reasons, sizes = Counter(), Counter(), Counter()
    roots, unresolved_proteins, calls, budget_passed = [], 0, 0, True
    for folder in folders:
        nodes = pq.ParquetFile(folder / "nodes.parquet").read().to_pylist()
        candidates = pq.ParquetFile(folder / "candidates.parquet").read().to_pylist()
        by_node = Counter(c["cluster_id"] for c in candidates)
        budget_passed &= all(by_node[c["cluster_id"]] <= c["evaluation_budget"] for c in candidates)
        metrics = json.loads((folder / "metrics.json").read_text())
        calls += metrics["leiden_calls"]
        for node in nodes:
            counts[node["split_status"]] += 1
            if node["terminal_reason"]:
                reasons[node["terminal_reason"]] += 1
                sizes[node["n_genes"]] += 1
            if node["split_status"] == "UNRESOLVED":
                unresolved_proteins += node["n_genes"]
        root = next(n for n in nodes if n["parent_id"] is None)
        roots.append(
            dict(
                root,
                component_id=int(folder.name.split("=")[1]),
                max_depth=max(n["depth"] for n in nodes),
                selected_candidates=[
                    c for c in candidates if c["cluster_id"] == root["cluster_id"] and c["selected"]
                ],
                root_candidates=[c for c in candidates if c["cluster_id"] == root["cluster_id"]]
                if int(folder.name.split("=")[1]) == 0
                else [],
            )
        )
    return dict(
        components=len(folders),
        node_statuses=dict(counts),
        terminal_reasons=dict(reasons),
        structural_leaf_size_histogram=dict(sizes),
        roots=roots,
        root_split=bool(roots) and all(r["split_status"] == "SPLIT" for r in roots),
        unresolved_nodes=counts["UNRESOLVED"],
        unresolved_proteins=unresolved_proteins,
        resolved=bool(folders) and not counts["UNRESOLVED"],
        leiden_calls=calls,
        candidate_budget_passed=budget_passed,
        accuracy_improvement_claimed=False,
    )


def score_campaign(root: Path, origin: Path):
    out = root / "h5-metrics"
    out.mkdir(exist_ok=False)
    results = []
    for label, run in [
        ("previous_repaired_mean", origin / "repaired-mean"),
        ("new_bounded_hierarchy", root / "new-hierarchy"),
    ]:
        for strategy, path in [
            ("terminal", "results/terminal-families/members.tsv"),
            ("v1_compatible", "results/members.tsv"),
        ]:
            results.append(
                score(
                    read_groups(run / path, exchange=True),
                    root / "benchmark",
                    out / f"{label}_{strategy}",
                    f"{label}_{strategy}",
                )
            )
    v1_manifest = dict(
        line.split("\t", 1)
        for line in (origin / "historical-v1-manifest.tsv").read_text().splitlines()
    )
    assert (
        v1_manifest["dataset_digest"]
        == "380c85d9548c607f5df6daf656d74b21597fe215ef78d4f8ff296776bb1fdd07"
    )
    results.append(
        score(
            read_groups(origin / "historical-v1-groups.tsv"),
            root / "benchmark",
            out / "historical_full_v1",
            "historical_full_v1",
        )
    )
    keys = (
        "official_precision",
        "official_recall",
        "official_f1",
        "coverage",
        "refog_raw_coverage",
    )
    deltas = {
        strategy: {key: results[i + 2][key] - results[i][key] for key in keys}
        for i, strategy in enumerate(("terminal", "v1_compatible"))
    }
    inputs = json.loads((root / "fixed-inputs.json").read_text())
    assert all(sha256_file(root / "new-hierarchy" / p) == h for p, h in inputs.items())
    benchmark = json.loads((root / "benchmark-inputs.json").read_text())
    assert all(sha256_file(root / "benchmark" / p) == h for p, h in benchmark.items())
    parity = json.loads((root / "h5-parity/report.json").read_text())
    assert parity["passed"]
    report = dict(
        evaluation_completed=True,
        strategy_parity_passed=True,
        fixed_ssn_unchanged=True,
        benchmark_unchanged=True,
        results=results,
        same_strategy_new_minus_previous=deltas,
        historical_full_v1_scope="independent full pipeline, not hierarchy-only control",
        hierarchy=json.loads((root / "h5-hierarchy-summary.json").read_text()),
        accuracy_improvement_claimed=False,
    )
    (out / "report.json").write_text(json.dumps(report, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("migrate", "summary", "score", "validate-h4"))
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--origin", type=Path)
    parser.add_argument("--component", type=int)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--parallel", type=Path)
    args = parser.parse_args()
    if args.command == "migrate":
        previous = yaml.safe_load((args.origin / "repaired-mean/run.yaml").read_text())
        migrated = migrate_config(previous)
        (args.run / "run.yaml").write_text(yaml.safe_dump(migrated, sort_keys=True))
        args.out.write_text(
            json.dumps(
                dict(previous=previous, migrated=migrated, migration_explicit=True), indent=2
            )
        )
    elif args.command == "summary":
        summary = summarize_hierarchy(
            args.run, None if args.component is None else (args.component,)
        )
        args.out.write_text(json.dumps(summary, indent=2) + "\n")
        if not summary["candidate_budget_passed"]:
            raise RuntimeError("Candidate budget exceeded")
    elif args.command == "validate-h4":
        report = validate_h4(args.run, args.parallel)
        args.out.write_text(json.dumps(report, indent=2) + "\n")
        if not report["passed"]:
            raise RuntimeError("H4 acceptance checks failed; H5 remains gated")
    else:
        score_campaign(args.run, args.origin)


if __name__ == "__main__":
    main()
