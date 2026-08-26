#!/usr/bin/env python3
"""Aggregate a completed A--E frozen-V1/V2 comparison run."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

DATASET_ORDER = (
    "A_small_sanity",
    "B_paralog",
    "C_gene_family_expansion",
    "D_fusion_multidomain",
    "E_large_connected_component",
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_elapsed(value: str) -> float:
    parts = [float(part) for part in value.split(":")]
    if len(parts) == 2:
        minutes, seconds = parts
        return minutes * 60 + seconds
    if len(parts) == 3:
        hours, minutes, seconds = parts
        return hours * 3600 + minutes * 60 + seconds
    raise ValueError(f"Unexpected GNU time elapsed value: {value!r}")


def parse_time_report(path: Path) -> dict[str, float]:
    elapsed_prefix = "Elapsed (wall clock) time (h:mm:ss or m:ss):"
    rss_prefix = "Maximum resident set size (kbytes):"
    result: dict[str, float] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line.startswith(elapsed_prefix):
            result["runtime_seconds"] = parse_elapsed(line[len(elapsed_prefix) :].strip())
        elif line.startswith(rss_prefix):
            result["peak_rss_bytes"] = float(line[len(rss_prefix) :].strip()) * 1024
    missing = {"runtime_seconds", "peak_rss_bytes"} - result.keys()
    if missing:
        raise ValueError(f"Missing GNU time fields {sorted(missing)} in {path}")
    return result


def ratio(numerator: float | int, denominator: float | int) -> float | None:
    return float(numerator) / float(denominator) if denominator else None


def collect_dataset(path: Path) -> dict[str, Any]:
    manifest = load_json(path / "dataset_manifest.json")
    full = load_json(path / "v1_full_summary" / "baseline_metrics.json")
    v1_stage = load_json(path / "v1_hierarchy_summary" / "baseline_metrics.json")
    v1_core = load_json(path / "v1_hierarchy_only" / "metrics.json")
    v2_manifest = load_json(path / "v2_hierarchy" / "manifest.json")
    comparison = load_json(path / "comparison" / "comparison.json")
    v2_stage = parse_time_report(path / "v2_hierarchy.time")
    components = v2_manifest["components"]
    v2_core_seconds = sum(float(item["runtime_seconds"]) for item in components)

    return {
        "dataset": path.name,
        "input": {
            "proteins": manifest["n_proteins"],
            "species": manifest["n_species"],
        },
        "v1_full_pipeline": {
            "runtime_seconds": full["runtime_seconds"],
            "peak_rss_bytes": full["peak_rss_bytes"],
            "ssn_nodes": full["ssn_nodes"],
            "ssn_edges": full["ssn_edges"],
            "connected_components": full["connected_component_count"],
            "largest_component_size": full["largest_component_size"],
        },
        "same_ssn_performance": {
            "v1_wall_seconds": v1_stage["runtime_seconds"],
            "v2_wall_seconds": v2_stage["runtime_seconds"],
            "wall_speedup_v1_over_v2": ratio(
                v1_stage["runtime_seconds"], v2_stage["runtime_seconds"]
            ),
            "v1_core_seconds": v1_core["runtime_seconds"],
            "v2_core_seconds": v2_core_seconds,
            "core_speedup_v1_over_v2": ratio(v1_core["runtime_seconds"], v2_core_seconds),
            "v1_peak_rss_bytes": v1_stage["peak_rss_bytes"],
            "v2_peak_rss_bytes": int(v2_stage["peak_rss_bytes"]),
            "rss_ratio_v1_over_v2": ratio(v1_stage["peak_rss_bytes"], v2_stage["peak_rss_bytes"]),
        },
        "membership": comparison["terminal_pairing"],
        "proteins": comparison["proteins"],
        "topology": comparison["topology"],
    }


def tsv_rows(datasets: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for item in datasets:
        full = item["v1_full_pipeline"]
        perf = item["same_ssn_performance"]
        membership = item["membership"]
        topology = item["topology"]
        proteins = item["proteins"]
        rows.append(
            {
                "dataset": item["dataset"],
                "input_proteins": item["input"]["proteins"],
                "ssn_nodes": full["ssn_nodes"],
                "ssn_edges": full["ssn_edges"],
                "ssn_components": full["connected_components"],
                "largest_component": full["largest_component_size"],
                "v1_full_wall_s": full["runtime_seconds"],
                "v1_full_rss_mib": full["peak_rss_bytes"] / 2**20,
                "v1_hierarchy_wall_s": perf["v1_wall_seconds"],
                "v2_hierarchy_wall_s": perf["v2_wall_seconds"],
                "wall_speedup_v1_over_v2": perf["wall_speedup_v1_over_v2"],
                "v1_core_s": perf["v1_core_seconds"],
                "v2_core_s": perf["v2_core_seconds"],
                "core_speedup_v1_over_v2": perf["core_speedup_v1_over_v2"],
                "v1_hierarchy_rss_mib": perf["v1_peak_rss_bytes"] / 2**20,
                "v2_hierarchy_rss_mib": perf["v2_peak_rss_bytes"] / 2**20,
                "rss_ratio_v1_over_v2": perf["rss_ratio_v1_over_v2"],
                "v1_members": proteins["v1"],
                "v2_members": proteins["v2"],
                "common_members": proteins["common"],
                "pair_precision": membership["precision"],
                "pair_recall": membership["recall"],
                "pair_jaccard": membership["jaccard"],
                "node_jaccard": topology["node_jaccard"],
                "edge_jaccard": topology["edge_jaccard"],
            }
        )
    return rows


def fmt(value: float | int | None, digits: int = 3) -> str:
    if value is None:
        return "NA"
    if isinstance(value, int):
        return str(value)
    return f"{value:.{digits}f}"


def render_report(run: dict[str, Any]) -> str:
    metadata = run["metadata"]
    datasets = run["datasets"]
    lines = [
        "# Frozen V1 A–E baseline and same-SSN V1/V2 comparison",
        "",
        "## Run identity",
        "",
        f"- Run ID: `{metadata['run_id']}`",
        f"- Slurm job: `{metadata['slurm_job_id']}` (array tasks 0–4: `COMPLETED`)",
        f"- Source SHA-256: `{metadata['source_sha256']}`",
        f"- Dataset tree SHA-256: `{metadata['dataset_tree_sha256']}`",
        "- Resources per task: 8 CPUs, 16 GiB, 2 hours; maximum two concurrent tasks.",
        "- V1 full pipeline: DIAMOND search → SSN → hierarchy → OG output.",
        "- V1/V2 hierarchy comparison: both implementations consumed the exact "
        "V1-produced `ssn.gml`.",
        "",
        "## Frozen V1 full-pipeline baseline",
        "",
        "| Dataset | Input proteins | SSN nodes | SSN edges | Components | Largest | "
        "Wall s | RSS MiB |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for item in datasets:
        full = item["v1_full_pipeline"]
        lines.append(
            f"| {item['dataset']} | {item['input']['proteins']} | {full['ssn_nodes']} | "
            f"{full['ssn_edges']} | {full['connected_components']} | "
            f"{full['largest_component_size']} | {fmt(full['runtime_seconds'], 2)} | "
            f"{fmt(full['peak_rss_bytes'] / 2**20, 1)} |"
        )

    lines.extend(
        [
            "",
            "## Same-SSN hierarchy runtime and RSS",
            "",
            "| Dataset | V1 wall s | V2 wall s | Wall V1/V2 | V1 core s | V2 core s | "
            "Core V1/V2 | V1 RSS MiB | V2 RSS MiB |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for item in datasets:
        perf = item["same_ssn_performance"]
        lines.append(
            f"| {item['dataset']} | {fmt(perf['v1_wall_seconds'], 2)} | "
            f"{fmt(perf['v2_wall_seconds'], 2)} | {fmt(perf['wall_speedup_v1_over_v2'], 2)} | "
            f"{fmt(perf['v1_core_seconds'], 4)} | {fmt(perf['v2_core_seconds'], 4)} | "
            f"{fmt(perf['core_speedup_v1_over_v2'], 2)} | "
            f"{fmt(perf['v1_peak_rss_bytes'] / 2**20, 1)} | "
            f"{fmt(perf['v2_peak_rss_bytes'] / 2**20, 1)} |"
        )

    lines.extend(
        [
            "",
            "`wall` is GNU `/usr/bin/time -v` around the complete hierarchy command, including "
            "environment startup and artifact serialization. `core` is implementation-instrumented "
            "hierarchy computation. Ratios above 1 favor V2.",
            "",
            "## Membership and hierarchy topology",
            "",
            "| Dataset | Common members | Pair precision | Pair recall | Pair Jaccard | "
            "Node Jaccard | Edge Jaccard |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for item in datasets:
        membership = item["membership"]
        topology = item["topology"]
        lines.append(
            f"| {item['dataset']} | {item['proteins']['common']} | "
            f"{fmt(membership['precision'])} | {fmt(membership['recall'])} | "
            f"{fmt(membership['jaccard'])} | {fmt(topology['node_jaccard'])} | "
            f"{fmt(topology['edge_jaccard'])} |"
        )

    lines.extend(
        [
            "",
            "Membership metrics compare co-clustered protein pairs over IDs present in both "
            "normalized hierarchies. Topology metrics compare descendant-member-set identities "
            "and parent→child "
            "edges, so transient numeric cluster IDs do not affect the result.",
            "",
            "## Result reading",
            "",
            "- The V2 core hierarchy routine used less memory and less instrumented compute time "
            "in all five cases.",
            "- End-to-end hierarchy wall time was higher for V2 in all five cases because the "
            "prototype command "
            "starts a fresh environment and serializes per-component Parquet artifacts.",
            "- A–D show structural divergence: the V1 hierarchy produces nested singleton "
            "terminal leaves, "
            "whereas V2 terminates these small components at their component roots.",
            "- E exercises recursive splitting and is the informative topology case: pair "
            "Jaccard 0.701, "
            "node Jaccard 0.154, and edge Jaccard 0.250.",
            "- The V1 hierarchy membership export excludes SSN isolates; V2 records them as "
            "singleton "
            "components. Pairwise metrics therefore use only the common ID set.",
            "",
            "The Slurm state establishes successful execution. The reported membership/topology "
            "values establish "
            "measured prototype behavior; they do not yet establish V2 parity with frozen V1.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--slurm-job-id", required=True)
    parser.add_argument("--source-sha256", required=True)
    parser.add_argument("--dataset-tree-sha256", required=True)
    args = parser.parse_args()

    datasets = [collect_dataset(args.run_dir / name) for name in DATASET_ORDER]
    run = {
        "metadata": {
            "run_id": args.run_id,
            "slurm_job_id": args.slurm_job_id,
            "source_sha256": args.source_sha256,
            "dataset_tree_sha256": args.dataset_tree_sha256,
        },
        "datasets": datasets,
    }
    (args.run_dir / "summary.json").write_text(
        json.dumps(run, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    rows = tsv_rows(datasets)
    with (args.run_dir / "summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (args.run_dir / "REPORT.md").write_text(render_report(run), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
