#!/usr/bin/env python3
"""Aggregate repeated Phase 2 reference-SSN C/E benchmarks."""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
from pathlib import Path
from typing import Any

DATASETS = ("C_gene_family_expansion", "E_large_connected_component")


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def normalized_hierarchy_hash(path: Path) -> str:
    value = load_json(path)
    nodes = sorted(
        (
            {
                "cluster_hash": node["cluster_hash"],
                "parent_hashes": sorted(node["parent_hashes"]),
                "terminal": bool(node["terminal"]),
            }
            for node in value["nodes"]
        ),
        key=lambda node: (
            node["cluster_hash"],
            node["parent_hashes"],
            node["terminal"],
        ),
    )
    semantic = {"nodes": nodes, "membership": value["membership"]}
    payload = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def parse_elapsed(value: str) -> float:
    fields = [float(field) for field in value.split(":")]
    if len(fields) == 2:
        return fields[0] * 60 + fields[1]
    if len(fields) == 3:
        return fields[0] * 3600 + fields[1] * 60 + fields[2]
    raise ValueError(f"Unexpected GNU time value: {value!r}")


def parse_time(path: Path) -> dict[str, float]:
    elapsed_label = "Elapsed (wall clock) time (h:mm:ss or m:ss):"
    rss_label = "Maximum resident set size (kbytes):"
    values: dict[str, float] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line.startswith(elapsed_label):
            values["wall_seconds"] = parse_elapsed(line[len(elapsed_label) :].strip())
        elif line.startswith(rss_label):
            values["peak_rss_bytes"] = float(line[len(rss_label) :].strip()) * 1024
    if set(values) != {"wall_seconds", "peak_rss_bytes"}:
        raise ValueError(f"Incomplete GNU time report: {path}")
    return values


def median(values: list[float]) -> float:
    return float(statistics.median(values))


def collect_repeat(path: Path) -> dict[str, Any]:
    comparison = load_json(path / "comparison" / "comparison.json")
    v1_summary = load_json(path / "v1_summary" / "baseline_metrics.json")
    v1_metrics = load_json(path / "v1_hierarchy" / "metrics.json")
    v2_manifest = load_json(path / "v2_hierarchy" / "manifest.json")
    components = v2_manifest["components"]
    return {
        "repeat": path.name,
        "v1_time": parse_time(path / "v1.time"),
        "v2_time": parse_time(path / "v2.time"),
        "v1_core_seconds": float(v1_metrics["runtime_seconds"]),
        "v2_core_seconds": sum(float(item["runtime_seconds"]) for item in components),
        "v1_fingerprint": normalized_hierarchy_hash(
            path / "v1_summary" / "normalized_hierarchy.json"
        ),
        "v2_fingerprint": normalized_hierarchy_hash(
            path / "comparison" / "v2_normalized_hierarchy.json"
        ),
        "ssn": {
            "nodes": v1_summary["ssn_nodes"],
            "edges": v1_summary["ssn_edges"],
            "components": v1_summary["connected_component_count"],
            "largest_component": v1_summary["largest_component_size"],
        },
        "v2_metrics": {
            "leiden_calls": sum(int(item["leiden_calls"]) for item in components),
            "subgraph_constructions": sum(
                int(item["subgraph_constructions"]) for item in components
            ),
            "resolution_candidates": sum(
                int(item.get("resolution_candidate_count", item["leiden_calls"]))
                for item in components
            ),
            "hierarchy_nodes": sum(int(item["hierarchy_node_count"]) for item in components),
            "terminal_families": sum(int(item["terminal_family_count"]) for item in components),
        },
        "proteins": comparison["proteins"],
        "membership": comparison["terminal_pairing"],
        "topology": comparison["topology"],
    }


def aggregate_dataset(path: Path) -> dict[str, Any]:
    repeats = [collect_repeat(item) for item in sorted(path.glob("repeat_*"))]
    if len(repeats) != 3:
        raise ValueError(f"Expected three repeats under {path}, found {len(repeats)}")
    first = repeats[0]
    v1_wall = [item["v1_time"]["wall_seconds"] for item in repeats]
    v2_wall = [item["v2_time"]["wall_seconds"] for item in repeats]
    v1_core = [item["v1_core_seconds"] for item in repeats]
    v2_core = [item["v2_core_seconds"] for item in repeats]
    v1_rss = [item["v1_time"]["peak_rss_bytes"] for item in repeats]
    v2_rss = [item["v2_time"]["peak_rss_bytes"] for item in repeats]
    return {
        "dataset": path.name,
        "repeat_count": len(repeats),
        "ssn": first["ssn"],
        "determinism": {
            "v1_unique_normalized_fingerprints": len({item["v1_fingerprint"] for item in repeats}),
            "v2_unique_normalized_fingerprints": len({item["v2_fingerprint"] for item in repeats}),
            "v1_fingerprints": [item["v1_fingerprint"] for item in repeats],
            "v2_fingerprints": [item["v2_fingerprint"] for item in repeats],
        },
        "performance": {
            "v1_wall_seconds_median": median(v1_wall),
            "v1_wall_seconds_range": [min(v1_wall), max(v1_wall)],
            "v2_wall_seconds_median": median(v2_wall),
            "v2_wall_seconds_range": [min(v2_wall), max(v2_wall)],
            "wall_ratio_v1_over_v2": median(v1_wall) / median(v2_wall),
            "v1_core_seconds_median": median(v1_core),
            "v2_core_seconds_median": median(v2_core),
            "core_ratio_v1_over_v2": median(v1_core) / median(v2_core),
            "v1_peak_rss_bytes_median": median(v1_rss),
            "v2_peak_rss_bytes_median": median(v2_rss),
            "rss_ratio_v1_over_v2": median(v1_rss) / median(v2_rss),
        },
        "v2_metrics": first["v2_metrics"],
        "proteins": first["proteins"],
        "membership": first["membership"],
        "topology": first["topology"],
        "repeats": repeats,
    }


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
        "# Phase 2 final reference-SSN regression and scale benchmark",
        "",
        "## Provenance",
        "",
        f"- Run ID: `{metadata['run_id']}`",
        f"- Slurm job: `{metadata['slurm_job_id']}`",
        f"- Source SHA-256: `{metadata['source_sha256']}`",
        "- Three repeats per dataset; frozen V1 and V2 consume the identical reference SSN.",
        "- Resources per repeat: 8 CPUs, 16 GiB, one hour.",
        "",
        "## Scale and performance",
        "",
        "| Dataset | Nodes | Edges | Components | Largest | V1 wall s | V2 wall s | "
        "Wall V1/V2 | V1 core s | V2 core s | Core V1/V2 | V1 RSS MiB | V2 RSS MiB |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for item in datasets:
        ssn = item["ssn"]
        perf = item["performance"]
        lines.append(
            f"| {item['dataset']} | {ssn['nodes']} | {ssn['edges']} | {ssn['components']} | "
            f"{ssn['largest_component']} | {fmt(perf['v1_wall_seconds_median'])} | "
            f"{fmt(perf['v2_wall_seconds_median'])} | "
            f"{fmt(perf['wall_ratio_v1_over_v2'])} | "
            f"{fmt(perf['v1_core_seconds_median'], 4)} | "
            f"{fmt(perf['v2_core_seconds_median'], 4)} | "
            f"{fmt(perf['core_ratio_v1_over_v2'])} | "
            f"{fmt(perf['v1_peak_rss_bytes_median'] / 2**20, 1)} | "
            f"{fmt(perf['v2_peak_rss_bytes_median'] / 2**20, 1)} |"
        )
    lines.extend(
        [
            "",
            "Ratios above 1 favor V2. Wall values include command startup and artifact writing; "
            "core values are implementation-instrumented hierarchy computation.",
            "",
            "## Regression and topology",
            "",
            "| Dataset | V1 deterministic | V2 deterministic | Common IDs | Pair Jaccard | "
            "Node Jaccard | Edge Jaccard | Leiden calls | Candidates | Subgraphs | V2 nodes | "
            "V2 terminals |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for item in datasets:
        v1_deterministic = item["determinism"]["v1_unique_normalized_fingerprints"] == 1
        v2_deterministic = item["determinism"]["v2_unique_normalized_fingerprints"] == 1
        metrics = item["v2_metrics"]
        lines.append(
            f"| {item['dataset']} | {str(v1_deterministic).lower()} | "
            f"{str(v2_deterministic).lower()} | "
            f"{item['proteins']['common']} | {fmt(item['membership']['jaccard'])} | "
            f"{fmt(item['topology']['node_jaccard'])} | "
            f"{fmt(item['topology']['edge_jaccard'])} | {metrics['leiden_calls']} | "
            f"{metrics['resolution_candidates']} | {metrics['subgraph_constructions']} | "
            f"{metrics['hierarchy_nodes']} | "
            f"{metrics['terminal_families']} |"
        )
    lines.extend(["", "## Phase 2 exit gate", ""])
    for item in datasets:
        checks = {
            "three completed repeats": item["repeat_count"] == 3,
            "V2 normalized topology deterministic": (
                item["determinism"]["v2_unique_normalized_fingerprints"] == 1
            ),
            "V2 median RSS below V1": (
                item["performance"]["v2_peak_rss_bytes_median"]
                < item["performance"]["v1_peak_rss_bytes_median"]
            ),
            "largest reference component exercised": item["ssn"]["largest_component"]
            == (60 if item["dataset"].startswith("C_") else 240),
        }
        lines.extend([f"### {item['dataset']}", ""])
        lines.extend(f"- [{'x' if passed else ' '}] {name}" for name, passed in checks.items())
        lines.append("")
    lines.append(
        "Successful V2 command completion also means component-level hierarchy invariants passed "
        "before each artifact was written. Membership and topology differences remain explicit "
        "compatibility measurements rather than an exact-parity claim."
    )
    lines.extend(
        [
            "",
            "C produced two frozen-V1 semantic topology fingerprints across three repeats, while "
            "V2 produced one. E was deterministic in both implementations. Candidate counts in "
            "this run equal Leiden calls because the prototype evaluates exactly one candidate per "
            "call; the post-benchmark implementation additionally persists `candidates.parquet`.",
        ]
    )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--slurm-job-id", required=True)
    parser.add_argument("--source-sha256", required=True)
    args = parser.parse_args()
    run = {
        "metadata": {
            "run_id": args.run_id,
            "slurm_job_id": args.slurm_job_id,
            "source_sha256": args.source_sha256,
        },
        "datasets": [aggregate_dataset(args.run_dir / dataset) for dataset in DATASETS],
    }
    (args.run_dir / "summary.json").write_text(
        json.dumps(run, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (args.run_dir / "REPORT.md").write_text(render_report(run), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
