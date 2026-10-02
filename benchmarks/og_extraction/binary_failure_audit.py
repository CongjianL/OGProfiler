"""Bounded diagnostic sweep of frozen failed subgraphs; not the 24-point search."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from dataclasses import asdict, replace
from pathlib import Path

from ogprofiler.core.manifest import sha256_file
from ogprofiler.graph.partition import load_component_edge_table
from ogprofiler.hierarchy.engine import HierarchyConfig
from ogprofiler.hierarchy.leiden import LeidenCallCounter
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig, search_resolution


def diagnostic_points(config, original):
    """129 global log points, 65 points per first raw count transition bracket."""
    points = [
        (config.gamma_min * (config.gamma_max / config.gamma_min) ** (i / 128), "IN_RANGE_LOG")
        for i in range(129)
    ]
    ordered = sorted(original, key=lambda c: c["gamma"])
    transitions = [
        (a, b)
        for a, b in zip(ordered, ordered[1:], strict=False)
        if a["child_count"] < 2 <= b["child_count"]
    ]
    if transitions:
        a, b = transitions[0]
        points.extend(
            (a["gamma"] + (b["gamma"] - a["gamma"]) * i / 64, "IN_RANGE_TRANSITION")
            for i in range(65)
        )
    if ordered[0]["child_count"] > 2:
        # Explicit out-of-policy diagnosis, never admitted into revised production bounds.
        points.extend(
            (config.gamma_min * 10 ** (-3 + 3 * i / 32), "BELOW_MIN_DIAGNOSTIC_ONLY")
            for i in range(33)
        )
    unique = []
    for gamma, scope in points:
        if not any(math.isclose(gamma, g, rel_tol=1e-12, abs_tol=0) for g, _ in unique):
            unique.append((gamma, scope))
    return unique


def audit(root, comparison):
    run = root / "new-hierarchy"
    report = {
        "diagnostic_protocol": "failure-coverage-log129-transition65-v1",
        "source_sha256": sha256_file(Path(__file__)),
        "comparison_sha256": sha256_file(comparison),
        "nodes": [],
    }
    frozen = json.loads(comparison.read_text())
    for cid, data in frozen["components"].items():
        cfg = dict(data["config"])
        cfg["resolution"] = ResolutionSearchConfig(**cfg["resolution"])
        cfg = HierarchyConfig(**cfg)
        table = load_component_edge_table(run / "components", int(cid))
        manifest = json.loads(
            (
                run
                / "hierarchy/components"
                / f"component={int(cid):08d}"
                / "hierarchy-manifest.json"
            ).read_text()
        )
        for relative, digest in manifest["input_checksums"].items():
            if relative.startswith(f"components/edges/component={int(cid):08d}/"):
                assert sha256_file(run / relative) == digest
        root_graph, ids = table.to_igraph()
        indices = {p: i for i, p in enumerate(ids)}
        for node in data["hard_binary"]["unresolved"]:
            proteins = sorted(
                p
                for p, leaf in data["hard_binary"]["terminal_membership"]
                if leaf == node["cluster_id"]
            )
            assert len(proteins) == node["n_genes"]
            graph = root_graph.induced_subgraph([indices[p] for p in proteins])
            original = [
                c
                for c in data["hard_binary"]["candidates"]
                if c["cluster_id"] == node["cluster_id"]
            ]
            counter = LeidenCallCounter(n_iterations=cfg.leiden_iterations)
            rows = []
            for gamma, scope in diagnostic_points(cfg.resolution, original):
                result = search_resolution(
                    graph,
                    replace(cfg.resolution, gamma_min=gamma, gamma_max=gamma),
                    method=cfg.method,
                    weights="weight",
                    seed=cfg.seed,
                    counter=counter,
                    stability_mode=cfg.stability_mode,
                )
                assert len(result.candidates) == 1
                rows.append({"scope": scope, **asdict(result.candidates[0])})
            binary = [r for r in rows if r["child_count"] == 2]
            accepted = [r for r in binary if r["valid"] and r["gamma"] >= cfg.resolution.gamma_min]
            report["nodes"].append(
                {
                    "component": int(cid),
                    "cluster_id": node["cluster_id"],
                    "protein_ids": proteins,
                    "edges": graph.ecount(),
                    "connected_components": len(graph.connected_components()),
                    "diagnostic_points": len(rows),
                    "leiden_calls": counter.count,
                    "count_distribution": dict(Counter(r["child_count"] for r in rows)),
                    "in_range_valid_binary": accepted,
                    "binary_candidates": binary,
                    "candidates": rows,
                }
            )
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-root", type=Path, required=True)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(audit(args.result_root, args.comparison), indent=2))


if __name__ == "__main__":
    main()
