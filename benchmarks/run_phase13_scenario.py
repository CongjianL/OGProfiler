#!/usr/bin/env python3
"""Generate and execute one Phase 13 full-factorial synthetic scenario."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

from ogprofiler.benchmark.synthetic import (
    generate_scenario_matrix,
    generate_synthetic_dataset,
    write_scenario_matrix,
)
from ogprofiler.core.manifest import sha256_file, write_json


def _run(command: list[str], cwd: Path) -> float:
    started = time.perf_counter()
    subprocess.run(command, cwd=cwd, check=True)
    return time.perf_counter() - started


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--index", required=True, type=int)
    parser.add_argument("--replicates", type=int, default=1)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--ancestral-families", type=int, default=12)
    parser.add_argument("--sequence-length", type=int, default=120)
    args = parser.parse_args()
    project_root = Path(__file__).resolve().parents[1]
    scenarios = generate_scenario_matrix(replicates=args.replicates)
    if not 0 <= args.index < len(scenarios):
        parser.error(f"--index must be between 0 and {len(scenarios) - 1}")
    scenario = scenarios[args.index]
    if args.index == 0:
        write_scenario_matrix(args.out / "matrix", scenarios)
    output = args.out / "runs" / scenario.scenario_id
    dataset = output / "dataset"
    run_root = output / "run"
    output.mkdir(parents=True, exist_ok=True)
    manifest = generate_synthetic_dataset(
        dataset,
        scenario,
        ancestral_families=args.ancestral_families,
        sequence_length=args.sequence_length,
    )
    python = sys.executable
    base = [python, "-m", "ogprofiler.cli"]
    seconds: dict[str, float] = {}
    commands = [
        (
            "prepare",
            base + ["prepare", "--proteomes", str(dataset / "proteomes"), "--out", str(run_root)],
        ),
        (
            "search",
            base
            + [
                "search",
                "--run",
                str(run_root),
                "--backend",
                "diamond",
                "--set",
                f"search.threads={args.workers}",
                "--set",
                "search.max_target_seqs=0",
            ],
        ),
        ("edges", base + ["edges", "--run", str(run_root)]),
        ("components", base + ["components", "--run", str(run_root)]),
        (
            "hierarchy",
            base
            + [
                "hierarchy-all",
                "--run",
                str(run_root),
                "--set",
                f"runtime.workers={args.workers}",
            ],
        ),
        ("network", base + ["annotate-network", "--run", str(run_root)]),
        ("export", base + ["export", "--run", str(run_root)]),
        (
            "orthology",
            base
            + ["orthologs", "--run", str(run_root), "--emit-pairwise-orthologs"],
        ),
        (
            "metrics",
            base
            + [
                "benchmark",
                "synthetic-metrics",
                "--run",
                str(run_root),
                "--dataset-root",
                str(dataset),
                "--out",
                str(output / "synthetic-metrics.json"),
            ],
        ),
    ]
    for name, command in commands:
        seconds[name] = _run(command, project_root)
    metrics = json.loads((output / "synthetic-metrics.json").read_text(encoding="utf-8"))
    write_json(
        output / "run-summary.json",
        {
            "scenario_id": scenario.scenario_id,
            "scenario": metrics["scenario"],
            "dataset_manifest_sha256": sha256_file(manifest),
            "stage_seconds": seconds,
            "total_seconds": sum(seconds.values()),
            "headline_metrics": {
                "family_f1": metrics["family"]["pairwise_clustering"]["f1"],
                "hierarchy_f1": metrics["hierarchy"]["f1"],
                "event_accuracy": metrics["events"]["end_to_end_accuracy"],
                "orthology_f1": (metrics.get("orthology") or {}).get("f1"),
            },
        },
    )
    print(output / "run-summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
