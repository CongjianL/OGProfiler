#!/usr/bin/env python3
"""Execute one Phase 12 matrix row from shared prepared inputs and search hits."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

from ogprofiler.benchmark.matrix import generate_ofat_matrix, write_matrix
from ogprofiler.core.manifest import sha256_file, write_json


def _run(command: list[str], cwd: Path) -> float:
    started = time.perf_counter()
    subprocess.run(command, cwd=cwd, check=True)
    return time.perf_counter() - started


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--proteomes", required=True, type=Path)
    parser.add_argument("--hits", required=True, type=Path)
    parser.add_argument("--ground-truth", required=True, type=Path)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--index", required=True, type=int)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    project_root = Path(__file__).resolve().parents[1]
    matrix_root = args.out / "matrix"
    runs = generate_ofat_matrix()
    write_matrix(matrix_root, runs)
    if not 0 <= args.index < len(runs):
        parser.error(f"--index must be between 0 and {len(runs) - 1}")
    matrix_run = runs[args.index]
    output = args.out / "runs" / matrix_run.run_id
    run_root = output / "run"
    output.mkdir(parents=True, exist_ok=True)
    python = sys.executable
    base = [python, "-m", "ogprofiler.cli"]
    overrides = [item for value in matrix_run.overrides for item in ("--set", value)]
    stage_seconds: dict[str, float] = {}
    stage_seconds["prepare"] = _run(
        base
        + [
            "prepare",
            "--proteomes",
            str(args.proteomes.resolve()),
            "--out",
            str(run_root.resolve()),
        ],
        project_root,
    )
    search_root = run_root / "search"
    search_root.mkdir(exist_ok=True)
    shutil.copyfile(args.hits, search_root / "hits.parquet")
    commands = [
        ("edges", base + ["edges", "--run", str(run_root.resolve()), *overrides]),
        ("components", base + ["components", "--run", str(run_root.resolve()), *overrides]),
        (
            "hierarchy",
            base
            + [
                "hierarchy-all",
                "--run",
                str(run_root.resolve()),
                "--set",
                f"runtime.workers={args.workers}",
                *overrides,
            ],
        ),
        (
            "evolution",
            base + ["annotate-network", "--run", str(run_root.resolve()), *overrides],
        ),
        ("export", base + ["export", "--run", str(run_root.resolve())]),
    ]
    for name, command in commands:
        stage_seconds[name] = _run(command, project_root)
    metrics_path = output / "scientific-metrics.json"
    stage_seconds["metrics"] = _run(
        base
        + [
            "benchmark",
            "metrics",
            "--run",
            str(run_root.resolve()),
            "--ground-truth",
            str(args.ground_truth.resolve()),
            "--dataset",
            args.dataset,
            "--method",
            f"ogprofiler2:{matrix_run.run_id}",
            "--out",
            str(metrics_path.resolve()),
        ],
        project_root,
    )
    metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    write_json(
        output / "run-summary.json",
        {
            "run_id": matrix_run.run_id,
            "axis": matrix_run.axis,
            "value": matrix_run.value,
            "overrides": list(matrix_run.overrides),
            "config_sha256": matrix_run.config_sha256,
            "input_checksums": {
                "hits.parquet": sha256_file(args.hits),
                "ground_truth.tsv": sha256_file(args.ground_truth),
            },
            "stage_seconds": stage_seconds,
            "total_seconds": sum(stage_seconds.values()),
            "headline_metrics": {
                "family_pair_f1": metrics["family"]["pairwise_clustering"]["f1"],
                "exact_family_recovery": metrics["family"]["exact_family_recovery"]["rate"],
                "single_copy_recovery": metrics["family"]["single_copy_family_recovery"]["rate"],
                "species_overlap_consistency": metrics["evolution"][
                    "species_overlap_consistency"
                ],
            },
        },
    )
    print(output / "run-summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
