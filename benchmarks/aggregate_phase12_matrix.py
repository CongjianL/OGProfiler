#!/usr/bin/env python3
"""Aggregate completed Phase 12 rows and select an evidence-backed default."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import mean
from typing import Any

from ogprofiler.core.manifest import write_json


def _score(metrics: dict[str, Any]) -> float:
    values = [
        metrics["family_pair_f1"],
        metrics["exact_family_recovery"],
        metrics["single_copy_recovery"],
        metrics["species_overlap_consistency"],
    ]
    available = [float(value) for value in values if value is not None]
    return mean(available) if available else 0.0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    rows: list[dict[str, Any]] = []
    for path in sorted((args.root / "runs").glob("*/run-summary.json")):
        value = json.loads(path.read_text(encoding="utf-8"))
        row = {
            "run_id": value["run_id"],
            "axis": value["axis"],
            "value": value["value"],
            **value["headline_metrics"],
            "total_seconds": value["total_seconds"],
        }
        row["scientific_score"] = _score(value["headline_metrics"])
        rows.append(row)
    rows.sort(key=lambda row: (-row["scientific_score"], row["run_id"]))
    baseline = next((row for row in rows if row["run_id"] == "baseline"), None)
    recommended = rows[0] if rows else None
    args.out.mkdir(parents=True, exist_ok=True)
    write_json(
        args.out / "matrix-summary.json",
        {
            "schema_version": "phase12-matrix-summary-v1",
            "completed_runs": len(rows),
            "baseline": baseline,
            "recommended": recommended,
            "rows": rows,
        },
    )
    fields = [
        "run_id",
        "axis",
        "value",
        "family_pair_f1",
        "exact_family_recovery",
        "single_copy_recovery",
        "species_overlap_consistency",
        "scientific_score",
        "total_seconds",
    ]
    with (args.out / "matrix-summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
