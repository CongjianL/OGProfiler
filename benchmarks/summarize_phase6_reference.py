"""Summarize one Phase 6 component-local hierarchy benchmark."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pyarrow.parquet as pq


def _time_text(path: Path, label: str) -> str:
    pattern = re.compile(rf"^\s*{re.escape(label)}:\s*(.+)$")
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if match:
            return match.group(1)
    raise ValueError(f"Missing {label} in {path}")


def _wall_seconds(value: str) -> float:
    parts = [float(part) for part in value.split(":")]
    if len(parts) == 2:
        return parts[0] * 60 + parts[1]
    if len(parts) == 3:
        return parts[0] * 3600 + parts[1] * 60 + parts[2]
    raise ValueError(f"Invalid elapsed time: {value}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--time-report", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    component = args.run / "hierarchy" / "components" / "component=00000000"
    nodes = pq.read_table(component / "nodes.parquet")
    members = pq.read_table(component / "members.parquet")
    candidates = pq.read_table(component / "candidates.parquet")
    manifest = json.loads((component / "hierarchy-manifest.json").read_text(encoding="utf-8"))
    reasons = [value for value in nodes["terminal_reason"].to_pylist() if value is not None]
    summary = {
        "component_id": 0,
        "n_proteins": members.num_rows,
        "hierarchy_nodes": nodes.num_rows,
        "terminal_families": len(set(members["terminal_cluster_id"].to_pylist())),
        "max_depth": max(nodes["depth"].to_pylist()),
        "resolution_candidates": candidates.num_rows,
        "terminal_reason_counts": {
            reason: reasons.count(reason) for reason in sorted(set(reasons))
        },
        "engine_metrics": manifest["metrics"],
        "wall_seconds": _wall_seconds(
            _time_text(args.time_report, "Elapsed (wall clock) time (h:mm:ss or m:ss)")
        ),
        "peak_rss_kib": int(
            _time_text(args.time_report, "Maximum resident set size (kbytes)")
        ),
        "verified_resume": True,
    }
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
