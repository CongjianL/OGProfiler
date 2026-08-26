"""Summarize a Phase 7 scheduler and verified-resume benchmark."""

from __future__ import annotations

import argparse
import json
import re
import sqlite3
from pathlib import Path

import pyarrow.parquet as pq


def _time_text(path: Path, label: str) -> str:
    pattern = re.compile(rf"^\s*{re.escape(label)}:\s*(.+)$")
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if match:
            return match.group(1)
    raise ValueError(f"Missing {label} in {path}")


def _wall_seconds(path: Path) -> float:
    value = _time_text(path, "Elapsed (wall clock) time (h:mm:ss or m:ss)")
    parts = [float(part) for part in value.split(":")]
    if len(parts) == 2:
        return parts[0] * 60 + parts[1]
    return parts[0] * 3600 + parts[1] * 60 + parts[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--initial-time", required=True, type=Path)
    parser.add_argument("--resume-time", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    components = sorted((args.run / "hierarchy" / "components").glob("component=*"))
    member_count = sum(pq.read_table(path / "members.parquet").num_rows for path in components)
    node_count = sum(pq.read_table(path / "nodes.parquet").num_rows for path in components)
    with sqlite3.connect(args.run / "run.db") as database:
        rows = database.execute(
            "SELECT task_id,status,attempts FROM tasks WHERE stage='hierarchy' "
            "ORDER BY CAST(task_id AS INTEGER)"
        ).fetchall()
        start_order = [
            int(row[0])
            for row in database.execute(
                "SELECT task_id FROM task_events WHERE stage='hierarchy' AND status='RUNNING' "
                "ORDER BY event_id"
            ).fetchall()
        ]
    initial = json.loads((args.run.parent / "scheduler-initial.json").read_text(encoding="utf-8"))
    resume = json.loads((args.run.parent / "scheduler-resume.json").read_text(encoding="utf-8"))
    summary = {
        "initial_counts": initial["counts"],
        "resume_counts": resume["counts"],
        "component_task_count": len(rows),
        "component_start_order": start_order,
        "task_statuses": [list(row) for row in rows],
        "hierarchy_nodes": node_count,
        "terminal_members": member_count,
        "initial_wall_seconds": _wall_seconds(args.initial_time),
        "resume_wall_seconds": _wall_seconds(args.resume_time),
        "initial_peak_rss_kib": int(
            _time_text(args.initial_time, "Maximum resident set size (kbytes)")
        ),
        "resume_peak_rss_kib": int(
            _time_text(args.resume_time, "Maximum resident set size (kbytes)")
        ),
    }
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
