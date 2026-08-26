#!/usr/bin/env python3
"""Build the immutable largest-first component manifest for a Phase 14 run."""

from __future__ import annotations

import argparse
from pathlib import Path

from ogprofiler.benchmark.extreme_scale import build_component_task_manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    tasks = build_component_task_manifest(
        args.run / "components" / "statistics.parquet", args.out
    )
    print(f"tasks={len(tasks)} manifest={args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
