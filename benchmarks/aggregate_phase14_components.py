#!/usr/bin/env python3
"""Verify and merge Phase 14 array outputs into a component result index."""

from __future__ import annotations

import argparse
from pathlib import Path

from ogprofiler.benchmark.extreme_scale import merge_component_task_results


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    table = merge_component_task_results(args.manifest, args.run, args.out)
    print(f"verified_components={table.num_rows} result_index={args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
