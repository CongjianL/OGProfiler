#!/usr/bin/env python3
"""Aggregate completed Phase 13 scenario metrics into an applicability map."""

from __future__ import annotations

import argparse
from pathlib import Path

from ogprofiler.benchmark.synthetic_metrics import aggregate_applicability


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    paths = sorted((args.root / "runs").glob("*/synthetic-metrics.json"))
    json_path, _ = aggregate_applicability(paths, args.out)
    print(json_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
