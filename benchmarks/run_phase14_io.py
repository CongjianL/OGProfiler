#!/usr/bin/env python3
"""Run the Phase 14 Parquet layout tracer on a deterministic edge-shaped table."""

from __future__ import annotations

import argparse
import random
from pathlib import Path

import pyarrow as pa

from ogprofiler.benchmark.extreme_scale import (
    ParquetLayout,
    profile_parquet_layouts,
    write_parquet_profiles,
)


def _table(rows: int, seed: int) -> pa.Table:
    rng = random.Random(seed)
    left = list(range(rows))
    return pa.table(
        {
            "u": pa.array(left, type=pa.int64()),
            "v": pa.array((value + rng.randrange(1, 1000) for value in left), type=pa.int64()),
            "weight": pa.array((rng.random() for _ in left), type=pa.float64()),
        }
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--rows", type=int, default=2_000_000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    profiles = profile_parquet_layouts(
        _table(args.rows, args.seed),
        args.out / "layouts",
        tuple(
            ParquetLayout(compression, row_group_size)
            for compression in ("none", "snappy", "zstd")
            for row_group_size in (16_384, 65_536, 262_144)
        ),
        seed=args.seed,
    )
    write_parquet_profiles(args.out / "profiles.parquet", profiles)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
