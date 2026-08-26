#!/usr/bin/env python3
"""Run one strided Phase 14 component-manifest shard with verified resume."""

from __future__ import annotations

import argparse
from pathlib import Path

import pyarrow.parquet as pq

from ogprofiler.cli import main as ogprofiler_main


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--shard", type=int, required=True)
    parser.add_argument("--shards", type=int, required=True)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--subtree-workers", type=int, default=1)
    parser.add_argument("--subtree-release-size", type=int, default=50_000)
    args = parser.parse_args()
    rows = pq.read_table(args.manifest, memory_map=True).to_pylist()
    selected = [row for row in rows if int(row["array_index"]) % args.shards == args.shard]
    for row in selected:
        command = [
            "hierarchy",
            "--run",
            str(args.run),
            "--component-id",
            str(int(row["component_id"])),
        ]
        if args.config is not None:
            command.extend(["--config", str(args.config)])
        command.extend(
            [
                "--set",
                f"hierarchy.subtree_workers={args.subtree_workers}",
                "--set",
                f"hierarchy.subtree_release_size={args.subtree_release_size}",
            ]
        )
        status = ogprofiler_main(command)
        if status != 0:
            return status
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
