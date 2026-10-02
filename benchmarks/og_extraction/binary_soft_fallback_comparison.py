"""Controlled three-policy by two-lower-bound full-component experiment."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import pyarrow.parquet as pq

from benchmarks.og_extraction.binary_coverage_schedule import (
    FAIR,
    SOFT,
    SOFT_V2,
    UPPER_GUARD,
    UPPER_GUARD_V2,
)
from benchmarks.og_extraction.binary_recursive_audit import compare_component
from ogprofiler.core.manifest import sha256_file


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expanded", action="store_true")
    args = parser.parse_args()
    metadata = (
        pq.ParquetFile(args.result_root / "new-hierarchy/input/proteins.parquet").read().to_pylist()
    )
    sources = [
        Path(__file__),
        Path(__file__).with_name("binary_coverage_schedule.py"),
        Path(__file__).with_name("binary_recursive_audit.py"),
        Path("src/ogprofiler/hierarchy/engine.py"),
        Path("src/ogprofiler/hierarchy/leiden.py"),
        Path("src/ogprofiler/hierarchy/resolution.py"),
        Path("src/ogprofiler/orthogroups/engine.py"),
    ]
    report = {
        "protocol": "soft-upper-guard-3way-by-lower-bound-v1",
        "source_hashes": {str(p): sha256_file(p) for p in sources},
        "results": {},
        "source_snapshots": {str(p): p.read_text() for p in sources},
        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "git_tracked_diff": subprocess.check_output(["git", "diff"], text=True),
    }
    report["protocol"] += "-expanded" if args.expanded else ""
    lowers = (0.01,) if args.expanded else (0.01, 0.001)
    policies = (FAIR, UPPER_GUARD_V2, SOFT_V2) if args.expanded else (FAIR, UPPER_GUARD, SOFT)
    for lower in lowers:
        for protocol in policies:
            key = f"{lower}/{protocol}"
            report["results"][key] = {
                str(cid): compare_component(
                    args.result_root, cid, metadata, protocol, gamma_min=lower
                )
                for cid in (2694, 470, 84, 844, 500)
            }
    args.output.write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
