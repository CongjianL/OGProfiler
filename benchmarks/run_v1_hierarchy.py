#!/usr/bin/env python3
"""Run only the frozen V1 hierarchy stage on an existing SSN."""

from __future__ import annotations

import argparse
import importlib.util
import json
import resource
import sys
import time
from pathlib import Path
from types import ModuleType

import igraph


def load_frozen_v1(path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location("ogprofiler_v1_frozen", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load frozen V1 module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def peak_rss_bytes() -> int:
    maximum = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return maximum if sys.platform == "darwin" else maximum * 1024


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--v1", required=True, type=Path)
    parser.add_argument("--ssn", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--method", default="rber")
    parser.add_argument("--weight", default="NBS")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--gamma-coefficient", type=float, default=1.0)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    module = load_frozen_v1(args.v1)
    graph = igraph.Graph.Read_GML(str(args.ssn))
    started = time.perf_counter()
    hierarchy = module.Call.ConstructedHnCC(
        graph,
        args.method,
        args.weight,
        args.threads,
        args.gamma_coefficient,
        str(args.out),
    )
    runtime = time.perf_counter() - started
    metrics = {
        "algorithm": "frozen-v1-hierarchy",
        "runtime_seconds": runtime,
        "peak_rss_bytes": peak_rss_bytes(),
        "ssn_nodes": graph.vcount(),
        "ssn_edges": graph.ecount(),
        "hierarchy_nodes": hierarchy.vcount(),
        "hierarchy_edges": hierarchy.ecount(),
        "method": args.method,
        "weight": args.weight,
        "threads": args.threads,
        "gamma_coefficient": args.gamma_coefficient,
    }
    (args.out / "metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(metrics, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

