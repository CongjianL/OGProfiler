"""Controlled root and full-tree audit of lower bounds, not a production default."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, replace
from pathlib import Path

import pyarrow.parquet as pq

from benchmarks.og_extraction.binary_coverage_schedule import FAIR
from benchmarks.og_extraction.binary_recursive_audit import compare_component, experimental_search
from ogprofiler.core.manifest import sha256_file
from ogprofiler.graph.partition import load_component_edge_table
from ogprofiler.hierarchy.engine import HierarchyConfig
from ogprofiler.hierarchy.leiden import LeidenCallCounter
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig

LOWER_BOUNDS = (0.01, 0.001, 0.0001)


def compare_root(root, lower):
    run = root / "new-hierarchy"
    manifest = json.loads(
        (run / "hierarchy/components/component=00000084" / "hierarchy-manifest.json").read_text()
    )
    for relative, digest in manifest["input_checksums"].items():
        if relative.startswith("components/edges/component=00000084/"):
            assert sha256_file(run / relative) == digest
    cfg = dict(manifest["parameters"]["hierarchy"])
    cfg["resolution"] = ResolutionSearchConfig(**cfg["resolution"])
    cfg = HierarchyConfig(**cfg)
    cfg = replace(cfg, resolution=replace(cfg.resolution, gamma_min=lower))
    cfg.validate()
    graph, ids = load_component_edge_table(run / "components", 84).to_igraph()
    outputs = []
    for _ in range(2):
        counter = LeidenCallCounter(n_iterations=cfg.leiden_iterations)
        result = experimental_search(
            graph,
            cfg.resolution,
            protocol=FAIR,
            method=cfg.method,
            weights="weight",
            seed=cfg.seed,
            counter=counter,
            stability_mode=cfg.stability_mode,
        )
        assert len(result.candidates) <= 24
        assert counter.count == 3 * len(result.candidates)
        assert all(lower <= c.gamma <= cfg.resolution.gamma_max for c in result.candidates)
        outputs.append(result)
    assert outputs[0] == outputs[1]
    return {
        "config": asdict(cfg),
        "protein_ids": list(ids),
        "result": asdict(outputs[0]),
        "repeat_identical": True,
        "calls": 3 * len(outputs[0].candidates),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
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
        "protocol": "component84-lower-bound-3way-v1",
        "search_protocol": FAIR,
        "source_hashes": {str(p): sha256_file(p) for p in sources},
        "root": {str(lower): compare_root(args.result_root, lower) for lower in LOWER_BOUNDS},
        "recursive": {
            str(lower): compare_component(args.result_root, 84, metadata, FAIR, gamma_min=lower)
            for lower in LOWER_BOUNDS
        },
    }
    args.output.write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
