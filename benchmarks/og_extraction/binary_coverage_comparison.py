"""Fixed-node then whole-component comparison of predeclared coverage schedules."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import pyarrow.parquet as pq

from benchmarks.og_extraction.binary_coverage_schedule import FAIR, STABILITY
from benchmarks.og_extraction.binary_recursive_audit import compare_component, experimental_search
from benchmarks.og_extraction.binary_schedule_v2 import PROTOCOL as V2
from ogprofiler.core.manifest import sha256_file
from ogprofiler.graph.partition import load_component_edge_table
from ogprofiler.hierarchy.engine import HierarchyConfig
from ogprofiler.hierarchy.leiden import LeidenCallCounter
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig


def fixed_node(root, comparison, cid, cluster):
    frozen = comparison["components"][str(cid)]
    node_by_id = {n["cluster_id"]: n for n in frozen["hard_binary"]["nodes"]}

    def inside(leaf):
        while leaf is not None:
            if leaf == cluster:
                return True
            leaf = node_by_id[leaf]["parent_id"]
        return False

    proteins = sorted(p for p, leaf in frozen["hard_binary"]["terminal_membership"] if inside(leaf))
    assert len(proteins) == node_by_id[cluster]["n_genes"]
    run = root / "new-hierarchy"
    table = load_component_edge_table(run / "components", cid)
    manifest = json.loads(
        (
            run / "hierarchy/components" / f"component={cid:08d}" / "hierarchy-manifest.json"
        ).read_text()
    )
    for relative, digest in manifest["input_checksums"].items():
        if relative.startswith(f"components/edges/component={cid:08d}/"):
            assert sha256_file(run / relative) == digest
    graph, ids = table.to_igraph()
    index = {p: i for i, p in enumerate(ids)}
    graph = graph.induced_subgraph([index[p] for p in proteins])
    cfg = dict(frozen["config"])
    cfg["resolution"] = ResolutionSearchConfig(**cfg["resolution"])
    cfg = HierarchyConfig(**cfg)
    results = {}
    for protocol in (V2, FAIR, STABILITY):
        outputs = []
        for _ in range(2):
            counter = LeidenCallCounter(n_iterations=cfg.leiden_iterations)
            result = experimental_search(
                graph,
                cfg.resolution,
                protocol=protocol,
                method=cfg.method,
                weights="weight",
                seed=cfg.seed,
                counter=counter,
                stability_mode=cfg.stability_mode,
            )
            assert len(result.candidates) <= 24
            assert counter.count == len(result.candidates) * 3
            outputs.append(result)
        assert outputs[0] == outputs[1]
        results[protocol] = {
            "result": asdict(outputs[0]),
            "repeat_identical": True,
            "refinement_truncated": outputs[0].selected is not None
            and len(outputs[0].candidates) == 24
            and sum(c.phase == "REFINE" for c in outputs[0].candidates) < 3,
            "leiden_calls": len(outputs[0].candidates) * 3,
        }
    return {
        "component": cid,
        "cluster_id": cluster,
        "protein_ids": proteins,
        "config": asdict(cfg),
        "protocols": results,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    comparison_path = args.result_root / "binary-recursive-comparison.json"
    comparison = json.loads(comparison_path.read_text())
    metadata = (
        pq.ParquetFile(args.result_root / "new-hierarchy/input/proteins.parquet").read().to_pylist()
    )
    sources = [
        Path(__file__),
        Path(__file__).with_name("binary_coverage_schedule.py"),
        Path(__file__).with_name("binary_recursive_audit.py"),
        Path(__file__).with_name("binary_schedule_v2.py"),
    ]
    sources += [
        Path("src/ogprofiler/hierarchy/engine.py"),
        Path("src/ogprofiler/hierarchy/resolution.py"),
        Path("src/ogprofiler/hierarchy/leiden.py"),
        Path("src/ogprofiler/orthogroups/engine.py"),
    ]
    report = {
        "source_hashes": {str(p): sha256_file(p) for p in sources},
        "comparison_sha256": sha256_file(comparison_path),
        "fixed_nodes": [
            fixed_node(args.result_root, comparison, cid, cluster)
            for cid, cluster in ((2694, 7), (2694, 0), (470, 25))
        ],
        "recursive": {
            protocol: {
                str(cid): compare_component(args.result_root, cid, metadata, protocol)
                for cid in (2694, 470, 84, 844, 500)
            }
            for protocol in (FAIR, STABILITY)
        },
    }
    args.output.write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
