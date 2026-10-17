"""Replay frozen small components through public serial/spawn production paths."""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from collections import Counter
from dataclasses import asdict
from pathlib import Path

import pyarrow.parquet as pq

from benchmarks.og_extraction.fixed_hierarchy import compare_component
from benchmarks.og_extraction.reference_v1 import load_reference
from ogprofiler.core.manifest import sha256_file
from ogprofiler.graph.components import Component
from ogprofiler.graph.partition import load_component_edge_table
from ogprofiler.hierarchy.engine import HierarchyConfig, infer_component_hierarchy
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig
from ogprofiler.hierarchy.subtree import infer_component_hierarchy_parallel
from ogprofiler.hierarchy.validation import validate_hierarchy
from ogprofiler.storage.hierarchy import write_hierarchy_result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    run = args.result_root / "new-hierarchy"
    frozen = json.loads(
        (args.result_root / "binary-soft-fallback-expanded-comparison.json").read_text()
    )
    frozen = frozen["results"]["0.01/adr0004-soft24-upper-guard-v2"]
    metadata = pq.ParquetFile(run / "input/proteins.parquet").read().to_pylist()
    species = {r["protein_id"]: r["species_id"] for r in metadata}
    originals = {r["protein_id"]: r["original_id"] for r in metadata}
    sources = list(Path("src/ogprofiler").rglob("*.py")) + [Path(__file__)]
    provenance = dict(
        git_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        git_status=subprocess.check_output(["git", "status", "--short"], text=True),
        source_hashes={str(p): sha256_file(p) for p in sources},
        input_hashes={
            str(p): sha256_file(p)
            for p in (
                run / "input/proteins.parquet",
                run / "components/index.parquet",
                args.result_root / "binary-soft-fallback-expanded-comparison.json",
            )
        },
        configurations={},
    )
    report = {}
    reference = load_reference()
    for cid in (2694, 470, 84, 844, 500):
        parameters = dict(frozen[str(cid)]["config"])
        parameters["resolution"] = ResolutionSearchConfig(
            **parameters["resolution"], topology_policy="soft_binary_24_v2"
        )
        config = HierarchyConfig(**parameters)
        provenance["configurations"][str(cid)] = asdict(config)
        for p in (run / "components/edges" / f"component={cid:08d}").rglob("*.parquet"):
            provenance["input_hashes"][str(p)] = sha256_file(p)
        table = load_component_edge_table(run / "components", cid)
        graph, ids = table.to_igraph()
        component = Component(cid, tuple(ids), table.edges)
        start = time.perf_counter()
        serial = infer_component_hierarchy(component, config, species, graph)
        serial_seconds = time.perf_counter() - start
        start = time.perf_counter()
        parallel = infer_component_hierarchy_parallel(component, config, species, graph, workers=2)
        parallel_seconds = time.perf_counter() - start
        for label, result in (("serial", serial), ("parallel", parallel)):
            validate_hierarchy(component, result)
            write_hierarchy_result(args.output_dir / label / f"component={cid:08d}", result)
        equal = (
            serial.nodes == parallel.nodes
            and serial.terminal_membership == parallel.terminal_membership
            and serial.resolution_candidates == parallel.resolution_candidates
            and serial.metrics.leiden_calls == parallel.metrics.leiden_calls
        )
        parquet_equal = all(
            pq.ParquetFile(args.output_dir / "serial" / f"component={cid:08d}" / name)
            .read()
            .equals(
                pq.ParquetFile(args.output_dir / "parallel" / f"component={cid:08d}" / name).read()
            )
            for name in ("nodes.parquet", "members.parquet", "candidates.parquet")
        )
        assert serial.metrics.leiden_calls == 3 * len(serial.resolution_candidates)
        assert len(members := serial.terminal_membership) == len(ids)
        assert {p for p, _ in members} == set(ids)
        baseline = frozen[str(cid)]["hard_binary"]

        def matches(rows, old):
            return len(rows) == len(old) and all(
                all(json.loads(json.dumps(row[k])) == v for k, v in prior.items())
                for row, prior in zip(rows, old, strict=True)
            )

        frozen_equal = (
            matches([asdict(n) for n in serial.nodes], baseline["nodes"])
            and matches([asdict(c) for c in serial.resolution_candidates], baseline["candidates"])
            and json.loads(json.dumps(serial.terminal_membership))
            == baseline["terminal_membership"]
        )
        members = [
            {"protein_id": p, "terminal_cluster_id": c} for p, c in serial.terminal_membership
        ]
        _, parity = compare_component(
            [asdict(n) for n in serial.nodes],
            members,
            species,
            originals,
            len(set(species.values())),
            set(),
            reference_env=reference,
        )
        counts = Counter(c.cluster_id for c in serial.resolution_candidates)
        report[str(cid)] = dict(
            serial_parallel_equal=equal,
            parquet_equal=parquet_equal,
            frozen_experiment_equal=frozen_equal,
            v1_parity=parity,
            serial_seconds=serial_seconds,
            parallel_seconds=parallel_seconds,
            nodes=len(serial.nodes),
            coverage=len(members),
            calls=serial.metrics.leiden_calls,
            max_candidates=max(counts.values(), default=0),
            fallbacks=sum(n.selection_kind == "FALLBACK_KWAY" for n in serial.nodes),
            unresolved=sum(n.split_status == "UNRESOLVED" for n in serial.nodes),
        )
    (args.output_dir / "provenance.json").write_text(json.dumps(provenance, indent=2))
    (args.output_dir / "report.json").write_text(json.dumps(report, indent=2))
    assert all(
        r["serial_parallel_equal"]
        and r["parquet_equal"]
        and r["frozen_experiment_equal"]
        and r["v1_parity"]["passed"]
        and not r["unresolved"]
        and r["max_candidates"] <= 24
        for r in report.values()
    ), report
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
