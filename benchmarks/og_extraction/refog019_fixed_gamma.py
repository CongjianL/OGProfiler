"""One fixed-gamma V2 robust evaluation on frozen mean SSN; no policy changes."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from dataclasses import asdict, replace
from itertools import combinations
from pathlib import Path
from unittest.mock import patch

from benchmarks.og_extraction.binary_target_audit import event
from benchmarks.og_extraction.refog_split_audit import rows
from ogprofiler.core.manifest import sha256_file
from ogprofiler.graph.partition import load_component_edge_table
from ogprofiler.hierarchy import resolution
from ogprofiler.hierarchy.engine import HierarchyConfig
from ogprofiler.hierarchy.leiden import LeidenCallCounter
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig

GAMMA = 0.953125
COMPONENT = 3976


def evaluate_fixed(graph, config, gamma=GAMMA):
    """Use the public production search seam with equal bounds; observe its calls."""
    traces = []
    delegate = resolution.run_leiden

    def observe(g, value, method, weights, seed, counter):
        result = delegate(g, value, method, weights, seed, counter)
        traces.append(
            dict(
                seed=seed,
                gamma=value,
                membership=list(result.membership),
                quality=result.quality,
                child_sizes=sorted(Counter(result.membership).values()),
            )
        )
        return result

    counter = LeidenCallCounter(n_iterations=config.leiden_iterations)
    fixed = replace(config.resolution, gamma_min=gamma, gamma_max=gamma)
    with patch.object(resolution, "run_leiden", observe):
        result = resolution.search_resolution(
            graph,
            fixed,
            method=config.method,
            weights="weight",
            seed=config.seed,
            counter=counter,
            stability_mode=config.stability_mode,
        )
    assert len(result.candidates) == 1 and counter.count == 3 and len(traces) == 3
    candidate = result.candidates[0]
    pairs = [
        dict(
            left_seed=a["seed"],
            right_seed=b["seed"],
            ari=resolution.adjusted_rand_index(tuple(a["membership"]), tuple(b["membership"])),
        )
        for a, b in combinations(traces, 2)
    ]
    assert abs(sum(p["ari"] for p in pairs) / 3 - candidate.stability) < 1e-12
    return candidate, traces, pairs, counter.count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--origin", type=Path, required=True)
    parser.add_argument("--v1-run", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    run = args.origin / "new-hierarchy"
    assert (args.origin / "provenance/job_id.txt").read_text().strip() == "1410868"
    assert (args.v1_run / "provenance/job_id.txt").read_text().strip() == "1411083"
    edge = json.loads((run / "edges/edge-manifest.json").read_text())
    assert edge["parameters"]["symmetrization"] == "mean"
    folder = run / f"hierarchy/components/component={COMPONENT:08d}"
    manifest = json.loads((folder / "hierarchy-manifest.json").read_text())
    checks = manifest["input_checksums"]
    for name, digest in checks.items():
        assert sha256_file(run / name) == digest
    for name, digest in manifest["output_checksums"].items():
        assert sha256_file(folder / name) == digest
    v1file = args.v1_run / f"v1-original/component-{COMPONENT}.json"
    v1 = json.loads(v1file.read_text())
    assert v1["input_checksums"] == checks
    v1report = json.loads((args.v1_run / "v1-original/report.json").read_text())
    parameters = dict(manifest["parameters"]["hierarchy"])
    parameters["resolution"] = ResolutionSearchConfig(**parameters["resolution"])
    config = HierarchyConfig(**parameters)
    assert config.stability_mode == "robust" and config.resolution.stability_threshold == 0.9
    assert config.seed == 42 and config.leiden_iterations == 10
    assert config.recursion_stop_size == 1 and config.resolution.max_child_fraction == 0.95
    frozen = rows(folder / "candidates.parquet", "cluster_id", [0])
    assert not any(c["gamma"] == GAMMA for c in frozen)
    table = load_component_edge_table(run / "components", COMPONENT)
    graph, ids = table.to_igraph()
    root = next(n for n in rows(folder / "nodes.parquet") if n["parent_id"] is None)
    assert root["cluster_id"] == 0 and len(ids) == root["n_genes"] == 12
    import igraph
    import leidenalg

    assert (
        manifest["parameters"]["environment"]["igraph"]
        == igraph.__version__
        == v1report["environment"]["igraph"]
    )
    assert (
        manifest["parameters"]["environment"]["leidenalg"]
        == leidenalg.__version__
        == v1report["environment"]["leidenalg"]
    )
    candidate, traces, pairs, calls = evaluate_fixed(graph, config)
    meta = {
        p["protein_id"]: p for p in rows(run / "input/proteins.parquet") if p["protein_id"] in ids
    }
    tracefile = args.v1_run / f"v1-original/component-{COMPONENT}-calls.jsonl"
    calls_v1 = [json.loads(line) for line in tracefile.read_text().splitlines()]
    last = next(c for c in reversed(calls_v1) if c["node"] == "0")
    assert last["gamma"] == GAMMA and last["child_count"] == 2
    v1map = {
        int(g.split("|g")[1]): label
        for g, label in zip(last["vertex_names"], last["membership"], strict=True)
    }
    v1membership = tuple(v1map[p] for p in ids)
    species = [meta[p]["species_id"] for p in ids]
    for trace in traces:
        membership = tuple(trace["membership"])
        trace["ari_to_v1_actual"] = resolution.adjusted_rand_index(membership, v1membership)
        trace["groups"] = [
            [
                meta[p]["original_id"]
                for p, label in zip(ids, membership, strict=True)
                if label == child
            ]
            for child in sorted(set(membership))
        ]
    classification = (
        "PASS_BINARY_AT_PREVIOUSLY_UNTESTED_POINT"
        if candidate.binary_eligible
        else "STABILITY_CONFLICT_AT_FIXED_POINT"
        if "UNSTABLE" in candidate.original_violations
        else "OTHER_GATE_OR_NON_BINARY_CONFLICT_AT_FIXED_POINT"
    )
    report = dict(
        completed=True,
        production_parameters_unchanged=True,
        component_id=COMPONENT,
        cluster_id=0,
        gamma=GAMMA,
        frozen_gamma_previously_tested=False,
        classification=classification,
        parameters=asdict(config),
        diagnostic_override={"gamma_min": GAMMA, "gamma_max": GAMMA},
        candidate=asdict(candidate),
        v1_event=event(candidate, species, False),
        traces=traces,
        pairwise_ari=pairs,
        leiden_calls=calls,
        v1_seed=None,
        v1_membership=list(v1membership),
        representative_ari_to_v1=resolution.adjusted_rand_index(candidate.membership, v1membership),
        input_checksums=checks,
        environment=manifest["parameters"]["environment"],
        v1_report_sha256=sha256_file(v1file),
        v1_trace_sha256=sha256_file(tracefile),
        source_sha256=sha256_file(Path(__file__)),
        interpretation=(
            "One fixed point under unchanged V2 gates; "
            "no proof about unsampled intervals or full descendants"
        ),
    )
    for name, digest in checks.items():
        assert sha256_file(run / name) == digest
    for name, digest in manifest["output_checksums"].items():
        assert sha256_file(folder / name) == digest
    (args.out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
