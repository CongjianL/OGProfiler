"""Finite fixed-node target comparison; not a production topology policy."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, replace
from pathlib import Path

import pyarrow.parquet as pq

from ogprofiler.core.manifest import sha256_file
from ogprofiler.graph.partition import load_component_edge_table
from ogprofiler.hierarchy.engine import HierarchyConfig
from ogprofiler.hierarchy.leiden import LeidenCallCounter
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig, search_resolution
from ogprofiler.orthogroups.legacy_events import v1_event_for_node

CASES = ((2694, 0), (2694, 3), (470, 0), (84, 5), (844, 2), (500, 0))


def candidate_bank(original, config):
    """Common finite bank: coarse, endpoint, rescue and frozen local points."""
    points = list(original)
    gamma = config.gamma_min
    for _ in range(config.max_coarse_candidates):
        if gamma > config.gamma_max:
            break
        points.append(gamma)
        gamma *= config.growth_factor
    points.append(config.gamma_max)
    points += [
        config.gamma_min
        * (config.gamma_max / config.gamma_min) ** (i / (config.rescue_grid_points + 1))
        for i in range(1, config.rescue_grid_points + 1)
    ]
    unique = []
    for value in sorted(points):
        if not unique or not math.isclose(value, unique[-1], rel_tol=1e-12, abs_tol=0):
            unique.append(value)
    if len(unique) > config.max_candidate_evaluations:
        raise ValueError("Common bank exceeds candidate budget; do not silently truncate")
    return unique


def event(candidate, species, has_parent):
    bitmaps = {}
    for label, sid in zip(candidate.membership, species, strict=True):
        bitmaps[label] = bitmaps.get(label, 0) | (1 << sid)
    parent = 0
    for bitmap in bitmaps.values():
        parent |= bitmap
    if candidate.child_count < 2:
        return None
    return v1_event_for_node(
        parent, list(bitmaps.values()), len(species), candidate.child_count + int(has_parent)
    )


def compare_node(result_root, component, cluster):
    run = result_root / "new-hierarchy"
    folder = run / "hierarchy/components" / f"component={component:08d}"
    manifest = json.loads((folder / "hierarchy-manifest.json").read_text())
    parameters = dict(manifest["parameters"]["hierarchy"])
    parameters["resolution"] = ResolutionSearchConfig(**parameters["resolution"])
    config = HierarchyConfig(**parameters)
    nodes = {
        n["cluster_id"]: n for n in pq.ParquetFile(folder / "nodes.parquet").read().to_pylist()
    }
    node = nodes[cluster]

    def inside(leaf):
        while leaf is not None:
            if leaf == cluster:
                return True
            leaf = nodes[leaf]["parent_id"]
        return False

    proteins = sorted(
        m["protein_id"]
        for m in pq.ParquetFile(folder / "members.parquet").read().to_pylist()
        if inside(m["terminal_cluster_id"])
    )
    assert len(proteins) == node["n_genes"] < 1000
    for relative, digest in manifest["input_checksums"].items():
        if relative.startswith(f"components/edges/component={component:08d}/"):
            assert sha256_file(run / relative) == digest
    table = load_component_edge_table(run / "components", component)
    graph, ids = table.to_igraph()
    indices = {p: i for i, p in enumerate(ids)}
    graph = graph.induced_subgraph([indices[p] for p in proteins])
    species_by_id = {
        p["protein_id"]: p["species_id"]
        for p in pq.ParquetFile(run / "input/proteins.parquet")
        .read(columns=["protein_id", "species_id"])
        .to_pylist()
    }
    species = [species_by_id[p] for p in proteins]
    original = [
        c
        for c in pq.ParquetFile(folder / "candidates.parquet").read().to_pylist()
        if c["cluster_id"] == cluster
    ]
    chosen = next(c for c in original if c["selected"])
    points = candidate_bank([c["gamma"] for c in original], config.resolution)

    def evaluate(gamma, counter):
        result = search_resolution(
            graph,
            replace(config.resolution, gamma_min=gamma, gamma_max=gamma),
            method=config.method,
            weights="weight",
            seed=config.seed,
            counter=counter,
            stability_mode=config.stability_mode,
        )
        assert len(result.candidates) == 1
        c = result.candidates[0]
        d = asdict(c)
        d.pop("membership")
        d.update(
            v1_event=event(c, species, node["parent_id"] is not None),
            binary_target_valid=c.valid and c.child_count == 2,
        )
        return c, d

    counter = LeidenCallCounter(n_iterations=config.leiden_iterations)
    bank = [evaluate(g, counter)[1] for g in points]
    baseline = next(c for c in bank if math.isclose(c["gamma"], chosen["gamma"], rel_tol=1e-12))
    assert baseline["child_count"] == chosen["child_count"]
    assert abs(baseline["stability"] - chosen["stability"]) < 1e-12
    valid_k = [c for c in bank if c["valid"]]
    valid_b = [c for c in bank if c["binary_target_valid"]]
    # Reproducible V1-style target probe for <1000 nodes, deliberately capped at24.
    # Uses V2 robust representative/gates; not literal seed-free V1 replay.
    lower, upper, probe = 0.0, 1.0, []
    probe_counter = LeidenCallCounter(n_iterations=config.leiden_iterations)
    for _ in range(config.resolution.max_candidate_evaluations):
        gamma = (lower + upper) / 2
        if any(math.isclose(gamma, c["gamma"], rel_tol=1e-12, abs_tol=0) for c in probe):
            break
        c, row = evaluate(gamma, probe_counter)
        probe.append(row)
        if c.child_count == 2:
            break
        if c.child_count < 2:
            lower = gamma
        else:
            upper = gamma
    assert counter.count == 3 * len(bank) <= 72 and probe_counter.count == 3 * len(probe) <= 72
    return dict(
        component_id=component,
        cluster_id=cluster,
        n_genes=len(proteins),
        n_species=len(set(species)),
        edge_count=graph.ecount(),
        parameters=asdict(config),
        frozen_selected={**chosen, "v1_event": baseline["v1_event"]},
        common_bank=bank,
        bank_leiden_calls=counter.count,
        lowest_bank_kway=min(valid_k, key=lambda c: c["gamma"], default=None),
        lowest_bank_binary=min(valid_b, key=lambda c: c["gamma"], default=None),
        binary_candidates=sum(c["child_count"] == 2 for c in bank),
        valid_binary_candidates=len(valid_b),
        v1_style_probe=probe,
        probe_leiden_calls=probe_counter.count,
        probe_found_binary=probe[-1]["child_count"] == 2,
        probe_binary_passes_original_gates=probe[-1]["binary_target_valid"],
        scope="Two separate finite diagnostics, not an extra production fallback or full V1 replay",
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    records = [compare_node(args.results, *case) for case in CASES]
    args.out.write_text(
        json.dumps(
            dict(
                production_unchanged=True,
                cases=records,
                audit_source_sha256=sha256_file(Path(__file__)),
            ),
            indent=2,
        )
        + "\n"
    )
    for r in records:
        print(
            json.dumps(
                {
                    k: v
                    for k, v in r.items()
                    if k not in ("parameters", "common_bank", "v1_style_probe")
                },
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
