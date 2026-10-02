"""ADR0004 experiment only: shared-budget hard-binary search and full DFS audit."""

from __future__ import annotations

import argparse
import json
import math
import time
from collections import Counter
from dataclasses import asdict, replace
from functools import partial
from pathlib import Path
from unittest.mock import patch

import pyarrow.parquet as pq

from benchmarks.og_extraction.binary_coverage_schedule import (
    FAIR,
    SOFT,
    SOFT_V2,
    STABILITY,
    UPPER_GUARD,
    UPPER_GUARD_V2,
    coverage_search,
)
from benchmarks.og_extraction.binary_schedule_v2 import (
    PROTOCOL as PROTOCOL_V2,
)
from benchmarks.og_extraction.binary_schedule_v2 import (
    finite_binary_search_v2,
)
from ogprofiler.core.manifest import sha256_file
from ogprofiler.graph.components import Component
from ogprofiler.graph.partition import load_component_edge_table
from ogprofiler.hierarchy.engine import HierarchyConfig, infer_component_hierarchy
from ogprofiler.hierarchy.resolution import (
    ResolutionSearchConfig,
    ResolutionSearchResult,
    search_resolution,
)
from ogprofiler.hierarchy.validation import validate_hierarchy
from ogprofiler.orthogroups.engine import extract_component_orthogroups

PROTOCOL = "adr0004-hard-binary-24-v1"


def finite_binary_search(config, evaluate):
    """Public experimental seam: evaluator returns one original-gated candidate."""
    config.validate()
    tested = []
    blocked = False

    def sample(gamma, phase):
        nonlocal blocked
        for c in tested:
            if math.isclose(c.gamma, gamma, rel_tol=1e-12, abs_tol=0):
                return c
        if len(tested) >= min(24, config.max_candidate_evaluations):
            blocked = True
            return None
        c = evaluate(gamma)
        violations = c.violations + (() if c.child_count == 2 else ("TARGET_CHILD_COUNT",))
        c = replace(
            c,
            valid=c.valid and c.child_count == 2,
            policy_valid=c.policy_valid and c.child_count == 2,
            violations=violations,
            rejection_reason=c.rejection_reason or (violations[0] if violations else None),
            phase=phase,
            evaluation_index=len(tested) + 1,
            evaluation_budget=min(24, config.max_candidate_evaluations),
        )
        tested.append(c)
        return c

    def accepted():
        return sorted((c for c in tested if c.valid), key=lambda c: c.gamma)

    gamma = config.gamma_min
    for _ in range(10):
        if gamma > config.gamma_max:
            break
        c = sample(gamma, "COARSE")
        if c is None or c.valid:
            break
        gamma *= config.growth_factor
    if not accepted() and not blocked:
        sample(config.gamma_max, "ENDPOINT")
    # Adjacent raw-count brackets only; invalid binary keeps both sides eligible.
    for _ in range(5):
        if accepted() or blocked:
            break
        ordered = sorted(tested, key=lambda c: c.gamma)
        brackets = [
            (a, b)
            for a, b in zip(ordered, ordered[1:], strict=False)
            if (a.child_count < 2 <= b.child_count) or a.child_count == 2 or b.child_count == 2
        ]
        if not brackets:
            break
        a, b = min(brackets, key=lambda ab: (ab[0].gamma, ab[1].gamma))
        midpoint = (a.gamma + b.gamma) / 2
        if any(math.isclose(midpoint, c.gamma, rel_tol=1e-12, abs_tol=0) for c in tested):
            break
        sample(midpoint, "TARGET_PROBE")
    if not accepted() and not blocked:
        for i in range(1, 6):
            gamma = config.gamma_min * (config.gamma_max / config.gamma_min) ** (i / 6)
            if sample(gamma, "RESCUE") is None:
                break
    if accepted():
        upper = accepted()[0].gamma
        lower = max((c.gamma for c in tested if c.gamma < upper), default=config.gamma_min)
        # Endpoints are already evaluated; exactly three possible new points.
        for i in range(1, 4):
            if sample(lower + (upper - lower) * i / 4, "REFINE") is None:
                break
    chosen = accepted()
    status = (
        "ACCEPTED"
        if chosen
        else "EVALUATION_BUDGET_EXHAUSTED"
        if blocked
        else "REJECTED_ALL_TESTED"
    )
    return ResolutionSearchResult(
        chosen[0] if chosen else None, tuple(tested), None if chosen else status, status
    )


def uses_binary_target(size):
    return 3 <= size < 10000


def experimental_search(graph, config, *, protocol=PROTOCOL, **kwargs):
    if protocol not in (
        PROTOCOL,
        PROTOCOL_V2,
        FAIR,
        STABILITY,
        UPPER_GUARD,
        SOFT,
        UPPER_GUARD_V2,
        SOFT_V2,
    ):
        raise ValueError("Unknown experimental protocol")
    if not uses_binary_target(graph.vcount()):
        return search_resolution(graph, config, **kwargs)

    def evaluate(gamma):
        result = search_resolution(
            graph, replace(config, gamma_min=gamma, gamma_max=gamma), **kwargs
        )
        assert len(result.candidates) == 1
        return result.candidates[0]

    if protocol in (FAIR, STABILITY, UPPER_GUARD, SOFT, UPPER_GUARD_V2, SOFT_V2):
        return coverage_search(config, evaluate, protocol=protocol)
    return (finite_binary_search_v2 if protocol == PROTOCOL_V2 else finite_binary_search)(
        config, evaluate
    )


def compare_component(root, cid, metadata, protocol=PROTOCOL, *, gamma_min=None):
    search = partial(experimental_search, protocol=protocol)
    run = root / "new-hierarchy"
    folder = run / "hierarchy/components" / f"component={cid:08d}"
    manifest = json.loads((folder / "hierarchy-manifest.json").read_text())
    for relative, digest in manifest["input_checksums"].items():
        if relative.startswith(f"components/edges/component={cid:08d}/"):
            assert sha256_file(run / relative) == digest
    parameters = dict(manifest["parameters"]["hierarchy"])
    parameters["resolution"] = ResolutionSearchConfig(**parameters["resolution"])
    config = HierarchyConfig(**parameters)
    experiment_config = (
        config
        if gamma_min is None
        else replace(config, resolution=replace(config.resolution, gamma_min=gamma_min))
    )
    experiment_config.validate()
    table = load_component_edge_table(run / "components", cid)
    graph, ids = table.to_igraph()
    component = Component(cid, tuple(ids), table.edges)
    species = {r["protein_id"]: r["species_id"] for r in metadata}
    originals = {r["protein_id"]: r["original_id"] for r in metadata}
    reports = {}
    partitions = {}
    for name in ("baseline", "hard_binary"):
        run_config = config if name == "baseline" else experiment_config
        started = time.perf_counter()
        # Scoped serial-only benchmark adapter; production sources/defaults unchanged.
        with patch(
            "ogprofiler.hierarchy.engine.search_resolution",
            search_resolution if name == "baseline" else search,
        ):
            result = infer_component_hierarchy(component, run_config, species, graph)
        validate_hierarchy(component, result)
        if name == "baseline":
            baseline_membership = result.terminal_membership
        nodes = [asdict(n) for n in result.nodes]
        unresolved = [n for n in nodes if n["split_status"] == "UNRESOLVED"]
        counts = Counter(c.cluster_id for c in result.resolution_candidates)
        report = {
            "seconds": time.perf_counter() - started,
            "metrics": asdict(result.metrics),
            "nodes": nodes,
            "candidates": [asdict(c) for c in result.resolution_candidates],
            "terminal_membership": result.terminal_membership,
            "unresolved": unresolved,
            "coverage": len(result.terminal_membership),
            "max_candidates_per_node": max(counts.values(), default=0),
            "og_status": "BLOCKED_UNRESOLVED" if unresolved else "EVALUATED",
            "refinement_truncated_nodes": [
                n.cluster_id
                for n in result.nodes
                if n.split_status == "SPLIT"
                and uses_binary_target(n.n_genes)
                and n.child_count == 2
                and counts[n.cluster_id] == 24
                and sum(
                    c.phase == "REFINE"
                    for c in result.resolution_candidates
                    if c.cluster_id == n.cluster_id
                )
                < 3
            ],
        }
        report["fallback_nodes"] = [
            n
            for n in nodes
            if n["selection_phase"] and n["selection_phase"].startswith("FALLBACK_KWAY/")
        ]
        report["candidate_eligibility"] = [
            {
                "cluster_id": c.cluster_id,
                "gamma": c.gamma,
                "original_violations": [v for v in c.violations if v != "TARGET_CHILD_COUNT"],
                "binary_eligible": c.structural_valid
                and c.child_count == 2
                and not [v for v in c.violations if v != "TARGET_CHILD_COUNT"],
                "kway_eligible": c.structural_valid
                and c.child_count >= 2
                and not [v for v in c.violations if v != "TARGET_CHILD_COUNT"],
                "selected": c.selected,
                "selection_phase": c.phase,
            }
            for c in result.resolution_candidates
        ]
        assert len(result.terminal_membership) == len(ids)
        assert {p for p, _ in result.terminal_membership} == set(ids)
        assert report["max_candidates_per_node"] <= 24
        assert result.metrics.leiden_calls == 3 * len(result.resolution_candidates)
        if name == "hard_binary":
            with patch("ogprofiler.hierarchy.engine.search_resolution", search):
                repeated = infer_component_hierarchy(component, run_config, species, graph)
            validate_hierarchy(component, repeated)
            report["repeat_identical"] = (
                repeated.nodes == result.nodes
                and repeated.terminal_membership == result.terminal_membership
                and repeated.resolution_candidates == result.resolution_candidates
            )
            assert report["repeat_identical"]
        if not unresolved:
            members = [
                {"protein_id": p, "terminal_cluster_id": c} for p, c in result.terminal_membership
            ]
            og = extract_component_orthogroups(
                nodes, members, species, originals, total_species=len(set(species.values()))
            )
            report["og"] = asdict(og)
            partitions[name] = sorted(sorted(g.protein_ids) for g in og.groups)
        reports[name] = report
    frozen_nodes = pq.ParquetFile(folder / "nodes.parquet").read().to_pylist()
    reports["baseline_matches_frozen_nodes"] = all(
        len(reports["baseline"]["nodes"]) == len(frozen_nodes)
        and all(
            math.isclose(actual[key], value, rel_tol=1e-12, abs_tol=1e-12)
            if key == "quality" and value is not None
            else json.loads(json.dumps(actual[key])) == value
            for key, value in frozen.items()
        )
        for actual, frozen in zip(reports["baseline"]["nodes"], frozen_nodes, strict=True)
    )
    frozen_membership = pq.ParquetFile(folder / "members.parquet").read().to_pylist()
    reports["baseline_matches_frozen_membership"] = sorted(
        result_pair for result_pair in baseline_membership
    ) == sorted((r["protein_id"], r["terminal_cluster_id"]) for r in frozen_membership)
    frozen_members = (
        pq.ParquetFile(
            root
            / "new-hierarchy/orthogroups/components"
            / f"component={cid:08d}"
            / "members.parquet"
        )
        .read()
        .to_pylist()
    )
    groups = {}
    for row in frozen_members:
        groups.setdefault(row["local_group_id"], []).append(row["protein_id"])
    reports["baseline_matches_frozen_og"] = partitions.get("baseline") == sorted(
        sorted(p) for p in groups.values()
    )
    reports["baseline_config"] = asdict(config)
    reports["config"] = asdict(experiment_config)
    return reports


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--protocol",
        choices=(
            PROTOCOL,
            PROTOCOL_V2,
            FAIR,
            STABILITY,
            UPPER_GUARD,
            SOFT,
            UPPER_GUARD_V2,
            SOFT_V2,
        ),
        default=PROTOCOL,
    )
    args = parser.parse_args()
    metadata = (
        pq.ParquetFile(args.result_root / "new-hierarchy/input/proteins.parquet").read().to_pylist()
    )
    report = {
        "protocol": args.protocol,
        "coverage_schedule_sha256": sha256_file(
            Path(__file__).with_name("binary_coverage_schedule.py")
        ),
        "schedule_v2_sha256": sha256_file(Path(__file__).with_name("binary_schedule_v2.py")),
        "source_sha256": sha256_file(Path(__file__)),
        "production_source_hashes": {
            str(path): sha256_file(path)
            for path in (
                Path("src/ogprofiler/hierarchy/engine.py"),
                Path("src/ogprofiler/hierarchy/resolution.py"),
                Path("src/ogprofiler/hierarchy/leiden.py"),
                Path("src/ogprofiler/orthogroups/engine.py"),
            )
        },
        "components": {
            str(c): compare_component(args.result_root, c, metadata, args.protocol)
            for c in (2694, 470, 84, 844, 500)
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
