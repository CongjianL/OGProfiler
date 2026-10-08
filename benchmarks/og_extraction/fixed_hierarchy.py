"""P5: independent V1 execution on a frozen, component-local hierarchy.

The reference alone materializes descendant lists (required by frozen V1).
No production event labels or extraction outputs construct the oracle.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import resource
import time
import tracemalloc
from collections import Counter, defaultdict
from pathlib import Path

import pyarrow.parquet as pq

from benchmarks.og_extraction.reference_v1 import REFERENCE_SHA256, load_reference, run_reference
from ogprofiler.core.manifest import sha256_file
from ogprofiler.exceptions import HierarchyError
from ogprofiler.orthogroups.engine import OrthogroupConflictError, extract_component_orthogroups
from ogprofiler.orthogroups.legacy_events import annotate_v1_events
from ogprofiler.orthogroups.models import OrthogroupConfig
from ogprofiler.orthogroups.stage import run_orthogroup_stage
from ogprofiler.output.stage import run_export_stage


def rows(path):
    return pq.ParquetFile(path).read().to_pylist()


def reference_spec(nodes, members, species, total_species, isolates):
    """Independent iterative postorder adapter; preserve original row order."""
    children = defaultdict(list)
    genes = defaultdict(list)
    for row in members:
        genes[row["terminal_cluster_id"]].append(row["protein_id"])
    for row in nodes:
        if row["parent_id"] is not None:
            children[row["parent_id"]].append(row["cluster_id"])
    root = next(row["cluster_id"] for row in nodes if row["parent_id"] is None)
    stack = [(root, False)]
    while stack:
        cluster, exiting = stack.pop()
        if exiting:
            if children[cluster]:
                genes[cluster] = [p for c in children[cluster] for p in genes[c]]
            continue
        stack.append((cluster, True))
        stack.extend((c, False) for c in reversed(children[cluster]))
    return dict(
        n_species=total_species,
        isolates=[str(p) for p in sorted(isolates)],
        vertices=[
            dict(
                name=str(n["cluster_id"]),
                genes=[str(p) for p in genes[n["cluster_id"]]],
                species=[str(s) for s in sorted({species[p] for p in genes[n["cluster_id"]]})],
            )
            for n in nodes
        ],
        edges=[
            (str(n["parent_id"]), str(n["cluster_id"])) for n in nodes if n["parent_id"] is not None
        ],
    )


def measured(call):
    tracemalloc.start()
    start = time.perf_counter()
    try:
        value = call()
        return value, dict(
            seconds=time.perf_counter() - start,
            python_peak_bytes=tracemalloc.get_traced_memory()[1],
        )
    finally:
        tracemalloc.stop()


def compare_component(
    nodes, members, species, originals, total_species, isolates, overlap=0, reference_env=None
):
    reference, ref_perf = measured(
        lambda: run_reference(
            reference_spec(nodes, members, species, total_species, isolates),
            overlap,
            reference_env=reference_env,
        )
    )

    def engine():
        try:
            return extract_component_orthogroups(
                nodes,
                members,
                species,
                originals,
                total_species=total_species,
                overlap_count=overlap,
                ssn_isolates=isolates,
            ), ()
        except OrthogroupConflictError as error:
            return error.result, error.duplicate_members

    (actual, duplicates), new_perf = measured(engine)
    annotations = annotate_v1_events(nodes, members, species, overlap)

    def canonical(proteins):
        return tuple(sorted((species[int(p)], originals[int(p)]) for p in proteins))

    ref_groups = [(g["level"], canonical(g["members"])) for g in reference["groups"]]
    new_groups = [(g.processing_level, canonical(g.protein_ids)) for g in actual.groups]
    checks = dict(
        events={str(a.cluster_id): a.v1_event for a in annotations} == reference["raw_events"],
        ordered_members=ref_groups == new_groups,
        member_multiset=Counter(ref_groups) == Counter(new_groups),
        unassigned=canonical(p.protein_id for p in actual.unassigned)
        == canonical(reference["unassigned"]),
        duplicates=canonical(duplicates) == canonical(reference["duplicate_members"]),
        remaining=sorted(map(str, actual.remaining_cluster_ids))
        == sorted(reference["remaining_nodes"]),
    )
    # Isolate source keys belong to the SSN namespace, not the hierarchy namespace.
    new_sources = [str(g.source_cluster_id) for g in actual.groups if g.processing_level != 0]
    ref_sources = [g["source"] for g in reference["groups"] if g["level"] != 0]
    checks["selection_sources"] = new_sources == ref_sources
    consumed = defaultdict(list)
    for trace in actual.trace:
        if trace.status == "DESCENDANT_CONSUMED":
            consumed[str(trace.consumed_by)].append(str(trace.cluster_id))
    checks["consumption"] = all(
        sorted(consumed[g["source"]] + [g["source"]]) == sorted(map(str, g["consumed"]))
        for g in reference["groups"]
        if g["level"] > 1
    )
    first = next(
        (
            dict(order=i, reference=ref, actual=new)
            for i, (ref, new) in enumerate(zip(ref_groups, new_groups, strict=False))
            if ref != new
        ),
        None,
    )
    if first is None and len(ref_groups) != len(new_groups):
        i = min(len(ref_groups), len(new_groups))
        first = dict(
            order=i,
            reference=ref_groups[i] if i < len(ref_groups) else None,
            actual=new_groups[i] if i < len(new_groups) else None,
        )
    first_event = next(
        (
            dict(
                cluster_id=a.cluster_id,
                reference=reference["raw_events"][str(a.cluster_id)],
                actual=a.v1_event,
            )
            for a in annotations
            if a.v1_event != reference["raw_events"][str(a.cluster_id)]
        ),
        None,
    )
    return actual, dict(
        passed=all(checks.values()) and not duplicates,
        checks=checks,
        duplicate_count=len(duplicates),
        classification=(
            "REFERENCE_OVERLAP"
            if duplicates
            else "MATCH"
            if all(checks.values())
            else "STRATEGY_DIVERGENCE"
        ),
        first_member_divergence=first,
        first_event_divergence=first_event,
        reference=ref_perf,
        engine=new_perf,
    )


def audit(run, out, overlap=0):
    out.mkdir(parents=True, exist_ok=False)
    metadata = rows(run / "input/proteins.parquet")
    species = {r["protein_id"]: r["species_id"] for r in metadata}
    originals = {r["protein_id"]: r["original_id"] for r in metadata}
    total_species = len(rows(run / "input/species.parquet"))
    by_component = defaultdict(list)
    for row in rows(run / "components/index.parquet"):
        by_component[row["component_id"]].append(row["protein_id"])
    isolates = {
        r["protein_id"] for r in rows(run / "components/singleton_terminal_families.parquet")
    }
    frozen = {
        str(p.relative_to(run)): sha256_file(p)
        for base in ("input", "edges", "components", "hierarchy")
        for p in sorted((run / base).rglob("*"))
        if p.is_file()
    }
    if (run / "run.yaml").is_file():
        frozen["run.yaml"] = sha256_file(run / "run.yaml")
    (out / "frozen-inputs.json").write_text(json.dumps(frozen, indent=2))
    failed, expected, unassigned = [], Counter(), set()
    sizes, covers = Counter(), Counter()
    peak_component, peak_python = 0, 0
    reference_env = load_reference()
    with (out / "components.jsonl").open("w") as stream:
        for component, proteins in sorted(by_component.items()):
            folder = run / "hierarchy/components" / f"component={component:08d}"
            if (folder / "nodes.parquet").is_file():
                nodes, members = rows(folder / "nodes.parquet"), rows(folder / "members.parquet")
            else:
                assert len(proteins) == 1 and proteins[0] in isolates
                nodes = [
                    dict(
                        component_id=component,
                        cluster_id=0,
                        parent_id=None,
                        depth=0,
                        n_genes=1,
                        n_species=1,
                    )
                ]
                members = [dict(protein_id=proteins[0], terminal_cluster_id=0)]
            assert {m["protein_id"] for m in members} == set(proteins)
            try:
                result, report = compare_component(
                    nodes,
                    members,
                    species,
                    originals,
                    total_species,
                    isolates & set(proteins),
                    overlap,
                    reference_env,
                )
            except (IndexError, ValueError, RuntimeError, HierarchyError) as error:
                report = dict(
                    passed=False,
                    classification="INPUT_OR_REFERENCE_EXCEPTION",
                    component_id=component,
                    error_type=type(error).__name__,
                    error=str(error),
                )
                failed.append(report)
                stream.write(json.dumps(report) + "\n")
                stream.flush()
                continue
            report.update(component_id=component, proteins=len(proteins), nodes=len(nodes))
            stream.write(json.dumps(report) + "\n")
            stream.flush()
            if not report["passed"]:
                failed.append(report)
            for g in result.groups:
                expected[g.membership_hash] += 1
                sizes[g.n_genes] += 1
                covers[g.n_species] += 1
            unassigned.update(p.protein_id for p in result.unassigned)
            peak_component = max(peak_component, len(proteins))
            peak_python = max(peak_python, report["engine"]["python_peak_bytes"])
    artifact_passed = False
    artifact_error = None
    if not failed:
        config = OrthogroupConfig(species_overlap_count=overlap)
        try:
            run_orthogroup_stage(run, config, ["p5-fixed-hierarchy"])
            run_export_stage(run, ["p5-fixed-hierarchy"], config=config)
        except Exception as error:
            artifact_error = dict(error_type=type(error).__name__, error=str(error))
    if not failed and artifact_error is None:
        persisted = Counter()
        for folder in sorted((run / "orthogroups/components").glob("component=*")):
            persisted.update(r["membership_hash"] for r in rows(folder / "groups.parquet"))
        exported = defaultdict(list)
        with (run / "results/members.tsv").open() as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                exported[row["family_id"]].append((int(row["species_id"]), row["original_id"]))
        hashes = Counter(
            hashlib.sha256(
                json.dumps(sorted(m), ensure_ascii=False, separators=(",", ":")).encode()
            ).hexdigest()
            for m in exported.values()
        )
        with (run / "results/unassigned.tsv").open() as stream:
            exported_unassigned = {
                int(r["protein_id"]) for r in csv.DictReader(stream, delimiter="\t")
            }
        artifact_passed = expected == persisted == hashes and unassigned == exported_unassigned
    unchanged = all(sha256_file(run / p) == digest for p, digest in frozen.items())
    summary = dict(
        passed=not failed and artifact_passed and unchanged,
        reference_sha256=REFERENCE_SHA256,
        components=len(by_component),
        first_divergence=failed[0] if failed else None,
        failed_components=len(failed),
        artifacts_passed=artifact_passed,
        artifact_error=artifact_error,
        frozen_inputs_unchanged=unchanged,
        group_size_distribution=dict(sorted(sizes.items())),
        species_coverage_distribution=dict(sorted(covers.items())),
        large_og_threshold=1000,
        large_ogs=sum(n for s, n in sizes.items() if s >= 1000),
        max_og_size=max(sizes, default=0),
        max_component_proteins=peak_component,
        max_engine_python_peak_bytes=peak_python,
        process_peak_rss_native_units=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        memory_scope=(
            "component-local engine tracemalloc; process RSS includes oracle "
            "and export; no isolated native engine RSS claim"
        ),
        ssn_reference_view=(
            "opaque protein tokens; degree-zero isolates exact, other "
            "vertices chained; extraction only queries SSN degree-zero"
        ),
        scope="fixed hierarchy V1 strategy regression, not Orthobench accuracy",
    )
    (out / "report.json").write_text(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.run, args.out)
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if report["passed"] else 1)
