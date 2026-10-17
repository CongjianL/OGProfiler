"""Read-only event eligibility and observed production paths for pure-LCA losses."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import pyarrow.parquet as pq

from benchmarks.og_extraction.embleya_reference import read_of
from benchmarks.qfo.fixed_tree import boundary_losses
from ogprofiler.core.manifest import sha256_file


def classify(node, event, trace):
    eligible = event is None or event == "None" or (event == "I" and node["n_species"] > 1)
    statuses = {t["status"] for t in trace}
    if not eligible:
        if "SELECTED" in statuses:
            raise ValueError("Ineligible node has selected trace")
        return "excluded_by_event_qualification"
    if "SELECTED" in statuses:
        return "eligible_selected_but_reference_fragmented"
    if statuses & {"DESCENDANT_CONSUMED", "SKIPPED_CONSUMED"}:
        return "eligible_consumed_without_selection"
    return "eligible_without_observed_selection"


def diagnose(run, reference_dir, baseline, out):
    baseline_hash = sha256_file(baseline)
    prior = json.loads(baseline.read_text())
    boundary_path = baseline.parent / "boundary-summary.json"
    boundary_hash = sha256_file(boundary_path)
    boundary = json.loads(boundary_path.read_text())
    frozen = dict(prior["input_hashes"])
    frozen[str(baseline)] = baseline_hash
    frozen[str(boundary_path)] = boundary_hash
    frozen[str(run / "edges/retained_edges.parquet")] = boundary["fixed_graph_sha256"]
    rows = []

    def read(path):
        digest = sha256_file(path)
        if str(path) in frozen and frozen[str(path)] != digest:
            raise ValueError("Baseline input changed")
        frozen[str(path)] = digest
        return pq.ParquetFile(path).read().to_pylist()

    proteins = read(run / "input/proteins.parquet")
    keys = {(r["species_id"], r["original_id"]): r["protein_id"] for r in proteins}
    names = {r["species_name"]: r["species_id"] for r in read(run / "input/species.parquet")}
    ref = {keys[k]: v for k, v in read_of(reference_dir / "Orthogroups.tsv", names, keys).items()}
    with (run / "results/members.tsv").open() as h:
        actual = {int(r["protein_id"]): r["family_id"] for r in csv.DictReader(h, delimiter="\t")}
    components = Counter(r["component_id"] for r in read(run / "components/index.parquet"))
    for cid, size in sorted(components.items()):
        if size == 1:
            continue
        folder = run / "hierarchy/components" / f"component={cid:08d}"
        nodes = read(folder / "nodes.parquet")
        membership = [
            (r["protein_id"], r["terminal_cluster_id"]) for r in read(folder / "members.parquet")
        ]
        losses = {}
        boundary_losses(nodes, membership, ref, actual, losses)
        if not losses:
            continue
        og = run / "orthogroups/components" / folder.name
        manifest_path = og / "og-manifest.json"
        manifest = json.loads(manifest_path.read_text())
        frozen[str(manifest_path)] = sha256_file(manifest_path)
        if manifest["status"] != "DONE":
            raise ValueError("Unfinished production component")
        for relative, digest in manifest["input_checksums"].items():
            path = run / relative
            if digest == "MISSING" or sha256_file(path) != digest:
                raise ValueError("Production input checksum mismatch")
            if str(path) in frozen and frozen[str(path)] != digest:
                raise ValueError("Production and hierarchy input identity differs")
            frozen[str(path)] = digest
        for name, digest in manifest["output_checksums"].items():
            path = og / name
            if sha256_file(path) != digest:
                raise ValueError("Production artifact checksum mismatch")
            frozen[str(path)] = digest
        events = {r["cluster_id"]: r["v1_event"] for r in read(og / "v1_events.parquet")}
        traces = {}
        for t in read(og / "selection_trace.parquet"):
            traces.setdefault(t["cluster_id"], []).append(t)
        member_counts, exported_families = Counter(), {}
        seen = set()
        for member in read(og / "members.parquet"):
            p, g = member["protein_id"], member["local_group_id"]
            if p in seen or p not in actual:
                raise ValueError("Invalid production membership")
            seen.add(p)
            member_counts[g] += 1
            if g in exported_families and exported_families[g] != actual[p]:
                raise ValueError("Local OG differs from exported production family")
            exported_families[g] = actual[p]
        groups = {}
        for g in read(og / "groups.parquet"):
            if member_counts[g["local_group_id"]] != g["n_genes"]:
                raise ValueError("Production group size mismatch")
            g["exported_family_id"] = exported_families[g["local_group_id"]]
            groups.setdefault(g["source_cluster_id"], []).append(g)
        for node in nodes:
            c = node["cluster_id"]
            if c not in losses:
                continue
            path = sorted(traces.get(c, []), key=lambda t: t["trace_order"])
            rows.append(
                dict(
                    component_id=cid,
                    cluster_id=c,
                    **losses[c],
                    n_genes=node["n_genes"],
                    n_species=node["n_species"],
                    v1_event=events[c],
                    category=classify(node, events[c], path),
                    production_trace=path,
                    sourced_groups=groups.get(c, []),
                )
            )
    counts, weighted, event_counts = Counter(), Counter(), Counter()
    for r in rows:
        counts[r["category"]] += 1
        weighted[r["category"]] += r["lost_pairs"]
        event_counts[str(r["v1_event"])] += r["lost_pairs"]
    # Independent previous audit total, not a newly redefined scoring universe.
    expected = prior["actual"]["lost_pairs"]
    if (
        sum(weighted.values()) > expected
        or sum(weighted.values()) != boundary["lca_losses"]["pure_assigned_lca_tree_available"]
    ):
        raise ValueError("Pure LCA losses exceed all production loss")
    for path, digest in frozen.items():
        if sha256_file(Path(path)) != digest:
            raise ValueError("Frozen inputs changed")
    out.mkdir(parents=True, exist_ok=False)
    (out / "pure-lca-nodes.json").write_text(json.dumps(rows, indent=2) + "\n")
    summary = dict(
        diagnostic_only=True,
        production_changed=False,
        immutable_inputs_verified=True,
        source_run=str(run),
        baseline_report_sha256=sha256_file(baseline),
        pure_lca_lost_pairs=sum(weighted.values()),
        pure_lca_nodes=len(rows),
        category_nodes=dict(counts),
        category_lost_pairs=dict(weighted),
        event_lost_pairs=dict(event_counts),
        interpretation=(
            "Observed qualification/trace paths, not causal attribution; "
            "selected nodes may emit active-view subsets"
        ),
        top_nodes=sorted(rows, key=lambda r: -r["lost_pairs"])[:30],
        input_hashes_file="input-hashes.json",
    )
    (out / "input-hashes.json").write_text(json.dumps(frozen, indent=2) + "\n")
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("run", "reference-dir", "baseline", "out"):
        p.add_argument("--" + key, type=Path, required=True)
    a = p.parse_args()
    diagnose(a.run, a.reference_dir, a.baseline, a.out)


if __name__ == "__main__":
    main()
