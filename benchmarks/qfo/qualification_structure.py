"""Fixed-tree label-free structural features; reference cohorts are diagnostic only."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from benchmarks.og_extraction.embleya_reference import choose2, read_of
from benchmarks.og_extraction.tree_cut import topology
from benchmarks.qfo.fixed_tree import boundary_losses
from ogprofiler.core.manifest import sha256_file


def structural_features(nodes, membership, species, edges):
    by_id, children, order = topology(nodes)
    leaf = dict(membership)
    if len(leaf) != len(membership) or any(c not in by_id or children[c] for c in leaf.values()):
        raise ValueError("Invalid terminal membership")
    sizes, strength, cross, within = Counter(), Counter(), Counter(), Counter()
    sets = {c: set() for c in by_id}
    depth = {}
    for c in order:
        parent = by_id[c]["parent_id"]
        depth[c] = 0 if parent is None else depth[parent] + 1
    for p, c in membership:
        sizes[c] += 1
        sets[c].add(species[p])
    seen = set()
    for u, v, w in edges:
        if (
            u not in leaf
            or v not in leaf
            or u >= v
            or (u, v) in seen
            or not math.isfinite(w)
            or w < 0
        ):
            raise ValueError("Invalid canonical edge")
        seen.add((u, v))
        a, b = leaf[u], leaf[v]
        strength[a] += w
        strength[b] += w
        while depth[a] > depth[b]:
            a = by_id[a]["parent_id"]
        while depth[b] > depth[a]:
            b = by_id[b]["parent_id"]
        while a != b:
            a, b = by_id[a]["parent_id"], by_id[b]["parent_id"]
        cross[a] += w
    result = {}
    for c in reversed(order):
        for ch in children[c]:
            sizes[c] += sizes[ch]
            strength[c] += strength[ch]
            sets[c].update(sets[ch])
        if not sizes[c] or sizes[c] != by_id[c]["n_genes"]:
            raise ValueError("Node size mismatch")
        within[c] = cross[c] + sum(within[ch] for ch in children[c])
        if not children[c]:
            continue
        wp = sum(choose2(sizes[ch]) for ch in children[c])
        cp = choose2(sizes[c]) - wp
        ww = sum(within[ch] for ch in children[c])
        wd, cd = ww / wp if wp else None, cross[c] / cp if cp else None
        external = max(0.0, strength[c] - 2 * within[c])
        intersection = set.intersection(*(sets[ch] for ch in children[c]))
        child_densities = [within[ch] / choose2(sizes[ch]) for ch in children[c] if sizes[ch] > 1]
        result[c] = dict(
            n_genes=sizes[c],
            n_species=len(sets[c]),
            child_count=len(children[c]),
            min_child_fraction=min(sizes[ch] for ch in children[c]) / sizes[c],
            species_intersection_fraction=len(intersection) / len(sets[c]),
            child_species_sum_ratio=sum(len(sets[ch]) for ch in children[c]) / len(sets[c]),
            cross_density=cd,
            within_density=wd,
            cross_to_within_density=cd / wd if wd is not None and wd > 0 else None,
            child_within_density_min=min(child_densities) if child_densities else None,
            child_within_density_max=max(child_densities) if child_densities else None,
            child_without_internal_pairs=sum(sizes[ch] < 2 for ch in children[c]),
            cross_weight=cross[c],
            within_child_weight=ww,
            external_weight_fraction=external / (external + 2 * within[c])
            if external + 2 * within[c]
            else None,
            density_status="no_within_pairs"
            if not wp
            else "zero_within_weight"
            if not ww
            else "positive",
        )
    return result


def diagnose(run, reference_dir, baseline, out):
    prior = json.loads(baseline.read_text())
    hashes = dict(prior["input_hashes"])
    hashes[str(baseline)] = sha256_file(baseline)
    boundary_path = baseline.parent / "boundary-summary.json"
    hashes[str(boundary_path)] = sha256_file(boundary_path)
    boundary = json.loads(boundary_path.read_text())
    if boundary["source_run"] != str(run):
        raise ValueError("Baseline source run mismatch")
    hashes[str(run / "edges/retained_edges.parquet")] = boundary["fixed_graph_sha256"]

    def read(path):
        digest = sha256_file(path)
        if str(path) in hashes and hashes[str(path)] != digest:
            raise ValueError("Frozen input mismatch")
        hashes[str(path)] = digest
        return pq.ParquetFile(path).read().to_pylist()

    proteins = read(run / "input/proteins.parquet")
    species = {r["protein_id"]: r["species_id"] for r in proteins}
    keys = {(r["species_id"], r["original_id"]): r["protein_id"] for r in proteins}
    names = {r["species_name"]: r["species_id"] for r in read(run / "input/species.parquet")}
    reference = {
        keys[k]: v for k, v in read_of(reference_dir / "Orthogroups.tsv", names, keys).items()
    }
    with (run / "results/members.tsv").open() as handle:
        actual = {
            int(r["protein_id"]): r["family_id"] for r in csv.DictReader(handle, delimiter="\t")
        }
    index = read(run / "components/index.parquet")
    components = Counter(r["component_id"] for r in index)
    rows = []
    for cid, size in sorted(components.items()):
        if size == 1:
            continue
        folder = run / "hierarchy/components" / f"component={cid:08d}"
        nodes = read(folder / "nodes.parquet")
        membership = [
            (r["protein_id"], r["terminal_cluster_id"]) for r in read(folder / "members.parquet")
        ]
        by_id, children, order = topology(nodes)
        refs = {c: Counter() for c in by_id}
        for p, c in membership:
            if p in reference:
                refs[c][reference[p]] += 1
        for c in reversed(order):
            for ch in children[c]:
                refs[c].update(refs[ch])
        pure_losses = {}
        boundary_losses(nodes, membership, reference, actual, pure_losses)
        og = run / "orthogroups/components" / folder.name
        manifest_path = og / "og-manifest.json"
        hashes[str(manifest_path)] = sha256_file(manifest_path)
        manifest = json.loads(manifest_path.read_text())
        event_path = og / "v1_events.parquet"
        if (
            manifest["status"] != "DONE"
            or sha256_file(event_path) != manifest["output_checksums"]["v1_events.parquet"]
        ):
            raise ValueError("Production event identity mismatch")
        events = {
            r["cluster_id"]: r["v1_event"]
            for r in read(run / "orthogroups/components" / folder.name / "v1_events.parquet")
        }
        paths = sorted((run / "components/edges" / folder.name).glob("*.parquet"))
        if not paths:
            raise ValueError("Missing component edges")

        def edges(paths=paths):
            for path in paths:
                for e in read(path):
                    yield e["u"], e["v"], e["weight"]

        features = structural_features(nodes, membership, species, edges())
        for c, f in features.items():
            if events[c] not in ("II", "III-1", "III-2", "III-3"):
                continue
            labels = len(refs[c])
            rows.append(
                dict(
                    component_id=cid,
                    cluster_id=c,
                    event=events[c],
                    **f,
                    pure_lca_lost_pairs=pure_losses.get(c, {}).get("lost_pairs", 0),
                    diagnostic_cohort="pure_loss"
                    if c in pure_losses
                    else "pure_no_loss"
                    if labels == 1
                    else "mixed"
                    if labels > 1
                    else "unassigned_only",
                    assigned_size=sum(refs[c].values()),
                    reference_families=labels,
                )
            )
    metrics = [
        "min_child_fraction",
        "species_intersection_fraction",
        "child_species_sum_ratio",
        "cross_to_within_density",
        "external_weight_fraction",
        "child_within_density_min",
        "child_within_density_max",
    ]
    pure_loss_total = sum(r["pure_lca_lost_pairs"] for r in rows)
    if pure_loss_total != boundary["lca_losses"]["pure_assigned_lca_tree_available"]:
        raise ValueError("Pure loss cohort differs from frozen diagnostic total")
    grouped = defaultdict(list)
    for r in rows:
        size_bin = int(math.log2(r["n_genes"]))
        species_bin = int(math.log2(r["n_species"]))
        grouped[
            (
                r["diagnostic_cohort"],
                r["event"],
                size_bin,
                species_bin,
                "component0" if r["component_id"] == 0 else "other",
            )
        ].append(r)
    strata = []
    for key, group in sorted(grouped.items()):
        stats = {}
        for m in metrics:
            values = [r[m] for r in group if r[m] is not None]
            stats[m] = dict(
                n=len(values),
                missing=len(group) - len(values),
                q10_q50_q90=np.quantile(values, [0.1, 0.5, 0.9]).tolist() if values else None,
            )
        strata.append(
            dict(
                cohort=key[0],
                event=key[1],
                log2_size_bin=key[2],
                log2_species_bin=key[3],
                component_scope=key[4],
                nodes=len(group),
                pure_lca_lost_pairs=sum(r["pure_lca_lost_pairs"] for r in group),
                features=stats,
            )
        )
    for path, digest in hashes.items():
        if sha256_file(Path(path)) != digest:
            raise ValueError("Frozen inputs changed")
    out.mkdir(parents=True, exist_ok=False)
    (out / "nodes.json").write_text(json.dumps(rows, indent=2) + "\n")
    (out / "input-hashes.json").write_text(json.dumps(hashes, indent=2) + "\n")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                diagnostic_only=True,
                production_changed=False,
                immutable_inputs_verified=True,
                feature_reference_labels_used=False,
                cohorts=dict(Counter(r["diagnostic_cohort"] for r in rows)),
                pure_lca_lost_pairs=pure_loss_total,
                strata=strata,
                scope=(
                    "All excluded internal II/III nodes, including no-loss controls; "
                    "nested nodes not independent"
                ),
                source_run=str(run),
            ),
            indent=2,
        )
        + "\n"
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("run", "reference-dir", "baseline", "out"):
        p.add_argument("--" + key, type=Path, required=True)
    a = p.parse_args()
    diagnose(a.run, a.reference_dir, a.baseline, a.out)


if __name__ == "__main__":
    main()
