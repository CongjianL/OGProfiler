"""Read-only frozen hierarchy audit: depth paths, clade identity and shrink bounds."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import pyarrow.parquet as pq

from ogprofiler.core.manifest import sha256_file


def shrink_steps(n, fraction):
    """Integer worst-case accepted lineage to singleton; no search success claim."""
    f = Fraction(str(fraction))
    assert 0 < f < 1 and n >= 1
    steps = 0
    while n > 1:
        n = n * f.numerator // f.denominator
        assert n >= 1
        steps += 1
    return steps


def _tree(rows, members):
    nodes = {r["cluster_id"]: r for r in rows}
    assert len(nodes) == len(rows)
    paths = {}

    def path(cid):
        if cid not in paths:
            node = nodes[cid]
            paths[cid] = (() if node["parent_id"] is None else path(node["parent_id"])) + (cid,)
            assert node["depth"] == len(paths[cid]) - 1
        return paths[cid]

    descendants = defaultdict(set)
    terminal = dict(members)
    assert len(terminal) == len(members)
    for pid, leaf in members:
        for cid in path(leaf):
            descendants[cid].add(pid)
    assert all(len(descendants[c]) == n["n_genes"] for c, n in nodes.items())
    return nodes, paths, descendants, terminal


def audit_trees(
    rows, members, baseline_rows, baseline_members, *, max_depth=20, fraction=0.95, selected=None
):
    nodes, paths, descendants, terminal = _tree(rows, members)
    old, old_paths, old_desc, old_terminal = _tree(baseline_rows, baseline_members)
    assert set(terminal) == set(old_terminal)
    # Real clade identity, never cross-run numeric cluster IDs.
    old_clades = {frozenset(ps): cid for cid, ps in old_desc.items()}
    selected = selected or {}
    children, old_children = defaultdict(list), defaultdict(list)
    for source, index in ((rows, children), (baseline_rows, old_children)):
        for row in source:
            if row["parent_id"] is not None:
                index[row["parent_id"]].append(row["cluster_id"])
    partition_equal = {}

    def same_partition(cid, old_id):
        key = (cid, old_id)
        if key not in partition_equal:
            partition_equal[key] = {frozenset(descendants[c]) for c in children[cid]} == {
                frozenset(old_desc[c]) for c in old_children[old_id]
            }
        return partition_equal[key]

    failures = []
    for node in rows:
        if node["split_status"] != "UNRESOLVED":
            continue
        cid = node["cluster_id"]
        ids = descendants[cid]
        old_leaves = Counter(old_terminal[p] for p in ids)
        lineage = paths[cid]
        ancestors = []
        for aid in lineage[:-1]:
            a = nodes[aid]
            child = nodes[lineage[lineage.index(aid) + 1]]
            same = old_clades.get(frozenset(descendants[aid]))
            evidence = selected.get(aid)
            ancestors.append(
                dict(
                    cluster_id=aid,
                    depth=a["depth"],
                    n_genes=a["n_genes"],
                    child_count=a["child_count"],
                    selection_kind=a.get("selection_kind") or "KWAY",
                    child_on_path=child["cluster_id"],
                    child_size=child["n_genes"],
                    retained_fraction=child["n_genes"] / a["n_genes"],
                    baseline_exact_cluster=same,
                    baseline_exact_depth=old[same]["depth"] if same is not None else None,
                    baseline_exact_child_count=old[same]["child_count"]
                    if same is not None
                    else None,
                    baseline_exact_resolution=old[same].get("resolution")
                    if same is not None
                    else None,
                    baseline_partition_equal=same_partition(aid, same)
                    if same is not None
                    else None,
                    selected_evidence=evidence,
                )
            )
        first = next(iter(ids))
        lca_path = old_paths[old_terminal[first]]
        for p in ids:
            other = old_paths[old_terminal[p]]
            count = next(
                (i for i, (a, b) in enumerate(zip(lca_path, other, strict=False)) if a != b),
                min(len(lca_path), len(other)),
            )
            lca_path = lca_path[:count]
        lca = old[lca_path[-1]]
        same = old_clades.get(frozenset(ids))
        failures.append(
            dict(
                cluster_id=cid,
                depth=node["depth"],
                n_genes=node["n_genes"],
                n_species=node.get("n_species"),
                terminal_reason=node["terminal_reason"],
                selected_candidate_at_failed_node=bool(cid in selected),
                baseline_exact_cluster=same,
                exact_clade_depth_delta=node["depth"] - old[same]["depth"]
                if same is not None
                else None,
                baseline_lca=dict(
                    cluster_id=lca["cluster_id"],
                    depth=lca["depth"],
                    n_genes=lca["n_genes"],
                    purity=len(ids) / lca["n_genes"],
                ),
                baseline_terminal_depths=dict(
                    Counter(old[leaf]["depth"] for leaf in (old_terminal[p] for p in ids))
                ),
                baseline_terminal_reasons=dict(
                    Counter(old[old_terminal[p]]["terminal_reason"] for p in ids)
                ),
                baseline_leaf_intersections=dict(old_leaves),
                ancestor_selection_kinds=dict(Counter(a["selection_kind"] for a in ancestors)),
                first_partition_divergence=next(
                    (
                        dict(
                            depth=a["depth"],
                            cluster_id=a["cluster_id"],
                            selection_kind=a["selection_kind"],
                            n_genes=a["n_genes"],
                            new_child_count=a["child_count"],
                            baseline_child_count=a["baseline_exact_child_count"],
                        )
                        for a in ancestors
                        if a["baseline_partition_equal"] is False
                    ),
                    None,
                ),
                ancestors=ancestors,
                protein_ids=sorted(ids),
                accepted_shrink_tail_bound=shrink_steps(len(ids), fraction),
                depth_limit_free_bound=node["depth"] + shrink_steps(len(ids), fraction),
            )
        )
    proteins = sum(f["n_genes"] for f in failures)

    def weighted_depths(ns, mapping):
        return dict(sorted(Counter(ns[leaf]["depth"] for leaf in mapping.values()).items()))

    summary = dict(
        unique_complete_members=True,
        unresolved_nodes=len(failures),
        unresolved_proteins=proteins,
        all_depth_limit=all(
            f["terminal_reason"] == "DEPTH_LIMIT" and f["depth"] == max_depth for f in failures
        ),
        max_failed_size=max((f["n_genes"] for f in failures), default=0),
        failure_size_histogram=dict(sorted(Counter(f["n_genes"] for f in failures).items())),
        exact_baseline_clade_count=sum(f["baseline_exact_cluster"] is not None for f in failures),
        exact_clade_depth_deltas=dict(
            Counter(
                f["exact_clade_depth_delta"]
                for f in failures
                if f["exact_clade_depth_delta"] is not None
            )
        ),
        current_protein_terminal_depths=weighted_depths(nodes, terminal),
        baseline_protein_terminal_depths=weighted_depths(old, old_terminal),
        failure_ancestor_binary_counts=dict(
            Counter(f["ancestor_selection_kinds"].get("BINARY", 0) for f in failures)
        ),
        depth_limit_free_bound=max(
            (f["depth_limit_free_bound"] for f in failures), default=max_depth
        ),
        extra_calls_upper_bound=72 * (proteins - len(failures)),
        extra_nodes_upper_bound=2 * (proteins - len(failures)),
        first_divergence_selection_kinds=dict(
            Counter(
                f["first_partition_divergence"]["selection_kind"]
                for f in failures
                if f["first_partition_divergence"] is not None
            )
        ),
        failed_protein_baseline_depths=dict(
            Counter(old[old_terminal[p]]["depth"] for f in failures for p in f["protein_ids"])
        ),
        mean_current_depth=sum(nodes[leaf]["depth"] for leaf in terminal.values()) / len(terminal),
        mean_baseline_depth=sum(old[leaf]["depth"] for leaf in old_terminal.values())
        / len(old_terminal),
        old_node_count=len(old),
        new_node_count=len(nodes),
        old_max_depth=max(n["depth"] for n in old.values()),
        new_max_depth=max(n["depth"] for n in nodes.values()),
    )
    return dict(summary=summary, failed_nodes=failures)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--current", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    provenance = {}

    def load(root):
        folder = root / "new-hierarchy/hierarchy/components/component=00000000"
        manifest = json.loads((folder / "hierarchy-manifest.json").read_text())
        for name, digest in manifest["output_checksums"].items():
            assert sha256_file(folder / name) == digest
        provenance[str(root)] = manifest
        rows = pq.ParquetFile(folder / "nodes.parquet").read().to_pylist()
        members = pq.ParquetFile(folder / "members.parquet").read().to_pylist()
        candidates = pq.ParquetFile(folder / "candidates.parquet").read()
        selected = candidates.filter(candidates["selected"]).to_pylist()
        return rows, [(r["protein_id"], r["terminal_cluster_id"]) for r in members], selected

    rows, members, chosen = load(args.current)
    old, old_members, old_chosen = load(args.baseline)
    new_manifest = provenance[str(args.current)]
    old_manifest = provenance[str(args.baseline)]
    assert new_manifest["input_checksums"] == old_manifest["input_checksums"] or all(
        old_manifest["input_checksums"].get(p) == h
        for p, h in new_manifest["input_checksums"].items()
        if p != "run.yaml"
    )
    new_cfg = new_manifest["parameters"]["hierarchy"]
    old_cfg = json.loads(json.dumps(old_manifest["parameters"]["hierarchy"]))
    old_cfg["resolution"].setdefault("topology_policy", "kway_v1")
    old_cfg["resolution"]["topology_policy"] = "soft_binary_24_v2"
    assert new_cfg == old_cfg
    report = audit_trees(
        rows,
        members,
        old,
        old_members,
        max_depth=new_cfg["max_depth"],
        fraction=new_cfg["resolution"]["max_child_fraction"],
        selected={c["cluster_id"]: c for c in chosen},
    )
    failures = {f["cluster_id"] for f in report["failed_nodes"]}
    table = pq.ParquetFile(
        args.current / "new-hierarchy/hierarchy/components/component=00000000/candidates.parquet"
    ).read(columns=["cluster_id"])
    report["summary"]["failed_nodes_have_zero_candidates"] = not any(
        cid in failures for cid in table["cluster_id"].to_pylist()
    )

    def selected_summary(cs):
        return dict(
            count=len(cs),
            max_child_fraction=max(c["max_child_fraction"] for c in cs),
            min_stability=min(c["stability"] for c in cs),
            all_valid=all(c["valid"] and not c["violations"] for c in cs),
        )

    report["selected_gates"] = dict(
        current=selected_summary(chosen), baseline=selected_summary(old_chosen)
    )
    report["metrics"] = {
        str(root): json.loads(
            (
                root / "new-hierarchy/hierarchy/components/component=00000000/metrics.json"
            ).read_text()
        )
        for root in (args.current, args.baseline)
    }
    report["source_hash"] = sha256_file(Path(__file__))
    report["manifests"] = provenance
    args.out.write_text(json.dumps(report, indent=2))
    print(json.dumps(report["summary"], indent=2))


if __name__ == "__main__":
    main()
