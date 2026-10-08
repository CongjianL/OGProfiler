"""Read-only all-group contamination delta and path audit on frozen artifacts.

Extra genes are benchmark-relative, not proof of biological non-orthology.
Global OG IDs are never compared between runs; original protein IDs are.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from benchmarks.og_extraction.refog_split_audit import load_component, rows
from ogprofiler.core.manifest import sha256_file


def common_path(a, b):
    result = []
    for x, y in zip(a, b, strict=False):
        if x != y:
            break
        result.append(x)
    if not result:
        raise ValueError("Pair has no common hierarchy root")
    return result


def export_members(root):
    groups, gene_group = defaultdict(set), {}
    with (root / "new-hierarchy/results/members.tsv").open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            gene = row["original_id"]
            assert gene not in gene_group
            groups[row["family_id"]].add(gene)
            gene_group[gene] = row["family_id"]
    return groups, gene_group


def audit(current, baseline):
    for relative in (
        "new-hierarchy/input/proteins.parquet",
        "new-hierarchy/components/index.parquet",
    ):
        assert sha256_file(current / relative) == sha256_file(baseline / relative)
    for ref in (current / "benchmark/RefOGs").rglob("*.txt"):
        assert sha256_file(ref) == sha256_file(baseline / ref.relative_to(current))
    groups, predicted = export_members(current)
    _, old_predicted = export_members(baseline)
    diag = list(
        csv.DictReader(
            (
                current / "h5-metrics/new_bounded_hierarchy_v1_compatible/refog-diagnostics.tsv"
            ).open(),
            delimiter="\t",
        )
    )
    proteins = {
        p["original_id"]: p["protein_id"]
        for p in rows(current / "new-hierarchy/input/proteins.parquet")
    }
    components = {
        p["protein_id"]: p["component_id"]
        for p in rows(current / "new-hierarchy/components/index.parquet")
    }
    cases = []
    views = {}
    old_views = {}
    for d in diag:
        name = d["refog"]
        raw = set((current / f"benchmark/RefOGs/{name}.txt").read_text().splitlines())
        low_path = current / f"benchmark/RefOGs/low_certainty_assignments/{name}.txt"
        low = set(low_path.read_text().splitlines()) if low_path.exists() else set()
        truth = raw - low
        for global_group in sorted({predicted[g] for g in truth}):
            members = groups[global_group] - low
            extras = members - truth
            overlap = members & truth
            if global_group == d["best_group"]:
                assert len(extras) == int(d["best_group_extra_genes"])
            pairs = [
                (a, b)
                for a in sorted(overlap)
                for b in sorted(extras)
                if old_predicted[a] != old_predicted[b]
            ]
            if not pairs:
                continue
            ids = {proteins[g] for g in members}
            cids = {components[p] for p in ids}
            assert len(cids) == 1
            cid = next(iter(cids))
            # Load all protein paths only once per relevant component; validates manifests.
            if cid not in views:
                hr = current / f"new-hierarchy/hierarchy/components/component={cid:08d}"
                all_ids = [r["protein_id"] for r in rows(hr / "members.parquet")]
                views[cid] = load_component(current / "new-hierarchy", cid, all_ids)
            view = views[cid]
            if cid not in old_views:
                old_hr = baseline / f"new-hierarchy/hierarchy/components/component={cid:08d}"
                old_ids = [r["protein_id"] for r in rows(old_hr / "members.parquet")]
                old_views[cid] = load_component(baseline / "new-hierarchy", cid, old_ids)
            old_view = old_views[cid]
            assert set(view["membership"]) == set(old_view["membership"])
            ogdir = current / f"new-hierarchy/orthogroups/components/component={cid:08d}"
            membership = rows(ogdir / "members.parquet")
            local_ids = {r["local_group_id"] for r in membership if r["protein_id"] in ids}
            assert len(local_ids) == 1
            local_id = next(iter(local_ids))
            assert {r["protein_id"] for r in membership if r["local_group_id"] == local_id} == {
                proteins[g] for g in groups[global_group]
            }
            group = next(
                r for r in rows(ogdir / "groups.parquet") if r["local_group_id"] == local_id
            )
            source = group["source_cluster_id"]
            source_members = {p for p, path in view["paths"].items() if source in path}
            old_clades = defaultdict(set)
            for pid, path in old_view["paths"].items():
                for node in path:
                    old_clades[node].add(pid)
            new_clades = defaultdict(set)
            for pid, path in view["paths"].items():
                for node in path:
                    new_clades[node].add(pid)

            def partitions(v, clades):
                children = defaultdict(list)
                for n in v["nodes"].values():
                    if n["parent_id"] is not None:
                        children[n["parent_id"]].append(n["cluster_id"])
                return {
                    node: frozenset(frozenset(clades[c]) for c in cs)
                    for node, cs in children.items()
                }

            new_parts, old_parts = partitions(view, new_clades), partitions(old_view, old_clades)
            matches = [node for node, clade in old_clades.items() if clade == source_members]
            baseline_source_clades = [
                dict(
                    node=old_view["nodes"][node],
                    event=old_view["events"][node]["v1_event"],
                    trace=old_view["traces"][node],
                )
                for node in matches
            ]
            records = []
            for left, right in pairs:
                a, b = view["paths"][proteins[left]], view["paths"][proteins[right]]
                lca = common_path(a, b)[-1]
                assert source in common_path(a, b)
                divergence = None
                for new_node, old_node in zip(a, old_view["paths"][proteins[left]], strict=False):
                    assert new_clades[new_node] == old_clades[old_node]
                    if new_parts.get(new_node) != old_parts.get(old_node):
                        divergence = dict(
                            current_node=view["nodes"][new_node],
                            current_event=view["events"][new_node]["v1_event"],
                            baseline_node=old_view["nodes"][old_node],
                            baseline_event=old_view["events"][old_node]["v1_event"],
                        )
                        break
                assert divergence is not None
                old_lca = common_path(
                    old_view["paths"][proteins[left]], old_view["paths"][proteins[right]]
                )[-1]
                records.append(
                    dict(
                        first_partition_divergence=divergence,
                        baseline_lca=old_lca,
                        baseline_lca_node=old_view["nodes"][old_lca],
                        baseline_lca_event=old_view["events"][old_lca]["v1_event"],
                        left=left,
                        extra=right,
                        baseline_left_og=old_predicted[left],
                        baseline_extra_og=old_predicted[right],
                        lca=lca,
                        lca_node=view["nodes"][lca],
                        lca_event=view["events"][lca]["v1_event"],
                        mechanism="SAME_TERMINAL"
                        if a[-1] == b[-1]
                        else "OG_REJOINS_DISTINCT_TERMINALS",
                    )
                )
            cases.append(
                dict(
                    refog=name,
                    component_id=cid,
                    group_id=global_group,
                    is_best_group=global_group == d["best_group"],
                    classification=d["classification"],
                    extra_genes=sorted(extras),
                    newly_joined_pairs=len(pairs),
                    group=group,
                    baseline_equivalent_source_clades=baseline_source_clades,
                    selected_source_node=view["nodes"][source],
                    selected_source_candidate=view["selected"].get(source),
                    source_trace=view["traces"][source],
                    pairs=records,
                )
            )
    return dict(
        read_only=True,
        scope=(
            "all current RefOG-overlapping groups; "
            "newly joined baseline-separated cross-truth pairs; "
            "not official precision decomposition"
        ),
        baseline_job="1410777",
        current_job="1410868",
        source_sha256=sha256_file(Path(__file__)),
        export_sha256=sha256_file(current / "new-hierarchy/results/members.tsv"),
        baseline_export_sha256=sha256_file(baseline / "new-hierarchy/results/members.tsv"),
        cases=cases,
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--current", type=Path, required=True)
    p.add_argument("--baseline", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    report = audit(a.current, a.baseline)
    a.out.write_text(json.dumps(report, indent=2) + "\n")
    counts = defaultdict(int)
    for c in report["cases"]:
        counts[c["refog"]] += c["newly_joined_pairs"]
    print(json.dumps(counts, indent=2))


if __name__ == "__main__":
    main()
