"""Read-only fixed-tree boundary opportunities; all reference-driven results are diagnostic."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import pyarrow.parquet as pq

from benchmarks.og_extraction.embleya_reference import choose2, compare, read_of
from benchmarks.og_extraction.embleya_representation import run_audit
from benchmarks.og_extraction.tree_cut import topology
from ogprofiler.core.manifest import sha256_file


def boundary_losses(nodes, membership, reference, actual, pure_losses=None):
    """Count lost same-reference pairs at their unique LCA, without enumerating pairs.

    A pure assigned-reference LCA is a tree-available grouping opportunity.
    A mixed LCA entails contamination when kept as a whole node. Neither
    category by itself proves the cause of the production selection.
    """
    by_id, children, order = topology(nodes)
    refs = {c: Counter() for c in by_id}
    joints = {c: Counter() for c in by_id}
    seen = set()
    for p, leaf in membership:
        if p in seen or leaf not in by_id or children[leaf] or p not in actual:
            raise ValueError("Invalid terminal membership/prediction")
        seen.add(p)
        if p in reference:
            refs[leaf][reference[p]] += 1
            joints[leaf][(reference[p], actual[p])] += 1
    tp, retained = {}, {}
    result = Counter(
        {
            k: 0
            for k in (
                "within_structural_leaf",
                "pure_assigned_lca_tree_available",
                "mixed_assigned_lca_requires_contamination",
            )
        }
    )
    for c in reversed(order):
        for child in children[c]:
            refs[c].update(refs[child])
            joints[c].update(joints[child])
        tp[c] = sum(map(choose2, refs[c].values()))
        retained[c] = sum(map(choose2, joints[c].values()))
        lost = (
            tp[c]
            - sum(tp[u] for u in children[c])
            - retained[c]
            + sum(retained[u] for u in children[c])
        )
        if lost < 0:
            raise ValueError("Negative LCA loss")
        category = (
            "within_structural_leaf"
            if not children[c]
            else "pure_assigned_lca_tree_available"
            if len(refs[c]) == 1
            else "mixed_assigned_lca_requires_contamination"
        )
        result[category] += lost
        if pure_losses is not None and category == "pure_assigned_lca_tree_available" and lost:
            pure_losses[c] = dict(
                lost_pairs=lost,
                reference_og=next(iter(refs[c])),
                assigned_size=sum(refs[c].values()),
                production_fragments=len(joints[c]),
            )
    result["within_component_lost_pairs"] = tp[order[0]] - retained[order[0]]
    if (
        sum(
            result[k]
            for k in (
                "within_structural_leaf",
                "pure_assigned_lca_tree_available",
                "mixed_assigned_lca_requires_contamination",
            )
        )
        != result["within_component_lost_pairs"]
    ):
        raise ValueError("LCA categories do not sum")
    return dict(result)


def diagnose(run, reference_dir, out):
    # Existing interface audits hashes, full coverage and global fractional oracle.
    graph = run / "edges/retained_edges.parquet"
    extra_hashes = {str(graph): sha256_file(graph)}
    run_audit(run, reference_dir, out)
    report = json.loads((out / "report.json").read_text())
    proteins = pq.ParquetFile(run / "input/proteins.parquet").read().to_pylist()
    keys = {(r["species_id"], r["original_id"]): r["protein_id"] for r in proteins}
    names = {
        r["species_name"]: r["species_id"]
        for r in pq.ParquetFile(run / "input/species.parquet").read().to_pylist()
    }
    ref = {keys[k]: v for k, v in read_of(reference_dir / "Orthogroups.tsv", names, keys).items()}
    with (run / "results/members.tsv").open() as h:
        actual = {int(r["protein_id"]): r["family_id"] for r in csv.DictReader(h, delimiter="\t")}
    index = pq.ParquetFile(run / "components/index.parquet").read().to_pylist()
    components = {r["protein_id"]: r["component_id"] for r in index}
    by_component = {}
    for p, c in components.items():
        by_component.setdefault(c, []).append(p)
    terminal, totals, details = {}, Counter(), []
    for cid, ids in sorted(by_component.items()):
        folder = run / "hierarchy/components" / f"component={cid:08d}"
        if len(ids) == 1 and not (folder / "nodes.parquet").exists():
            nodes = [dict(cluster_id=0, parent_id=None)]
            membership = [(ids[0], 0)]
        else:
            nodes = pq.ParquetFile(folder / "nodes.parquet").read().to_pylist()
            membership = [
                (r["protein_id"], r["terminal_cluster_id"])
                for r in pq.ParquetFile(folder / "members.parquet").read().to_pylist()
            ]
        loss = boundary_losses(nodes, membership, ref, actual)
        totals.update(loss)
        details.append(dict(component_id=cid, **loss))
        terminal.update({p: (cid, leaf) for p, leaf in membership})
    if set(terminal) != set(components):
        raise ValueError("Terminal cut coverage mismatch")
    terminal_score = compare(ref, terminal, components)
    cross_lost = report["actual"]["lost_pairs"] - totals["within_component_lost_pairs"]
    if cross_lost < 0:
        raise ValueError("Negative cross-component loss")
    summary = dict(
        diagnostic_only=True,
        production_changed=False,
        source_run=str(run),
        reference_dir=str(reference_dir),
        fixed_graph_sha256=extra_hashes[str(graph)],
        audit_report_sha256=sha256_file(out / "report.json"),
        diagnostic_source_sha256=sha256_file(Path(__file__)),
        attribution_scope=(
            "Actual lost reference pairs by LCA; opportunities/constraints, "
            "not causal policy attribution"
        ),
        purity_scope="OF3 assigned universe only; unassigned excluded",
        cross_component_lost_pairs=cross_lost,
        lca_losses=dict(totals),
        actual=report["actual"],
        terminal_cut=terminal_score,
        unrestricted_oracle=report["oracle_cut"],
        eligible_oracle=report["eligible_oracle_cut"],
        cut_family_comparison=report["cut_family_comparison"],
        oracle_scope=(
            "Complete fixed-tree cuts; not an upper bound for arbitrary active-view OG extraction"
        ),
        representation_summary={
            k: report[k]
            for k in [
                "exact_any_node",
                "exact_eligible_node",
                "positive_eligibility_gap",
                "positive_selection_gap",
                "whole_hierarchy_macro_best_f1",
                "eligible_macro_best_f1",
                "actual_macro_best_f1",
            ]
        },
        top_components=sorted(details, key=lambda r: -r["within_component_lost_pairs"])[:30],
    )
    for p, digest in {**report["input_hashes"], **extra_hashes}.items():
        if sha256_file(Path(p)) != digest:
            raise ValueError("Frozen source changed during diagnostics")
    summary["immutable_inputs_verified"] = True
    (out / "boundary-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("run", "reference-dir", "out"):
        p.add_argument("--" + key, type=Path, required=True)
    a = p.parse_args()
    diagnose(a.run, a.reference_dir, a.out)


if __name__ == "__main__":
    main()
