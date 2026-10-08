"""Read-only induced-subtree merge/cut objective; no production partition changes."""

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
from benchmarks.og_extraction.fixed_tree_graph import graph_cut
from benchmarks.og_extraction.tree_cut import topology
from ogprofiler.core.manifest import sha256_file


def diagnostic_pairs(reference, prediction):
    counts = Counter(reference.values())
    groups = Counter(prediction[p] for p in reference)
    joint = Counter((reference[p], prediction[p]) for p in reference)
    tp = sum(map(choose2, joint.values()))
    pp = sum(map(choose2, groups.values()))
    return dict(tp=tp, fp=pp - tp, lost=sum(map(choose2, counts.values())) - tp)


def inspect(nodes, membership, edges, *, species=None):
    if species is None:
        cut = graph_cut(nodes, membership, edges, strength_null=True)
    else:
        from benchmarks.og_extraction.species_pair_graph import species_pair_cut

        cut = species_pair_cut(nodes, membership, edges, species)
    root = next(n["cluster_id"] for n in nodes if n["parent_id"] is None)
    root_score = next(n["keep_score"] for n in cut["node_scores"] if n["cluster_id"] == root)
    return dict(
        merge_score=root_score,
        best_cut_score=cut["score"],
        cut_gain=cut["score"] - root_score,
        selected=cut["selected"],
        split=root not in cut["selected"],
        total_weight=cut["component_total_weight"],
        groups=cut["diagnostics"]["groups"],
    ), {r["protein_id"]: r["cluster_id"] for r in cut["members"]}


def diagnose(run, reference_dir, trace_dir, out, *, condition_species_pairs=False):
    hashes = json.loads((trace_dir / "input-hashes.json").read_text())
    for name in ("summary.json", "groups.json", "input-hashes.json"):
        path = trace_dir / name
        hashes[str(path)] = sha256_file(path)
    trace = json.loads((trace_dir / "summary.json").read_text())
    if not trace["replay_verified"] or not trace["immutable_inputs_verified"]:
        raise ValueError("Unverified source replay")
    rows = json.loads((trace_dir / "groups.json").read_text())
    if any(r["scope"] != "complete_source" or not r["msa_admitted"] for r in rows):
        raise ValueError("Diagnostic targets no longer disjoint complete admitted sources")

    def read(path):
        digest = sha256_file(path)
        if str(path) in hashes and hashes[str(path)] != digest:
            raise ValueError("Frozen identity mismatch")
        hashes[str(path)] = digest
        return pq.ParquetFile(path).read().to_pylist()

    proteins = read(run / "input/proteins.parquet")
    species = {r["protein_id"]: r["species_id"] for r in proteins}
    names = {r["species_name"]: r["species_id"] for r in read(run / "input/species.parquet")}
    keys = {(r["species_id"], r["original_id"]): r["protein_id"] for r in proteins}
    with (run / "results/members.tsv").open() as h:
        actual = {int(r["protein_id"]): r["family_id"] for r in csv.DictReader(h, delimiter="\t")}
    targets = defaultdict(dict)
    for r in rows:
        targets[r["component_id"]][r["source_cluster_id"]] = r
    cases = []
    cut_files = []
    seen_proteins = set()
    for cid, source_rows in sorted(targets.items()):
        folder = run / "hierarchy/components" / f"component={cid:08d}"
        nodes = read(folder / "nodes.parquet")
        by_id, children, order = topology(nodes)
        owner = {}
        subnodes = defaultdict(list)
        members = defaultdict(list)
        edge_groups = defaultdict(list)
        for c in order:
            parent = by_id[c]["parent_id"]
            if c in source_rows and owner.get(parent) is not None:
                raise ValueError("Nested selected target")
            owner[c] = c if c in source_rows else owner.get(parent)
            if owner[c] is not None:
                subnodes[owner[c]].append(
                    dict(by_id[c], parent_id=None if c == owner[c] else parent)
                )
        protein_owner = {}
        for m in read(folder / "members.parquet"):
            p, c = m["protein_id"], m["terminal_cluster_id"]
            if owner[c] is not None:
                members[owner[c]].append((p, c))
                protein_owner[p] = owner[c]
                if p in seen_proteins:
                    raise ValueError("Duplicate target protein")
                seen_proteins.add(p)
        paths = sorted((run / "components/edges" / folder.name).glob("*.parquet"))
        if not paths:
            raise ValueError("Missing component edges")
        for path in paths:
            for e in read(path):
                a, b = protein_owner.get(e["u"]), protein_owner.get(e["v"])
                if a is not None and a == b:
                    edge_groups[a].append((e["u"], e["v"], e["weight"]))
        # Selection uses full proteins and induced weights, never reference labels.
        for c, r in source_rows.items():
            baseline_result, baseline_prediction = inspect(subnodes[c], members[c], edge_groups[c])
            result, prediction = (
                inspect(subnodes[c], members[c], edge_groups[c], species=species)
                if condition_species_pairs
                else (baseline_result, baseline_prediction)
            )
            if condition_species_pairs:
                from benchmarks.qfo.species_pair_null import conditioned_null

                saved_gain = conditioned_null(baseline_prediction, species, edge_groups[c])["gain"]
                result["unconditional_cut_under_conditioned_objective"] = saved_gain
                result["optimization_gain_over_unconditional_cut"] = (
                    result["best_cut_score"] - saved_gain
                )
                if result["best_cut_score"] < saved_gain and not math.isclose(
                    result["best_cut_score"], saved_gain, abs_tol=1e-9, rel_tol=1e-8
                ):
                    raise ValueError("Conditioned DP worse than feasible unconditional cut")
                result["objective_identity_verified"] = True
            if len(prediction) != r["source_full_size"]:
                raise ValueError("Target source coverage mismatch")
            cases.append(
                dict(
                    component_id=cid,
                    source_cluster_id=c,
                    original_event=r["original_event"],
                    full_size=len(prediction),
                    merge_transitions=r["output_transitions"],
                    objective=result,
                    unconditional_objective=baseline_result,
                )
            )
            cut_files.append(
                dict(
                    component_id=cid,
                    source_cluster_id=c,
                    prediction=prediction,
                    unconditional_prediction=baseline_prediction,
                )
            )
    # Reference labels are first loaded only after all cuts have been fixed.
    path = reference_dir / "Orthogroups.tsv"
    digest = sha256_file(path)
    if str(path) not in hashes or hashes[str(path)] != digest:
        raise ValueError("Reference identity mismatch")
    reference = {keys[k]: v for k, v in read_of(path, names, keys).items()}
    for case, cut in zip(cases, cut_files, strict=True):
        pred = cut["prediction"]
        ref = {p: reference[p] for p in pred if p in reference}
        merged = diagnostic_pairs(ref, {p: 0 for p in pred})
        if any(merged[k] != case["merge_transitions"][k] for k in ("tp", "fp")):
            raise ValueError("Source merge scoring differs from frozen trace")
        score = diagnostic_pairs(ref, pred)
        retained = Counter((actual[p], ref[p], pred[p]) for p in ref)
        shared_tp = sum(map(choose2, retained.values()))
        shared_pairs = sum(map(choose2, Counter((actual[p], pred[p]) for p in ref).values()))
        score["new_tp"] = score["tp"] - shared_tp
        score["new_fp"] = score["fp"] - (shared_pairs - shared_tp)
        case["cut_pairs"] = score
        case["unconditional_cut_pairs"] = diagnostic_pairs(ref, cut["unconditional_prediction"])
        case["conditioned_minus_unconditional_tp"] = (
            score["tp"] - case["unconditional_cut_pairs"]["tp"]
        )
        case["conditioned_minus_unconditional_fp"] = (
            score["fp"] - case["unconditional_cut_pairs"]["fp"]
        )
        case["reference_cohort"] = (
            "polluted" if case["merge_transitions"]["new_fp"] else "clean_tp_gain"
        )
        case["removed_merge_tp"] = case["merge_transitions"]["tp"] - score["tp"]
        case["removed_merge_fp"] = case["merge_transitions"]["fp"] - score["fp"]
    strata = []
    grouped = defaultdict(list)
    for c in cases:
        grouped[(c["reference_cohort"], c["original_event"])].append(c)
    for (cohort, event), group in sorted(grouped.items()):
        strata.append(
            dict(
                cohort=cohort,
                event=event,
                nodes=len(group),
                conditioned_minus_unconditional_tp=sum(
                    c["conditioned_minus_unconditional_tp"] for c in group
                ),
                conditioned_minus_unconditional_fp=sum(
                    c["conditioned_minus_unconditional_fp"] for c in group
                ),
                unconditional_split_nodes=sum(c["unconditional_objective"]["split"] for c in group),
                split_nodes=sum(c["objective"]["split"] for c in group),
                cut_gain_q10_q50_q90=np.quantile(
                    [c["objective"]["cut_gain"] for c in group], [0.1, 0.5, 0.9]
                ).tolist(),
                removed_merge_tp=sum(c["removed_merge_tp"] for c in group),
                removed_merge_fp=sum(c["removed_merge_fp"] for c in group),
                retained_new_tp=sum(c["cut_pairs"]["new_tp"] for c in group),
                retained_new_fp=sum(c["cut_pairs"]["new_fp"] for c in group),
            )
        )
    for path, digest in hashes.items():
        if sha256_file(Path(path)) != digest:
            raise ValueError("Frozen input changed")
    out.mkdir(parents=True, exist_ok=False)
    (out / "cases.json").write_text(json.dumps(cases, indent=2) + "\n")
    # Independent diagnostic cuts, not a replacement full production partition.
    (out / "cuts-DIAGNOSTIC-ONLY.json").write_text(json.dumps(cut_files) + "\n")
    (out / "input-hashes.json").write_text(json.dumps(hashes, indent=2) + "\n")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                diagnostic_only=True,
                production_changed=False,
                immutable_inputs_verified=True,
                reference_labels_used_for_objective=False,
                species_pair_conditioned=condition_species_pairs,
                objective=(
                    "sum_v (Win_v - sum_s k_v,ss^2/(4W_ss) - sum_s<t k_v,st*k_v,ts/W_st)/W; "
                    "source-fixed denominators, zero blocks omitted, resolution1, tie keep parent"
                    if condition_species_pairs
                    else "sum Win/W - (strength/(2W))^2; induced source graph, "
                    "resolution1, tie keep parent"
                ),
                baseline_merge_score="root zero (or zero-weight case)",
                target_nodes=len(cases),
                coverage_verified=True,
                strata=strata,
                scope=(
                    "Only changed complete-source clades; "
                    "no full partition strategy or biological inference"
                ),
                top_polluted=sorted(
                    [c for c in cases if c["reference_cohort"] == "polluted"],
                    key=lambda c: -c["merge_transitions"]["new_fp"],
                )[:10],
                top_clean=sorted(
                    [c for c in cases if c["reference_cohort"] == "clean_tp_gain"],
                    key=lambda c: -c["merge_transitions"]["new_tp"],
                )[:10],
            ),
            indent=2,
        )
        + "\n"
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("run", "reference-dir", "trace-dir", "out"):
        p.add_argument("--" + key, type=Path, required=True)
    p.add_argument("--condition-species-pairs", action="store_true")
    a = p.parse_args()
    diagnose(
        a.run,
        a.reference_dir,
        a.trace_dir,
        a.out,
        condition_species_pairs=a.condition_species_pairs,
    )


if __name__ == "__main__":
    main()
