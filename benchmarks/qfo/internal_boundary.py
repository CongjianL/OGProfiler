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


def inspect(nodes, membership, edges, *, species=None, exact=False):
    if species is None:
        cut = graph_cut(nodes, membership, edges, strength_null=True)
    else:
        from benchmarks.og_extraction.species_pair_graph import species_pair_cut

        cut = species_pair_cut(nodes, membership, edges, species, exact=exact)
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


def diagnose(
    run,
    reference_dir,
    trace_dir,
    out,
    *,
    condition_species_pairs=False,
    tie_replay_dir=None,
    exact_conditioned=False,
    raw_replay_dir=None,
    singleton_flow=False,
    block_degeneracy=False,
    context_flow=False,
):
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

    if context_flow and not block_degeneracy:
        raise ValueError("Context flow requires fixed legal-node block profiles")
    if block_degeneracy and not exact_conditioned:
        raise ValueError("Degeneracy requires fixed exact conditioned cuts")
    if singleton_flow and tie_replay_dir is None:
        raise ValueError("Singleton flow requires fixed focused tie replay")
    proteins = read(run / "input/proteins.parquet")
    protein_metadata = {r["protein_id"]: r for r in proteins}
    species = {r["protein_id"]: r["species_id"] for r in proteins}
    names = {r["species_name"]: r["species_id"] for r in read(run / "input/species.parquet")}
    keys = {(r["species_id"], r["original_id"]): r["protein_id"] for r in proteins}
    with (run / "results/members.tsv").open() as h:
        actual = {int(r["protein_id"]): r["family_id"] for r in csv.DictReader(h, delimiter="\t")}
    replay_cuts = {}
    if raw_replay_dir is not None and not exact_conditioned:
        raise ValueError("Raw replay requires exact conditioned mode")
    if tie_replay_dir is not None and exact_conditioned:
        raise ValueError("Tie audit must replay historical raw mode")
    replay_dir = tie_replay_dir or raw_replay_dir
    if exact_conditioned and not condition_species_pairs:
        raise ValueError("Exact mode requires conditioned objective")
    if replay_dir is not None:
        if not condition_species_pairs:
            raise ValueError("Tie audit requires conditioned objective")
        for name in ("summary.json", "cuts-DIAGNOSTIC-ONLY.json", "input-hashes.json"):
            hashes[str(replay_dir / name)] = sha256_file(replay_dir / name)
        replay = json.loads((replay_dir / "summary.json").read_text())
        if not replay["coverage_verified"] or not replay["immutable_inputs_verified"]:
            raise ValueError("Unverified tie replay")
        for r in json.loads((replay_dir / "cuts-DIAGNOSTIC-ONLY.json").read_text()):
            replay_cuts[r["component_id"], r["source_cluster_id"]] = {
                int(p): c for p, c in r["prediction"].items()
            }
        if tie_replay_dir is not None:
            rows = [
                r
                for r in rows
                if (r["component_id"], r["source_cluster_id"]) in {(0, 10040), (0, 21394)}
            ]
            if len(rows) != 2:
                raise ValueError("Missing focused tie targets")
    targets = defaultdict(dict)
    for r in rows:
        targets[r["component_id"]][r["source_cluster_id"]] = r
    cases = []
    cut_files = []
    seen_proteins = set()
    degeneracy_inputs = {}
    for cid, source_rows in sorted(targets.items()):
        folder = run / "hierarchy/components" / f"component={cid:08d}"
        nodes = read(folder / "nodes.parquet")
        by_id, children, order = topology(nodes)
        owner = {}
        subnodes = defaultdict(list)
        members = defaultdict(list)
        edge_groups = defaultdict(list)
        singleton_incident = []
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
                if (
                    singleton_flow
                    and cid == 0
                    and (e["u"] in {31672, 31901} or e["v"] in {31672, 31901})
                ):
                    singleton_incident.append((e["u"], e["v"], e["weight"]))
                a, b = protein_owner.get(e["u"]), protein_owner.get(e["v"])
                if a is not None and a == b:
                    edge_groups[a].append((e["u"], e["v"], e["weight"]))
        # Selection uses full proteins and induced weights, never reference labels.
        for c, r in source_rows.items():
            baseline_result, baseline_prediction = inspect(subnodes[c], members[c], edge_groups[c])
            result, prediction = (
                inspect(
                    subnodes[c],
                    members[c],
                    edge_groups[c],
                    species=species,
                    exact=exact_conditioned,
                )
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
            if exact_conditioned:
                raw_result, raw_prediction = inspect(
                    subnodes[c], members[c], edge_groups[c], species=species, exact=False
                )
                result["raw_conditioned_objective"] = raw_result
                if raw_replay_dir is not None and raw_prediction != replay_cuts[cid, c]:
                    raise ValueError("Frozen raw conditioned replay differs")
            else:
                raw_prediction = prediction
            if tie_replay_dir is not None:
                from benchmarks.og_extraction.conditioned_tie_audit import audit_tree

                if prediction != replay_cuts[cid, c]:
                    raise ValueError("Frozen conditioned cut replay differs")
                audit = audit_tree(subnodes[c], members[c], edge_groups[c], species)
                if audit["raw_prediction"] != prediction:
                    raise ValueError("Local raw path replay differs")
                result["tie_audit"] = audit
            if singleton_flow and cid == 0 and c == 21394:
                from benchmarks.og_extraction.singleton_edge_flow import trace_flow

                flow = trace_flow(
                    subnodes[c], members[c], edge_groups[c], singleton_incident, species, 21396
                )
                if set(flow["singletons"]) != {31672, 31901}:
                    raise ValueError("Frozen singleton identity differs")
                for edge in flow["incident_edges"]:
                    edge["singleton_original_id"] = protein_metadata[edge["singleton"]][
                        "original_id"
                    ]
                    edge["other_original_id"] = protein_metadata[edge["other"]]["original_id"]
                split_gain = flow["candidate_partitions"][1]["gain_vs_keep"]
                audit_row = next(
                    row for row in result["tie_audit"]["local_rows"] if row["cluster_id"] == 21396
                )
                if not math.isclose(
                    split_gain, audit_row["exact_split_minus_keep"], rel_tol=1e-12, abs_tol=1e-15
                ):
                    raise ValueError("Raw child-pair flow and saved local objective differ")
                flow["local_objective_identity_verified"] = True
                result["singleton_flow"] = flow
            if block_degeneracy:
                from benchmarks.og_extraction.species_block_degeneracy import profile

                result["block_degeneracy"] = profile(
                    subnodes[c],
                    members[c],
                    edge_groups[c],
                    species,
                    set(result["selected"]),
                    context_flow=context_flow,
                )
                degeneracy_inputs[cid, c] = (subnodes[c], members[c])
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
                    raw_conditioned_prediction=raw_prediction,
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
        if tie_replay_dir is not None:
            audit = case["objective"]["tie_audit"]
            audit["exact_cut_pairs"] = diagnostic_pairs(ref, audit["exact_prediction"])
            for row in audit["local_rows"]:
                local_ref = {p: ref[p] for p in row["proteins"] if p in ref}
                row["keep_pairs"] = diagnostic_pairs(local_ref, dict.fromkeys(local_ref, 0))
                if row["raw_child_prediction"]:
                    row["raw_child_pairs"] = diagnostic_pairs(
                        local_ref, row["raw_child_prediction"]
                    )
                    row["exact_child_pairs"] = diagnostic_pairs(
                        local_ref, row["exact_child_prediction"]
                    )
        if "singleton_flow" in case["objective"]:
            for candidate in case["objective"]["singleton_flow"]["candidate_partitions"]:
                local_ref = {p: ref[p] for p in candidate["prediction"] if p in ref}
                candidate["reference_pairs"] = diagnostic_pairs(local_ref, candidate["prediction"])
        case["raw_conditioned_cut_pairs"] = diagnostic_pairs(ref, cut["raw_conditioned_prediction"])
        case["exact_minus_raw_tp"] = score["tp"] - case["raw_conditioned_cut_pairs"]["tp"]
        case["exact_minus_raw_fp"] = score["fp"] - case["raw_conditioned_cut_pairs"]["fp"]
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
        if block_degeneracy:
            from benchmarks.og_extraction.species_block_degeneracy import annotate_reference

            rs = case["objective"]["block_degeneracy"]["legal_internal_nodes"]
            ns, ms = degeneracy_inputs[case["component_id"], case["source_cluster_id"]]
            annotate_reference(rs, ns, ms, reference)
            if (
                sum(r["direct_removed_tp"] for r in rs if r["dp_state"] == "active_split")
                != case["removed_merge_tp"]
                or sum(r["direct_removed_fp"] for r in rs if r["dp_state"] == "active_split")
                != case["removed_merge_fp"]
            ):
                raise ValueError("LCA attribution differs from complete-cut pair losses")
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
                exact_minus_raw_tp=sum(c["exact_minus_raw_tp"] for c in group),
                exact_minus_raw_fp=sum(c["exact_minus_raw_fp"] for c in group),
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
    degeneracy_strata = []
    degeneracy_source_summary = []
    if block_degeneracy:
        cohorts = defaultdict(list)
        for case in cases:
            for r in case["objective"]["block_degeneracy"]["legal_internal_nodes"]:
                cohorts[
                    case["reference_cohort"], r["local_reference_cohort"], r["dp_state"]
                ].append(r)
        for (cohort, local, state), rs in sorted(cohorts.items()):
            fractions = [
                r["single_side_equal_weight_fraction"]
                for r in rs
                if r["single_side_equal_weight_fraction"] is not None
            ]
            degeneracy_strata.append(
                dict(
                    source_cohort=cohort,
                    local_reference_cohort=local,
                    dp_state=state,
                    context_turns_node_positive=sum(r["context_turns_node_positive"] for r in rs)
                    if context_flow
                    else None,
                    internal_gain_sum=sum(r["internal_only_gain"] or 0 for r in rs)
                    if context_flow
                    else None,
                    external_gain_sum=sum(r["external_gain"] or 0 for r in rs)
                    if context_flow
                    else None,
                    legal_nodes=len(rs),
                    connected_nodes=sum(r["cross_connected_blocks"] > 0 for r in rs),
                    degenerate_nodes=sum(r["single_side_connected_equal_blocks"] > 0 for r in rs),
                    connected_blocks=sum(r["cross_connected_blocks"] for r in rs),
                    degenerate_blocks=sum(r["single_side_connected_equal_blocks"] for r in rs),
                    one_side_opposite_distributed_blocks=sum(
                        r["one_side_opposite_distributed_blocks"] for r in rs
                    ),
                    both_sides_concentrated_blocks=sum(
                        r["both_sides_concentrated_blocks"] for r in rs
                    ),
                    connected_cross_weight=sum(r["connected_cross_weight"] for r in rs),
                    degenerate_cross_weight=sum(r["single_side_equal_weight"] for r in rs),
                    weight_fraction_q10_q50_q90=np.quantile(fractions, [0.1, 0.5, 0.9]).tolist()
                    if fractions
                    else None,
                    lca_direct_tp=sum(r["direct_removed_tp"] for r in rs),
                    lca_direct_fp=sum(r["direct_removed_fp"] for r in rs),
                )
            )
    if block_degeneracy:
        source_cohorts = defaultdict(list)
        for case in cases:
            rs = case["objective"]["block_degeneracy"]["legal_internal_nodes"]
            connected = sum(r["cross_connected_blocks"] > 0 for r in rs)
            degenerate = sum(r["single_side_connected_equal_blocks"] > 0 for r in rs)
            source_cohorts[case["reference_cohort"]].append(
                dict(
                    legal_nodes=len(rs),
                    connected_nodes=connected,
                    degenerate_nodes=degenerate,
                    connected_node_fraction=degenerate / connected if connected else None,
                )
            )
        for cohort, rs in sorted(source_cohorts.items()):
            values = [
                r["connected_node_fraction"] for r in rs if r["connected_node_fraction"] is not None
            ]
            degeneracy_source_summary.append(
                dict(
                    source_cohort=cohort,
                    sources=len(rs),
                    connected_sources=sum(r["connected_nodes"] > 0 for r in rs),
                    degenerate_sources=sum(r["degenerate_nodes"] > 0 for r in rs),
                    per_source_node_fraction_q10_q50_q90=np.quantile(
                        values, [0.1, 0.5, 0.9]
                    ).tolist()
                    if values
                    else None,
                )
            )
    context_strata = []
    if context_flow:
        context_cohorts = defaultdict(list)
        for case in cases:
            for row in case["objective"]["block_degeneracy"]["legal_internal_nodes"]:
                for group in row["context_groups"]:
                    context_cohorts[
                        case["reference_cohort"],
                        row["local_reference_cohort"],
                        row["dp_state"],
                        group["kind"],
                    ].append(group)
        for (cohort, local, state, kind), groups in sorted(context_cohorts.items()):
            fields = (
                "blocks",
                "internal_gain",
                "external_gain",
                "external_present_blocks",
                "source_context_changes_sign",
                "internal_nonpositive_total_positive",
                "observed",
                "expected",
                "expected_internal_internal",
                "expected_internal_external",
                "expected_external_external",
                "node_internal_weight",
                "node_boundary_weight",
                "outside_node_weight",
            )
            context_strata.append(
                dict(
                    source_cohort=cohort,
                    local_reference_cohort=local,
                    dp_state=state,
                    kind=kind,
                    legal_nodes=len(groups),
                    nodes_with_blocks=sum(g["blocks"] > 0 for g in groups),
                    **{key: sum(g[key] for g in groups) for key in fields},
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

    def compact_case(case):
        if not block_degeneracy:
            return case
        result = dict(case, objective=dict(case["objective"]))
        profile = case["objective"]["block_degeneracy"]
        rs = profile["legal_internal_nodes"]
        result["objective"]["block_degeneracy"] = dict(
            criterion=profile["criterion"],
            scope=profile["scope"],
            legal_nodes=len(rs),
            degenerate_nodes=sum(r["single_side_connected_equal_blocks"] > 0 for r in rs),
            key_nodes=[r for r in rs if r["cluster_id"] in (21396, 10052)],
        )
        return result

    (out / "summary.json").write_text(
        json.dumps(
            dict(
                diagnostic_only=True,
                production_changed=False,
                immutable_inputs_verified=True,
                reference_labels_used_for_objective=False,
                context_flow=context_flow,
                context_strata=context_strata,
                context_scope=(
                    "Fixed source denominator; "
                    "internal/internal + internal/external + external/external. "
                    "Summed block weights across nested nodes are repeated incidence, "
                    "not unique edge fractions."
                ),
                block_degeneracy=block_degeneracy,
                degeneracy_strata=degeneracy_strata,
                degeneracy_source_summary=degeneracy_source_summary,
                singleton_flow=singleton_flow,
                exact_conditioned=exact_conditioned,
                raw_conditioned_replay_verified=raw_replay_dir is not None,
                tie_audit_replay_verified=tie_replay_dir is not None,
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
                top_polluted=[
                    compact_case(c)
                    for c in sorted(
                        [c for c in cases if c["reference_cohort"] == "polluted"],
                        key=lambda c: -c["merge_transitions"]["new_fp"],
                    )[:10]
                ],
                top_clean=[
                    compact_case(c)
                    for c in sorted(
                        [c for c in cases if c["reference_cohort"] == "clean_tp_gain"],
                        key=lambda c: -c["merge_transitions"]["new_tp"],
                    )[:10]
                ],
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
    p.add_argument("--tie-replay-dir", type=Path)
    p.add_argument("--exact-conditioned", action="store_true")
    p.add_argument("--raw-replay-dir", type=Path)
    p.add_argument("--singleton-flow", action="store_true")
    p.add_argument("--block-degeneracy", action="store_true")
    p.add_argument("--context-flow", action="store_true")
    a = p.parse_args()
    diagnose(
        a.run,
        a.reference_dir,
        a.trace_dir,
        a.out,
        condition_species_pairs=a.condition_species_pairs,
        tie_replay_dir=a.tie_replay_dir,
        exact_conditioned=a.exact_conditioned,
        raw_replay_dir=a.raw_replay_dir,
        singleton_flow=a.singleton_flow,
        block_degeneracy=a.block_degeneracy,
        context_flow=a.context_flow,
    )


if __name__ == "__main__":
    main()
