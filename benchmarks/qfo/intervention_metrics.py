"""Posthoc unique LCA transitions and explicit full/scoped reference-universe scores."""

from collections import Counter

from benchmarks.og_extraction.species_block_degeneracy import annotate_reference


def annotate_intervention(result, nodes, membership, reference):
    from benchmarks.qfo.internal_boundary import diagnostic_pairs

    rows = result["node_transitions"]
    annotate_reference(rows, nodes, membership, reference)
    counts = Counter(
        dict.fromkeys(("recovered_tp", "new_merge_fp", "new_split_tp", "removed_fp"), 0)
    )
    for r in rows:
        a = r["baseline_state"] == "active_split"
        b = r["candidate_state"] == "active_split"
        if a and not b:
            counts["recovered_tp"] += r["direct_removed_tp"]
            counts["new_merge_fp"] += r["direct_removed_fp"]
        if b and not a:
            counts["new_split_tp"] += r["direct_removed_tp"]
            counts["removed_fp"] += r["direct_removed_fp"]
    ref = {p: reference[p] for p in result["baseline_prediction"] if p in reference}
    baseline = diagnostic_pairs(ref, result["baseline_prediction"])
    candidate = diagnostic_pairs(ref, result["candidate_prediction"])
    if candidate["tp"] - baseline["tp"] != counts["recovered_tp"] - counts["new_split_tp"]:
        raise ValueError("Unique LCA TP transitions differ from complete partition")
    if candidate["fp"] - baseline["fp"] != counts["new_merge_fp"] - counts["removed_fp"]:
        raise ValueError("Unique LCA FP transitions differ from complete partition")
    if counts["new_split_tp"] or counts["removed_fp"]:
        raise ValueError("Same-null ablation violated proven coarsening pair direction")
    result["pair_transitions"] = dict(counts)
    result["baseline_pairs"] = baseline
    result["candidate_pairs"] = candidate


def metrics(reference, prediction):
    from benchmarks.qfo.internal_boundary import diagnostic_pairs

    if set(reference) - set(prediction):
        raise ValueError("Reference universe missing from prediction")
    score = diagnostic_pairs(reference, prediction)
    rp = sum(n * (n - 1) // 2 for n in Counter(reference.values()).values())
    pp = score["tp"] + score["fp"]
    return dict(
        assigned_proteins=len(reference),
        reference_pairs=rp,
        predicted_pairs=pp,
        tp=score["tp"],
        fp=score["fp"],
        fn=rp - score["tp"],
        precision=score["tp"] / pp if pp else 0.0,
        recall=score["tp"] / rp if rp else 0.0,
        pair_f1=2 * score["tp"] / (pp + rp) if pp + rp else 0.0,
    )


def pair_transitions(reference, baseline, candidate):
    shared = Counter((reference[p], baseline[p], candidate[p]) for p in reference)
    shared_tp = sum(n * (n - 1) // 2 for n in shared.values())
    joint = Counter((baseline[p], candidate[p]) for p in reference)
    shared_pp = sum(n * (n - 1) // 2 for n in joint.values())
    a = metrics(reference, baseline)
    b = metrics(reference, candidate)
    return dict(
        recovered_tp=b["tp"] - shared_tp,
        new_split_tp=a["tp"] - shared_tp,
        new_merge_fp=b["fp"] - (shared_pp - shared_tp),
        removed_fp=a["fp"] - (shared_pp - shared_tp),
    )


def summarize(cases, reference, actual):
    baseline = {p: ("outside", g) for p, g in actual.items()}
    candidate = dict(baseline)
    targets = set()
    direct = Counter()
    changed = Counter()
    cohorts = {}
    for case in cases:
        result = case["objective"]["same_attraction_intervention"]
        identity = (case["component_id"], case["source_cluster_id"])
        if targets.intersection(result["baseline_prediction"]):
            raise ValueError("Overlapping intervention sources")
        targets.update(result["baseline_prediction"])
        for label, pred in (
            (baseline, result["baseline_prediction"]),
            (candidate, result["candidate_prediction"]),
        ):
            label.update({p: ("source", *identity, g) for p, g in pred.items()})
        direct.update(result["pair_transitions"])
        cohort = case["reference_cohort"]
        acc = cohorts.setdefault(cohort, Counter())
        acc.update(result["pair_transitions"])
        acc["sources"] += 1
        acc["changed_sources"] += result["baseline_prediction"] != result["candidate_prediction"]
        for r in result["node_transitions"]:
            changed["local_decision_changes"] += r["local_decision_changed"]
            changed["active_state_changes"] += r["active_state_changed"]
            changed["ancestor_only_state_changes"] += r["ancestor_only_state_change"]
    if targets - set(actual):
        raise ValueError("Intervention target absent from frozen production")
    ref = {p: g for p, g in reference.items() if p in actual}
    scoped = {p: g for p, g in ref.items() if p in targets}
    global_transitions = pair_transitions(ref, baseline, candidate)
    for k in ("recovered_tp", "new_split_tp", "new_merge_fp", "removed_fp"):
        if global_transitions[k] != direct[k]:
            raise ValueError("Full hybrid transitions differ from unique source LCA attribution")
    return dict(
        algorithm="same_attraction_cross_null_v1",
        production_changed=False,
        full_reference_universe=dict(
            baseline=metrics(ref, baseline),
            candidate=metrics(ref, candidate),
            frozen_production=metrics(ref, actual),
        ),
        target_reference_universe=dict(
            baseline=metrics(scoped, baseline), candidate=metrics(scoped, candidate)
        ),
        pair_transitions=global_transitions,
        node_changes=dict(changed),
        cohorts=[dict(cohort=k, **v) for k, v in sorted(cohorts.items())],
        total_proteins=len(actual),
        target_proteins=len(targets),
        scope="Full hybrid: 1205 diagnostic source clades replaced by their complete cuts; "
        "outside proteins retain frozen production groups; source/outside namespaces are separated "
        "identically in both arms. Frozen production is a separate comparator. "
        "Target-universe recall includes cross-source reference pairs as unrecovered. "
        "Not a new end-to-end production run or official QFO ground truth.",
    )
