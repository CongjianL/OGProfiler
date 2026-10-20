"""Same/cross-species x actual/null x direct/descendant frozen DP decomposition."""

from collections import Counter, defaultdict
from fractions import Fraction

from benchmarks.og_extraction.tree_cut import topology


def attach_decomposition(rows, nodes, scores, total, strengths, weights, raw_internal, selected):
    _, children, order = topology(nodes)
    # Two additive channels; same source denominators in keep and every frontier node.
    actual, expected = {}, {}
    for c in order:
        actual[c] = [
            sum((w for (s, t), w in raw_internal[c].items() if (s == t) == same), Fraction())
            for same in (True, False)
        ]
        k = strengths[c]
        expected[c] = [
            sum(
                (
                    k[s, s] ** 2 / (4 * w) if same else k[s, t] * k[t, s] / w
                    for (s, t), w in weights.items()
                    if w and (s == t) == same
                ),
                Fraction(),
            )
            for same in (True, False)
        ]
        if sum(actual[c]) - sum(expected[c]) != scores[c] * total:
            raise ValueError("Stored edge channels differ from exact keep score")
    best, keep, frontier_actual, frontier_expected, count = {}, {}, {}, {}, {}
    for c in reversed(order):
        split = sum((best[ch] for ch in children[c]), Fraction())
        keep[c] = not children[c] or scores[c] >= split
        best[c] = scores[c] if keep[c] else split
        frontier_actual[c] = (
            actual[c]
            if keep[c]
            else [sum(frontier_actual[ch][i] for ch in children[c]) for i in range(2)]
        )
        frontier_expected[c] = (
            expected[c]
            if keep[c]
            else [sum(frontier_expected[ch][i] for ch in children[c]) for i in range(2)]
        )
        count[c] = 1 if keep[c] else sum(count[ch] for ch in children[c])
    chosen, stack = set(), [order[0]]
    while stack:
        c = stack.pop()
        if keep[c]:
            chosen.add(c)
        else:
            stack.extend(children[c])
    if chosen != set(selected):
        raise ValueError("Decomposition changed frozen DP frontier")
    for r in rows:
        c = r["cluster_id"]
        parts = []
        for i, kind in enumerate(("same_species", "cross_species")):
            da = actual[c][i] - sum(actual[ch][i] for ch in children[c])
            de = expected[c][i] - sum(expected[ch][i] for ch in children[c])
            # Descendant optimization sums children's local optimal cuts, even when c keeps.
            xa = sum(actual[ch][i] - frontier_actual[ch][i] for ch in children[c])
            xe = sum(expected[ch][i] - frontier_expected[ch][i] for ch in children[c])
            ca, ce = (0, 0) if keep[c] else (da + xa, de + xe)
            for stage, a, e in (
                ("direct", da, de),
                ("descendant", xa, xe),
                ("forced_split", da + xa, de + xe),
                ("optimal", ca, ce),
            ):
                if a < 0 or e < 0:
                    raise ValueError("Negative boundary weight in exact decomposition")
                gain = (e - a) / total if total else Fraction()
                parts.append(
                    dict(
                        species_channel=kind,
                        stage=stage,
                        actual=float(a),
                        expected=float(e),
                        actual_over_source=float(a / total) if total else 0.0,
                        expected_over_source=float(e / total) if total else 0.0,
                        gain=float(gain),
                        gain_exact=str(gain),
                        gain_sign=(gain > 0) - (gain < 0),
                    )
                )
            if (
                ca != actual[c][i] - frontier_actual[c][i]
                or ce != expected[c][i] - frontier_expected[c][i]
            ):
                raise ValueError("Complete optimal boundary channel identity failed")

        def stage_gain(stage, parts=parts):
            return sum(
                (Fraction(p["gain_exact"]) for p in parts if p["stage"] == stage), Fraction()
            )

        direct = sum((scores[ch] for ch in children[c]), Fraction()) - scores[c]
        forced = sum((best[ch] for ch in children[c]), Fraction()) - scores[c]
        if (
            stage_gain("direct") != direct
            or stage_gain("forced_split") != forced
            or stage_gain("optimal") != best[c] - scores[c]
        ):
            raise ValueError("Exact DP gain/channel identity failed")
        if stage_gain("direct") + stage_gain("descendant") != stage_gain("forced_split"):
            raise ValueError("Direct/descendant identity failed")
        category = (
            "internal_positive"
            if r["internal_only_gain_sign"] == 1
            else "context_reversal"
            if r["context_turns_node_positive"]
            else "nonpositive_direct_with_descendant_optimum"
            if not keep[c]
            else "local_keep"
        )
        if category == "nonpositive_direct_with_descendant_optimum" and (direct > 0 or forced <= 0):
            raise ValueError("Descendant-driven category mismatch")
        r["objective_decomposition"] = dict(
            category=category,
            local_keep=keep[c],
            frontier_count=count[c],
            direct_children=len(children[c]),
            parts=parts,
            source_weight=float(total),
            direct_gain_exact=str(direct),
            forced_split_gain_exact=str(forced),
            optimal_gain_exact=str(best[c] - scores[c]),
            identities_verified=True,
        )
    # Independent per-channel telescoping over actual split LCA boundaries.
    for i in range(2):
        a = sum(
            actual[r["cluster_id"]][i] - sum(actual[ch][i] for ch in children[r["cluster_id"]])
            for r in rows
            if r["dp_state"] == "active_split"
        )
        e = sum(
            expected[r["cluster_id"]][i] - sum(expected[ch][i] for ch in children[r["cluster_id"]])
            for r in rows
            if r["dp_state"] == "active_split"
        )
        root = order[0]
        if (a, e) != (
            actual[root][i] - frontier_actual[root][i],
            expected[root][i] - frontier_expected[root][i],
        ):
            raise ValueError("Unique LCA channel telescoping failed")


def annotate_losses(rows, nodes, membership, reference):
    _, children, order = topology(nodes)
    by_row = {r["cluster_id"]: r for r in rows}
    counts = {c: Counter() for c in order}
    for p, leaf in membership:
        if p in reference:
            counts[leaf][reference[p]] += 1
    retained, losses = {}, {}
    for c in reversed(order):
        for ch in children[c]:
            counts[c].update(counts[ch])
        n = sum(counts[c].values())
        tp = sum(v * (v - 1) // 2 for v in counts[c].values())
        merged = (tp, n * (n - 1) // 2 - tp)
        row = by_row.get(c)
        retained[c] = (
            merged
            if not row or row["objective_decomposition"]["local_keep"]
            else tuple(sum(retained[ch][i] for ch in children[c]) for i in range(2))
        )
        own = (
            (row["direct_removed_tp"], row["direct_removed_fp"])
            if row and row["dp_state"] == "active_split"
            else (0, 0)
        )
        losses[c] = tuple(own[i] + sum(losses[ch][i] for ch in children[c]) for i in range(2))
        if row:
            local = tuple(merged[i] - retained[c][i] for i in range(2))
            if row["dp_state"] == "active_split" and local != losses[c]:
                raise ValueError("Complete pair losses differ from actual LCA subtree")
            row["objective_decomposition"].update(
                actual_direct_lca_tp=own[0],
                actual_direct_lca_fp=own[1],
                local_complete_tp=local[0],
                local_complete_fp=local[1],
                actual_subtree_tp=losses[c][0],
                actual_subtree_fp=losses[c][1],
            )


def summarize(cases):
    strata = defaultdict(list)
    sources = []
    for case in cases:
        rows = case["objective"]["block_degeneracy"]["legal_internal_nodes"]
        for r in rows:
            d = r["objective_decomposition"]
            strata[
                case["reference_cohort"], r["local_reference_cohort"], r["dp_state"], d["category"]
            ].append(r)
        sources.append(
            dict(
                component_id=case["component_id"],
                source_cluster_id=case["source_cluster_id"],
                actual_lca_channels=[
                    dict(
                        species_channel=kind,
                        actual=sum(
                            p["actual"]
                            for r in rows
                            if r["dp_state"] == "active_split"
                            for p in r["objective_decomposition"]["parts"]
                            if p["stage"] == "direct" and p["species_channel"] == kind
                        ),
                        expected=sum(
                            p["expected"]
                            for r in rows
                            if r["dp_state"] == "active_split"
                            for p in r["objective_decomposition"]["parts"]
                            if p["stage"] == "direct" and p["species_channel"] == kind
                        ),
                    )
                    for kind in ("same_species", "cross_species")
                ],
            )
        )
    result = []
    for (source, local, state, category), rs in sorted(strata.items()):
        parts = []
        for stage in ("direct", "descendant", "forced_split", "optimal"):
            for kind in ("same_species", "cross_species"):
                ps = [
                    p
                    for r in rs
                    for p in r["objective_decomposition"]["parts"]
                    if p["stage"] == stage and p["species_channel"] == kind
                ]
                parts.append(
                    dict(
                        stage=stage,
                        species_channel=kind,
                        actual=sum(p["actual"] for p in ps),
                        expected=sum(p["expected"] for p in ps),
                        actual_over_source=sum(p["actual_over_source"] for p in ps),
                        expected_over_source=sum(p["expected_over_source"] for p in ps),
                        gain=sum(p["gain"] for p in ps),
                        positive_nodes=sum(p["gain_sign"] == 1 for p in ps),
                        negative_nodes=sum(p["gain_sign"] == -1 for p in ps),
                        zero_nodes=sum(p["gain_sign"] == 0 for p in ps),
                    )
                )
        result.append(
            dict(
                source_cohort=source,
                local_reference_cohort=local,
                dp_state=state,
                category=category,
                nodes=len(rs),
                actual_direct_lca_tp=sum(
                    r["objective_decomposition"]["actual_direct_lca_tp"] for r in rs
                ),
                actual_direct_lca_fp=sum(
                    r["objective_decomposition"]["actual_direct_lca_fp"] for r in rs
                ),
                parts=parts,
            )
        )
    return dict(
        strata=result,
        sources=sources,
        matching_performed=False,
        objective_changed=False,
        scope=(
            "Only actual active_split direct LCA terms are additive within a source. "
            "Descendant/optimal node totals overlap; gain sums use source-specific normalizers. "
            "No intervention effect estimated."
        ),
    )
