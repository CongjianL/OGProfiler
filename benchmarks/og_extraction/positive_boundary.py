"""Fixed complete boundary objective for internally positive nodes, diagnostic only."""

from collections import Counter, defaultdict
from fractions import Fraction

from benchmarks.og_extraction.tree_cut import topology


def attach_positive(
    rows, nodes, membership, species, scores, total, strengths, weights, selected, raw_internal
):
    by_id, children, order = topology(nodes)
    composition = {c: Counter() for c in order}
    for p, leaf in membership:
        composition[leaf][species[p]] += 1
    depth, expected, internal, best, keep = {}, {}, {}, {}, {}
    for c in order:
        parent = by_id[c]["parent_id"]
        depth[c] = 0 if parent is None else depth[parent] + 1
        k = strengths[c]
        expected[c] = sum(
            (k[s, s] ** 2 / (4 * w) if s == t else k[s, t] * k[t, s] / w)
            for (s, t), w in weights.items()
            if w
        )
        internal[c] = sum(raw_internal[c].values(), Fraction())
        if internal[c] - expected[c] != scores[c] * total:
            raise ValueError("Raw internal flow differs from exact node score")
    frontier_internal, frontier_expected, frontier_count = {}, {}, {}
    for c in reversed(order):
        for ch in children[c]:
            composition[c].update(composition[ch])
        split = sum((best[ch] for ch in children[c]), Fraction())
        keep[c] = not children[c] or scores[c] >= split
        best[c] = scores[c] if keep[c] else split
        frontier_internal[c] = (
            internal[c] if keep[c] else sum(frontier_internal[ch] for ch in children[c])
        )
        frontier_expected[c] = (
            expected[c] if keep[c] else sum(frontier_expected[ch] for ch in children[c])
        )
        frontier_count[c] = 1 if keep[c] else sum(frontier_count[ch] for ch in children[c])
    actual, stack = set(), [order[0]]
    while stack:
        c = stack.pop()
        if keep[c]:
            actual.add(c)
        else:
            stack.extend(children[c])
    if actual != set(selected):
        raise ValueError("Positive-boundary exact DP differs from fixed cut")
    for r in rows:
        c = r["cluster_id"]
        r["complete_dp_keep"] = keep[c]
        if r["internal_only_gain_sign"] != 1:
            continue
        direct = sum((scores[ch] for ch in children[c]), Fraction()) - scores[c]
        descendant = sum((best[ch] - scores[ch] for ch in children[c]), Fraction())
        gain = best[c] - scores[c]
        actual_cross = internal[c] - frontier_internal[c]
        expected_cross = expected[c] - frontier_expected[c]
        if keep[c] or direct <= 0 or direct + descendant != gain:
            raise ValueError("Internally positive node DP decomposition failed")
        if (expected_cross - actual_cross) / total != gain:
            raise ValueError("Complete frontier boundary identity failed")
        sizes = tuple(sorted(sum(composition[ch].values()) for ch in children[c]))
        comp = composition[c]
        r["positive_boundary"] = dict(
            proteins=sum(comp.values()),
            species_count=len(comp),
            max_copy=max(comp.values()),
            child_sizes=sizes,
            depth=depth[c],
            source_size=sum(composition[order[0]].values()),
            merge_score=float(scores[c]),
            direct_cut_score=float(scores[c] + direct),
            optimal_cut_score=float(best[c]),
            direct_gain_exact=str(direct),
            descendant_gain=float(descendant),
            descendant_gain_exact=str(descendant),
            complete_gain=float(gain),
            complete_gain_exact=str(gain),
            frontier_count=frontier_count[c],
            deeper_frontier=frontier_count[c] > len(children[c]),
            complete_cross_actual=float(actual_cross),
            complete_cross_expected=float(expected_cross),
            complete_boundary_identity_verified=True,
        )


def annotate_positive(rows, nodes, membership, reference):
    _, children, order = topology(nodes)
    by_row = {r["cluster_id"]: r for r in rows}
    counts = {c: Counter() for c in order}
    retained, actual_loss = {}, {}
    for p, leaf in membership:
        if p in reference:
            counts[leaf][reference[p]] += 1

    def pairs(hist):
        n = sum(hist.values())
        tp = sum(v * (v - 1) // 2 for v in hist.values())
        return tp, n * (n - 1) // 2 - tp

    for c in reversed(order):
        for ch in children[c]:
            counts[c].update(counts[ch])
        r = by_row.get(c)
        merged = pairs(counts[c])
        retained[c] = (
            merged
            if not r or r["complete_dp_keep"]
            else tuple(sum(retained[ch][i] for ch in children[c]) for i in range(2))
        )
        own = (
            (r["direct_removed_tp"], r["direct_removed_fp"])
            if r and r["dp_state"] == "active_split"
            else (0, 0)
        )
        actual_loss[c] = tuple(
            own[i] + sum(actual_loss[ch][i] for ch in children[c]) for i in range(2)
        )
        if r and "positive_boundary" in r:
            b = r["positive_boundary"]
            local = tuple(merged[i] - retained[c][i] for i in range(2))
            if r["dp_state"] == "active_split" and local != actual_loss[c]:
                raise ValueError("Positive complete-cut losses differ from actual LCA subtree")
            b.update(
                local_complete_removed_tp=local[0],
                local_complete_removed_fp=local[1],
                actual_direct_lca_tp=own[0],
                actual_direct_lca_fp=own[1],
                actual_subtree_removed_tp=actual_loss[c][0],
                actual_subtree_removed_fp=actual_loss[c][1],
            )


def compare_positive(cases):
    """Reference-posthoc strata and exact structural common support; no scoring thresholds."""
    groups, cells = defaultdict(list), defaultdict(lambda: defaultdict(list))
    coverage = defaultdict(list)
    targets = []
    for case in cases:
        for r in case["objective"]["block_degeneracy"]["legal_internal_nodes"]:
            if r["dp_state"] == "active_split":
                category = (
                    "internal_positive"
                    if r["internal_only_gain_sign"] == 1
                    else "context_reversal"
                    if r["context_turns_node_positive"]
                    else "nonpositive_direct_with_descendant_optimum"
                )
                coverage[case["reference_cohort"], r["local_reference_cohort"], category].append(r)
            if "positive_boundary" not in r:
                continue
            b = r["positive_boundary"]
            target = dict(
                component_id=case["component_id"],
                source_cluster_id=case["source_cluster_id"],
                source_cohort=case["reference_cohort"],
                node=r,
            )
            targets.append(target)
            key = case["reference_cohort"], r["local_reference_cohort"], r["dp_state"]
            groups[key].append(r)
            signature = (
                case["reference_cohort"],
                r["dp_state"],
                b["proteins"],
                b["species_count"],
                b["max_copy"],
                tuple(b["child_sizes"]),
                b["depth"],
                b["source_size"],
            )
            cells[signature][r["local_reference_cohort"]].append(r)

    def stats(rs):
        n = len(rs)
        return dict(
            nodes=n,
            deeper_frontier_nodes=sum(r["positive_boundary"]["deeper_frontier"] for r in rs),
            actual_direct_lca_tp=sum(r["positive_boundary"]["actual_direct_lca_tp"] for r in rs),
            actual_direct_lca_fp=sum(r["positive_boundary"]["actual_direct_lca_fp"] for r in rs),
            mean_direct_gain=sum(r["direct_gain"] for r in rs) / n,
            mean_complete_gain=sum(r["positive_boundary"]["complete_gain"] for r in rs) / n,
            mean_descendant_gain=sum(r["positive_boundary"]["descendant_gain"] for r in rs) / n,
        )

    strata = [
        dict(source_cohort=a, local_reference_cohort=b, dp_state=c, **stats(rs))
        for (a, b, c), rs in sorted(groups.items())
    ]
    matched = [
        dict(signature=list(sig), pure=stats(cell["pure"]), mixed=stats(cell["mixed"]))
        for sig, cell in sorted(cells.items())
        if cell["pure"] and cell["mixed"]
    ]
    return dict(
        coverage=[
            dict(
                source_cohort=a,
                local_reference_cohort=b,
                category=c,
                nodes=len(rs),
                actual_direct_lca_tp=sum(r["direct_removed_tp"] for r in rs),
                actual_direct_lca_fp=sum(r["direct_removed_fp"] for r in rs),
            )
            for (a, b, c), rs in sorted(coverage.items())
        ],
        matched_pure_nodes=sum(cell["pure"]["nodes"] for cell in matched),
        matched_mixed_nodes=sum(cell["mixed"]["nodes"] for cell in matched),
        strata=strata,
        matched_cells=matched,
        targets=targets,
        scope="Unique direct LCA losses additive; nested complete/subtree losses not additive. "
        "Exact structural common support is observational, not a causal correction.",
    )
