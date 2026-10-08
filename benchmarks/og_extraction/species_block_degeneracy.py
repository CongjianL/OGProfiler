"""Fixed legal immediate-child cuts: exact one-sided block degeneracy, no filter."""

from collections import Counter
from fractions import Fraction

from benchmarks.og_extraction.exact_species_pair_graph import exact_context
from benchmarks.og_extraction.tree_cut import topology


def profile(nodes, membership, edges, species, selected):
    scores, total, strengths, weights = exact_context(nodes, membership, edges, species)
    by_id, children, order = topology(nodes)
    leaves = dict(membership)
    depths = {}
    for c in order:
        parent = by_id[c]["parent_id"]
        depths[c] = depths[parent] + 1 if parent is not None else 0
    observed = Counter()
    for u, v, w in edges:
        a, b = leaves[u], leaves[v]
        while depths[a] > depths[b]:
            a = by_id[a]["parent_id"]
        while depths[b] > depths[a]:
            b = by_id[b]["parent_id"]
        while a != b:
            a, b = by_id[a]["parent_id"], by_id[b]["parent_id"]
        if children[a]:
            s, t = sorted((species[u], species[v]))
            observed[a, s, t] += Fraction.from_float(float(w))
    rows = []
    for c in order:
        cs = children[c]
        if not cs:
            continue
        blocks = []
        delta = Fraction()
        for (s, t), w in sorted(weights.items()):
            if not w:
                continue
            a, b = [strengths[ch][s, t] for ch in cs], [strengths[ch][t, s] for ch in cs]
            expected = (
                (sum(a) ** 2 - sum(x * x for x in a)) / (4 * w)
                if s == t
                else (sum(a) * sum(b) - sum(x * y for x, y in zip(a, b, strict=True))) / w
            )
            actual = observed[c, s, t]
            if not expected and not actual:
                continue
            left = [ch for ch, x in zip(cs, a, strict=True) if x == w] if s != t else []
            right = [ch for ch, x in zip(cs, b, strict=True) if x == w] if s != t else []
            carrier = bool(left or right)
            equal = actual > 0 and expected == actual
            if carrier and actual > 0 and not equal:
                raise ValueError("One-sided concentration algebra identity failed")
            delta += expected - actual
            blocks.append(
                dict(
                    species_a=s,
                    species_b=t,
                    same_species=s == t,
                    source_block_weight=float(w),
                    observed=float(actual),
                    expected=float(expected),
                    nonzero_deficit=expected != actual,
                    connected=actual > 0,
                    connected_equal=equal,
                    one_side_carriers_left=left,
                    one_side_carriers_right=right,
                    occupied_children_left=sum(x > 0 for x in a),
                    occupied_children_right=sum(x > 0 for x in b),
                    single_side_connected_equal=carrier and equal,
                    both_sides_concentrated=bool(left and right) and equal,
                )
            )
        direct_gain = sum((scores[ch] for ch in cs), Fraction()) - scores[c]
        if direct_gain != (delta / total if total else Fraction()):
            raise ValueError("Immediate-child objective and block decomposition differ")
        ancestors, parent = [], by_id[c]["parent_id"]
        while parent is not None:
            ancestors.append(parent)
            parent = by_id[parent]["parent_id"]
        cross = [r for r in blocks if not r["same_species"]]
        connected = [r for r in cross if r["connected"]]
        degenerate = [r for r in cross if r["single_side_connected_equal"]]
        weight = sum(r["observed"] for r in connected)
        deg_weight = sum(r["observed"] for r in degenerate)
        examples = (
            blocks
            if c == 21396
            else (
                [r for r in blocks if r["single_side_connected_equal"]][:2]
                + [r for r in blocks if r["nonzero_deficit"]][:2]
            )
        )
        rows.append(
            dict(
                cluster_id=c,
                child_count=len(cs),
                children=cs,
                dp_state="inactive"
                if any(p in selected for p in ancestors)
                else "active_keep"
                if c in selected
                else "active_split",
                direct_gain=float(direct_gain),
                direct_gain_sign=1 if direct_gain > 0 else -1 if direct_gain < 0 else 0,
                cross_active_blocks=len(cross),
                cross_connected_blocks=len(connected),
                cross_connected_equal_blocks=sum(r["connected_equal"] for r in cross),
                single_side_connected_equal_blocks=len(degenerate),
                one_side_opposite_distributed_blocks=sum(
                    bool(
                        r["one_side_carriers_left"]
                        and r["occupied_children_right"] >= 2
                        or r["one_side_carriers_right"]
                        and r["occupied_children_left"] >= 2
                    )
                    for r in degenerate
                ),
                both_sides_concentrated_blocks=sum(r["both_sides_concentrated"] for r in cross),
                decisive_blocks=sum(r["nonzero_deficit"] for r in blocks),
                connected_cross_weight=weight,
                single_side_equal_weight=deg_weight,
                single_side_equal_weight_fraction=deg_weight / weight if weight else None,
                objective_identity_verified=True,
                blocks=examples,
                blocks_truncated=len(examples) < len(blocks),
                block_examples_only=c != 21396,
            )
        )
    return dict(
        legal_internal_nodes=rows,
        objective_identity_verified=True,
        criterion=(
            "cross block: actual>0, exact expected=actual; "
            "one source endpoint side entirely in one direct child"
        ),
        scope=(
            "immediate-child legal cuts, not child-optimal DP partitions; "
            "no thresholds or selection changes"
        ),
    )


def annotate_reference(rows, nodes, membership, reference):
    _, children, order = topology(nodes)
    counts = {c: Counter() for c in order}
    for p, leaf in membership:
        if p in reference:
            counts[leaf][reference[p]] += 1
    for c in reversed(order):
        for ch in children[c]:
            counts[c].update(counts[ch])

    def pairs(hist):
        n = sum(hist.values())
        tp = sum(v * (v - 1) // 2 for v in hist.values())
        return tp, n * (n - 1) // 2 - tp

    for r in rows:
        c = r["cluster_id"]
        tp, fp = pairs(counts[c])
        child_counts = [pairs(counts[ch]) for ch in children[c]]
        r.update(
            local_reference_cohort="pure"
            if len(counts[c]) == 1
            else "mixed"
            if counts[c]
            else "unassigned",
            assigned_proteins=sum(counts[c].values()),
            reference_families=len(counts[c]),
            direct_removed_tp=tp - sum(t for t, f in child_counts),
            direct_removed_fp=fp - sum(f for t, f in child_counts),
        )
