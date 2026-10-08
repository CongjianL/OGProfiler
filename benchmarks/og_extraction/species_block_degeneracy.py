"""Fixed legal immediate-child cuts: exact one-sided block degeneracy, no filter."""

from collections import Counter
from fractions import Fraction

from benchmarks.og_extraction.exact_species_pair_graph import exact_context
from benchmarks.og_extraction.tree_cut import topology


def profile(
    nodes, membership, edges, species, selected, *, context_flow=False, endpoint_trace=False
):
    if endpoint_trace and not context_flow:
        raise ValueError("Endpoint tracing requires context flow")
    scores, total, strengths, weights = exact_context(nodes, membership, edges, species)
    by_id, children, order = topology(nodes)
    leaves = dict(membership)
    depths = {}
    for c in order:
        parent = by_id[c]["parent_id"]
        depths[c] = depths[parent] + 1 if parent is not None else 0
    observed = Counter()
    at_lca_edges = {c: [] for c in order}
    internal_blocks = {c: Counter() for c in order}
    for u, v, w in edges:
        a, b = leaves[u], leaves[v]
        while depths[a] > depths[b]:
            a = by_id[a]["parent_id"]
        while depths[b] > depths[a]:
            b = by_id[b]["parent_id"]
        while a != b:
            a, b = by_id[a]["parent_id"], by_id[b]["parent_id"]
        if context_flow:
            at_lca_edges[a].append((u, v, Fraction.from_float(float(w))))
            internal_blocks[a][min(species[u], species[v]), max(species[u], species[v])] += (
                Fraction.from_float(float(w))
            )
        if children[a]:
            s, t = sorted((species[u], species[v]))
            observed[a, s, t] += Fraction.from_float(float(w))
    if context_flow:
        for c in reversed(order):
            for ch in children[c]:
                internal_blocks[c].update(internal_blocks[ch])
    rows = []
    for c in order:
        cs = children[c]
        if not cs:
            continue
        blocks = []
        pair_actual = {}
        inner = {ch: Counter() for ch in cs}
        if context_flow:
            # Within-child edges contribute both endpoints; between-child edges one each.
            for ch in cs:
                for (s, t), iw in internal_blocks[ch].items():
                    inner[ch][s, t] += 2 * iw if s == t else iw
                    if s != t:
                        inner[ch][t, s] += iw

            def child_of(p, ancestor=c):
                x = leaves[p]
                while by_id[x]["parent_id"] != ancestor:
                    x = by_id[x]["parent_id"]
                return x

            for u, v, ew in at_lca_edges[c]:
                if endpoint_trace:
                    key = tuple(sorted((species[u], species[v])))
                    cs_pair = tuple(sorted((child_of(u), child_of(v))))
                    pair_actual.setdefault(key, Counter())[cs_pair] += ew
                inner[child_of(u)][species[u], species[v]] += ew
                inner[child_of(v)][species[v], species[u]] += ew
        delta = Fraction()
        delta_internal = Fraction()
        delta_external = Fraction()
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
            context = {}
            if context_flow:
                ai, bi = [inner[ch][s, t] for ch in cs], [inner[ch][t, s] for ch in cs]
                ax, bx = (
                    [x - y for x, y in zip(a, ai, strict=True)],
                    [x - y for x, y in zip(b, bi, strict=True)],
                )
                if any(x < 0 for x in (*ax, *bx)):
                    raise ValueError("Negative external endpoint strength")

                def pair(x, y):
                    return sum(x) * sum(y) - sum(u * v for u, v in zip(x, y, strict=True))

                denom = 4 * w if s == t else w
                eii = pair(ai, bi) / denom
                eix = (pair(ai, bx) + pair(ax, bi)) / denom
                exx = pair(ax, bx) / denom
                delta_internal += eii - actual
                delta_external += eix + exx
                if eii + eix + exx != expected:
                    raise ValueError("Context expectation identity mismatch")
                iw = internal_blocks[c][s, t]
                boundary = (sum(a) + sum(b)) / 2 - 2 * iw if s == t else sum(a) + sum(b) - 2 * iw
                outside = w - iw - boundary
                if min(boundary, outside) < 0:
                    raise ValueError("Invalid internal/boundary/outside weights")
                context = dict(
                    expected_internal_internal=float(eii),
                    expected_internal_external=float(eix),
                    expected_external_external=float(exx),
                    internal_only_deficit=float(eii - actual),
                    external_expected=float(eix + exx),
                    internal_only_sign=1 if eii > actual else -1 if eii < actual else 0,
                    total_deficit_sign=1 if expected > actual else -1 if expected < actual else 0,
                    source_context_changes_sign=((eii > actual) - (eii < actual))
                    != ((expected > actual) - (expected < actual)),
                    node_internal_weight=float(iw),
                    node_boundary_weight=float(boundary),
                    outside_node_weight=float(outside),
                    external_left_strength=float(sum(ax)),
                    external_right_strength=float(sum(bx)),
                    context_identity_verified=True,
                )
            if endpoint_trace:
                from benchmarks.og_extraction.context_path_trace import endpoint_products

                trace, totals = endpoint_products(
                    cs, ai, bi, ax, bx, denom, pair_actual.get((s, t), Counter())
                )
                if totals != [eii, eix, exx, actual]:
                    raise ValueError("Child-pair endpoint decomposition mismatch")
                context["endpoint_trace"] = trace
            delta += expected - actual
            blocks.append(
                dict(
                    **context,
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
        context_groups = []
        if context_flow:
            for kind in ("connected_equal", "decisive_positive", "decisive_negative"):
                rs = [
                    r
                    for r in blocks
                    if (
                        r["connected_equal"]
                        if kind == "connected_equal"
                        else r["total_deficit_sign"] == (1 if kind == "decisive_positive" else -1)
                    )
                ]
                context_groups.append(
                    dict(
                        kind=kind,
                        internal_gain=sum(r["internal_only_deficit"] for r in rs) / float(total)
                        if total
                        else 0.0,
                        external_gain=sum(r["external_expected"] for r in rs) / float(total)
                        if total
                        else 0.0,
                        blocks=len(rs),
                        external_present_blocks=sum(
                            r["external_left_strength"] > 0 or r["external_right_strength"] > 0
                            for r in rs
                        ),
                        source_context_changes_sign=sum(
                            r["source_context_changes_sign"] for r in rs
                        ),
                        internal_nonpositive_total_positive=sum(
                            r["internal_only_sign"] <= 0 and r["total_deficit_sign"] > 0 for r in rs
                        ),
                        node_internal_weight=sum(r["node_internal_weight"] for r in rs),
                        node_boundary_weight=sum(r["node_boundary_weight"] for r in rs),
                        outside_node_weight=sum(r["outside_node_weight"] for r in rs),
                        observed=sum(r["observed"] for r in rs),
                        expected=sum(r["expected"] for r in rs),
                        expected_internal_internal=sum(r["expected_internal_internal"] for r in rs),
                        expected_internal_external=sum(r["expected_internal_external"] for r in rs),
                        expected_external_external=sum(r["expected_external_external"] for r in rs),
                    )
                )
        if context_flow and delta_internal + delta_external != delta:
            raise ValueError("Complete node context identity mismatch")
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
            if c == 21396 or (endpoint_trace and delta_internal <= 0 and delta > 0)
            else (
                [r for r in blocks if r["single_side_connected_equal"]][:2]
                + [r for r in blocks if r["nonzero_deficit"]][:2]
            )
        )
        if endpoint_trace and not (delta_internal <= 0 and delta > 0 or c == 21396):
            for block in blocks:
                block.pop("endpoint_trace", None)
        rows.append(
            dict(
                context_groups=context_groups,
                context_identity_verified=context_flow,
                internal_only_gain=float(delta_internal / total)
                if context_flow and total
                else None,
                external_gain=float(delta_external / total) if context_flow and total else None,
                internal_only_gain_sign=(
                    1 if delta_internal > 0 else -1 if delta_internal < 0 else 0
                )
                if context_flow
                else None,
                context_turns_node_positive=delta_internal <= 0 and delta > 0
                if context_flow
                else None,
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
                block_examples_only=len(examples) < len(blocks),
            )
        )
    if endpoint_trace:
        from benchmarks.og_extraction.context_path_trace import attach_paths

        attach_paths(rows, by_id, children, order, scores, selected)
        composition = {c: Counter() for c in order}
        for protein, leaf in membership:
            composition[leaf][species[protein]] += 1
        for c in reversed(order):
            for child in children[c]:
                composition[c].update(composition[child])
        for row in rows:
            if "path_trace" in row:
                row["path_trace"]["child_composition"] = [
                    dict(
                        child_id=ch,
                        proteins=sum(composition[ch].values()),
                        species_counts=dict(sorted(composition[ch].items())),
                    )
                    for ch in children[row["cluster_id"]]
                ]
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
