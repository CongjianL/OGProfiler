"""Read-only endpoint products and exact DP paths for context sign reversals."""

from fractions import Fraction


def endpoint_products(children, ai, bi, ax, bx, denominator, actual):
    """Unordered child pairs; both species orientations, including same-species factors."""
    pairs = []
    totals = [Fraction(), Fraction(), Fraction(), Fraction()]
    for i, left in enumerate(children):
        for j in range(i + 1, len(children)):
            right = children[j]
            ii = (ai[i] * bi[j] + ai[j] * bi[i]) / denominator
            ix = (ai[i] * bx[j] + ax[i] * bi[j] + ai[j] * bx[i] + ax[j] * bi[i]) / denominator
            xx = (ax[i] * bx[j] + ax[j] * bx[i]) / denominator
            obs = actual[min(left, right), max(left, right)]
            for k, value in enumerate((ii, ix, xx, obs)):
                totals[k] += value
            if ii or ix or xx or obs:
                pairs.append(
                    dict(
                        left_child=left,
                        right_child=right,
                        observed=float(obs),
                        expected_ii=float(ii),
                        expected_ix=float(ix),
                        expected_xx=float(xx),
                        deficit=float(ii + ix + xx - obs),
                        deficit_exact=str(ii + ix + xx - obs),
                    )
                )
    within_ix = (
        sum((x * y + u * v for x, y, u, v in zip(ai, bx, ax, bi, strict=True)), Fraction())
        / denominator
    )
    return dict(
        endpoints=[
            dict(
                child_id=c,
                internal_left=float(a),
                internal_right=float(b),
                external_left=float(x),
                external_right=float(y),
            )
            for c, a, b, x, y in zip(children, ai, bi, ax, bx, strict=True)
        ],
        child_pairs=pairs,
        within_child_ix_excluded=float(within_ix),
    ), totals


def attach_paths(rows, by_id, children, order, scores, selected):
    """Reconstruct the existing exact DP and prove the root frontier matches supplied cut."""
    best, keep = {}, {}
    for c in reversed(order):
        split = sum((best[ch] for ch in children[c]), Fraction())
        keep[c] = not children[c] or scores[c] >= split
        best[c] = scores[c] if keep[c] else split

    def frontier(root):
        result, stack = [], [root]
        while stack:
            c = stack.pop()
            if keep[c]:
                result.append(c)
            else:
                stack.extend(reversed(sorted(children[c])))
        return result

    if set(frontier(order[0])) != set(selected):
        raise ValueError("Endpoint diagnostic DP frontier differs from fixed cut")
    for row in rows:
        if not row["context_turns_node_positive"]:
            continue
        c = row["cluster_id"]
        path, parent = [c], by_id[c]["parent_id"]
        while parent is not None:
            path.append(parent)
            parent = by_id[parent]["parent_id"]
        path.reverse()
        selected_ancestor = next((p for p in path[:-1] if p in selected), None)
        local = frontier(c)
        if c in local:
            raise ValueError("Strict positive direct gain must choose local DP split")
        descendant_gain = sum((best[ch] - scores[ch] for ch in children[c]), Fraction())
        direct_gain = sum((scores[ch] for ch in children[c]), Fraction()) - scores[c]
        optimized_gain = best[c] - scores[c]
        if direct_gain + descendant_gain != optimized_gain:
            raise ValueError("Direct/descendant DP margin decomposition mismatch")
        visited, stack = [], [c]
        while stack:
            p = stack.pop()
            visited.append(
                dict(
                    cluster_id=p,
                    parent_id=by_id[p]["parent_id"],
                    local_decision="keep" if keep[p] else "split",
                    actual_selected=p in selected,
                    keep_score=float(scores[p]),
                    best_score=float(best[p]),
                    split_minus_keep_exact=str(
                        sum((best[ch] for ch in children[p]), Fraction()) - scores[p]
                    )
                    if children[p]
                    else None,
                )
            )
            if not keep[p]:
                stack.extend(reversed(sorted(children[p])))
        row["path_trace"] = dict(
            local_dp_traversal=visited,
            root_path=[
                dict(
                    cluster_id=p,
                    selected=p in selected,
                    local_decision="keep" if keep[p] else "split",
                    keep_score=float(scores[p]),
                    best_score=float(best[p]),
                    split_minus_keep_exact=str(
                        sum((best[ch] for ch in children[p]), Fraction()) - scores[p]
                    )
                    if children[p]
                    else None,
                )
                for p in path
            ],
            selected_ancestor=selected_ancestor,
            local_frontier=local,
            actual_frontier=local if selected_ancestor is None else [selected_ancestor],
            child_frontiers=[dict(child_id=ch, local_frontier=frontier(ch)) for ch in children[c]],
            direct_gain_exact=str(direct_gain),
            descendant_optimization_gain_exact=str(descendant_gain),
            optimized_split_minus_keep_exact=str(optimized_gain),
            dp_replay_verified=True,
        )


def annotate_path_losses(rows, nodes):
    """Unique LCA losses for traced subset; subtree sums separate from direct attribution."""
    from benchmarks.og_extraction.tree_cut import topology

    _, children, order = topology(nodes)
    by_row = {r["cluster_id"]: r for r in rows}
    actual_losses = {}
    for c in reversed(order):
        r = by_row.get(c)
        own = (
            (r["direct_removed_tp"], r["direct_removed_fp"])
            if (r and r["dp_state"] == "active_split")
            else (0, 0)
        )
        actual_losses[c] = tuple(
            own[k] + sum(actual_losses[ch][k] for ch in children[c]) for k in range(2)
        )
    for r in rows:
        if "path_trace" not in r:
            continue
        c = r["cluster_id"]
        trace = r["path_trace"]
        trace["actual_direct_lca_tp"] = (
            r["direct_removed_tp"] if r["dp_state"] == "active_split" else 0
        )
        trace["actual_direct_lca_fp"] = (
            r["direct_removed_fp"] if r["dp_state"] == "active_split" else 0
        )
        trace["actual_subtree_removed_tp"], trace["actual_subtree_removed_fp"] = actual_losses[c]
        trace["actual_lca_nodes"] = [
            dict(
                cluster_id=v["cluster_id"],
                direct_removed_tp=by_row[v["cluster_id"]]["direct_removed_tp"],
                direct_removed_fp=by_row[v["cluster_id"]]["direct_removed_fp"],
            )
            for v in trace["local_dp_traversal"]
            if v["cluster_id"] in by_row and by_row[v["cluster_id"]]["dp_state"] == "active_split"
        ]
        if (
            tuple(
                sum(v["direct_removed_" + k] for v in trace["actual_lca_nodes"])
                for k in ("tp", "fp")
            )
            != actual_losses[c]
        ):
            raise ValueError("Actual DP traversal LCA losses disagree with subtree")
        trace["subtree_loss_scope"] = "Nested subtree totals are not additive across traced nodes"
