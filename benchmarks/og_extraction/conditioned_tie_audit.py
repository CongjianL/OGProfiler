"""Exact-rational local keep/split audit; diagnostic only, no score changes."""

from collections import Counter
from fractions import Fraction

from benchmarks.og_extraction.species_pair_graph import species_pair_cut
from benchmarks.og_extraction.tree_cut import topology


def audit_tree(nodes, membership, edges, species):
    raw = species_pair_cut(nodes, membership, edges, species)
    by_id, children, order = topology(nodes)
    proteins = {c: set() for c in order}
    for p, leaf in membership:
        proteins[leaf].add(p)
    for c in reversed(order):
        for child in children[c]:
            proteins[c].update(proteins[child])
    # Exact stored binary64 weights, not rounded decimal strings or fitted labels.
    rational_edges = [(u, v, Fraction.from_float(float(w))) for u, v, w in edges]
    weights = Counter()
    for u, v, w in rational_edges:
        s, t = species[u], species[v]
        weights[min(s, t), max(s, t)] += w
    total = sum(weights.values(), Fraction())
    scores = {}
    for c in order:
        ps = proteins[c]
        k = Counter()
        internal = Fraction()
        for u, v, w in rational_edges:
            s, t = species[u], species[v]
            if u in ps:
                k[s, t] += w
            if v in ps:
                k[t, s] += w
            if u in ps and v in ps:
                internal += w
        expected = sum(
            (k[s, s] ** 2 / (4 * w) if s == t else k[s, t] * k[t, s] / w)
            for (s, t), w in weights.items()
            if w
        )
        scores[c] = (internal - expected) / total if total else Fraction()
    assert scores[order[0]] == 0
    raw_scores = {r["cluster_id"]: r["keep_score"] for r in raw["node_scores"]}
    exact_best, raw_best, exact_keep, raw_keep, rows = {}, {}, {}, {}, []
    for c in reversed(order):
        exact_split = (
            sum((exact_best[ch] for ch in children[c]), Fraction()) if children[c] else None
        )
        raw_split = sum(raw_best[ch] for ch in children[c]) if children[c] else None
        exact_keep[c] = exact_split is None or scores[c] >= exact_split
        raw_keep[c] = raw_split is None or raw_scores[c] >= raw_split
        exact_best[c] = scores[c] if exact_keep[c] else exact_split
        raw_best[c] = raw_scores[c] if raw_keep[c] else raw_split
        gap = exact_split - scores[c] if exact_split is not None else None
        rows.append(
            dict(
                cluster_id=c,
                parent_id=by_id[c]["parent_id"],
                proteins=sorted(proteins[c]),
                raw_keep_score=raw_scores[c],
                raw_split_best=raw_split,
                raw_split_minus_keep=raw_split - raw_scores[c] if raw_split is not None else None,
                exact_keep_score=float(scores[c]),
                exact_split_best=float(exact_split) if exact_split is not None else None,
                exact_split_minus_keep=float(gap) if gap is not None else None,
                exact_gap_sign=(1 if gap > 0 else -1 if gap < 0 else 0)
                if gap is not None
                else None,
                exact_tie=gap == 0 if gap is not None else False,
                raw_keep=raw_keep[c],
                exact_keep=exact_keep[c],
                choice_changed=raw_keep[c] != exact_keep[c],
            )
        )

    def partition(keep):
        owner, selected = {}, []
        for c in order:
            parent = by_id[c]["parent_id"]
            if parent is not None and owner[parent] is not None:
                owner[c] = owner[parent]
            elif keep[c]:
                owner[c] = c
                selected.append(c)
            else:
                owner[c] = None
        return selected, {p: owner[leaf] for p, leaf in membership}

    def child_partition(c, keep):
        chosen, stack = [], list(children[c])
        while stack:
            child = stack.pop()
            if keep[child]:
                chosen.append(child)
            else:
                stack.extend(children[child])
        return {p: ch for ch in chosen for p in proteins[ch]}

    for row in rows:
        c = row["cluster_id"]
        row["raw_child_prediction"] = child_partition(c, raw_keep)
        row["exact_child_prediction"] = child_partition(c, exact_keep)

    def boundary(child_pred):
        k, actual = Counter(), Counter()
        for u, v, w in rational_edges:
            s, t = species[u], species[v]
            if u in child_pred:
                k[child_pred[u], s, t] += w
            if v in child_pred:
                k[child_pred[v], t, s] += w
            if u in child_pred and v in child_pred and child_pred[u] != child_pred[v]:
                actual["same" if s == t else "cross"] += w
        expected = Counter()
        gs = sorted(set(child_pred.values()))
        for (s, t), w in weights.items():
            if not w:
                continue
            if s == t:
                vals = [k[g, s, s] for g in gs]
                expected["same"] += (sum(vals) ** 2 - sum(x * x for x in vals)) / (4 * w)
            else:
                a, b = [k[g, s, t] for g in gs], [k[g, t, s] for g in gs]
                expected["cross"] += (
                    sum(a) * sum(b) - sum(x * y for x, y in zip(a, b, strict=True))
                ) / w
        return {
            kind: dict(
                observed=float(actual[kind]),
                expected=float(expected[kind]),
                gain=float((expected[kind] - actual[kind]) / total) if total else 0.0,
            )
            for kind in ("same", "cross")
        }

    for row in rows:
        row["raw_child_boundary"] = boundary(row["raw_child_prediction"])
        row["exact_child_boundary"] = boundary(row["exact_child_prediction"])
    selected, pred = partition(exact_keep)
    raw_selected, raw_pred = partition(raw_keep)
    for row in rows:
        ancestors, parent = [], row["parent_id"]
        while parent is not None:
            ancestors.append(parent)
            parent = by_id[parent]["parent_id"]
        row["raw_active"] = not any(c in raw_selected for c in ancestors)
        row["exact_active"] = not any(c in selected for c in ancestors)
    if raw_selected != raw["selected"]:
        raise ValueError("Floating-point selection trace mismatch")
    # Independent exact cross-group deficit identity.
    endpoint = Counter()
    observed = Fraction()
    for u, v, w in rational_edges:
        endpoint[pred[u], species[u], species[v]] += w
        endpoint[pred[v], species[v], species[u]] += w
        if pred[u] != pred[v]:
            observed += w
    expected_cut = Fraction()
    for (s, t), w in weights.items():
        if not w:
            continue
        if s == t:
            vals = [endpoint[g, s, s] for g in selected]
            expected_cut += (sum(vals) ** 2 - sum(x * x for x in vals)) / (4 * w)
        else:
            a = [endpoint[g, s, t] for g in selected]
            b = [endpoint[g, t, s] for g in selected]
            expected_cut += (sum(a) * sum(b) - sum(x * y for x, y in zip(a, b, strict=True))) / w
    independent = (expected_cut - observed) / total if total else Fraction()
    if independent != exact_best[order[0]]:
        raise ValueError("Exact objective identity mismatch")
    raw_exact_score = sum((scores[c] for c in raw_selected), Fraction())
    if independent < raw_exact_score:
        raise ValueError("Exact DP worse than replayed raw cut")
    return dict(
        raw_cut_exact_score=float(raw_exact_score),
        exact_gain_over_raw_cut=float(independent - raw_exact_score),
        raw_cut_exact_score_equal=independent == raw_exact_score,
        exact_selected=selected,
        exact_prediction=pred,
        raw_prediction=raw_pred,
        exact_score=float(independent),
        exact_objective_identity_verified=True,
        arithmetic="Fraction.from_float(binary64 edge weights); exact equality, no tolerance",
        local_rows=sorted(rows, key=lambda r: r["cluster_id"]),
    )
