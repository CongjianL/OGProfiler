"""Exact stored-edge arithmetic for experimental conditioned complete cuts."""

from collections import Counter
from fractions import Fraction

from benchmarks.og_extraction.fixed_tree_graph import graph_cut
from benchmarks.og_extraction.tree_cut import optimal_cut, topology


def exact_scores(nodes, membership, edges, species):
    by_id, children, order = topology(nodes)
    leaves = dict(membership)
    strengths = {c: Counter() for c in order}
    internal = Counter()
    weights = Counter()
    depth = {}
    for c in order:
        parent = by_id[c]["parent_id"]
        depth[c] = depth[parent] + 1 if parent is not None else 0
    for u, v, w in edges:
        w = Fraction.from_float(float(w))
        s, t = species[u], species[v]
        weights[min(s, t), max(s, t)] += w
        strengths[leaves[u]][s, t] += w
        strengths[leaves[v]][t, s] += w
        a, b = leaves[u], leaves[v]
        while depth[a] > depth[b]:
            a = by_id[a]["parent_id"]
        while depth[b] > depth[a]:
            b = by_id[b]["parent_id"]
        while a != b:
            a, b = by_id[a]["parent_id"], by_id[b]["parent_id"]
        internal[a] += w
    for c in reversed(order):
        for child in children[c]:
            strengths[c].update(strengths[child])
            internal[c] += internal[child]
    total = sum(weights.values(), Fraction())
    scores = {}
    for c in order:
        k = strengths[c]
        expected = sum(
            (k[s, s] ** 2 / (4 * w) if s == t else k[s, t] * k[t, s] / w)
            for (s, t), w in weights.items()
            if w
        )
        scores[c] = (internal[c] - expected) / total if total else Fraction()
    if scores[order[0]] != 0:
        raise ValueError("Exact root score is not zero")
    return scores, total


def exact_cut(nodes, membership, edges, species):
    base = graph_cut(nodes, membership, edges, strength_null=True)
    scores, total = exact_scores(nodes, membership, edges, species)
    cut = optimal_cut(nodes, scores)
    by_id, _, order = topology(nodes)
    owner = {}
    chosen = set(cut.selected)
    for c in order:
        owner[c] = c if c in chosen else owner.get(by_id[c]["parent_id"])
    prediction = {p: owner[leaf] for p, leaf in membership}
    if any(g is None for g in prediction.values()) or len(prediction) != len(membership):
        raise ValueError("Incomplete exact cut")
    # Independently recompute the cross-group null deficit with exact arithmetic.
    weights, k, observed = Counter(), Counter(), Fraction()
    for u, v, w in edges:
        w = Fraction.from_float(float(w))
        s, t = species[u], species[v]
        weights[min(s, t), max(s, t)] += w
        k[prediction[u], s, t] += w
        k[prediction[v], t, s] += w
        if prediction[u] != prediction[v]:
            observed += w
    expected = Fraction()
    for (s, t), w in weights.items():
        if not w:
            continue
        a = [k[g, s, t] for g in cut.selected]
        if s == t:
            expected += (sum(a) ** 2 - sum(x * x for x in a)) / (4 * w)
        else:
            b = [k[g, t, s] for g in cut.selected]
            expected += (sum(a) * sum(b) - sum(x * y for x, y in zip(a, b, strict=True))) / w
    independent = (expected - observed) / total if total else Fraction()
    if independent != cut.score:
        raise ValueError("Exact DP and independent objective differ")
    return dict(
        algorithm="fixed-tree-species-pair-strength-null-exact-v1",
        selected=list(cut.selected),
        score=float(cut.score),
        component_total_weight=float(total),
        members=[dict(protein_id=p, cluster_id=g) for p, g in sorted(prediction.items())],
        node_scores=[
            dict(r, keep_score=float(scores[r["cluster_id"]])) for r in base["node_scores"]
        ],
        reference_labels_used=False,
        tie_policy="exact equality keeps parent",
        objective_identity_verified=True,
        diagnostics=dict(groups=len(cut.selected), proteins=len(prediction)),
    )
