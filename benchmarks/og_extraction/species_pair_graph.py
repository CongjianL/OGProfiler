"""Additive species-pair-conditioned complete fixed-tree cut, diagnostic only."""

import math
from collections import Counter

from benchmarks.og_extraction.fixed_tree_graph import graph_cut
from benchmarks.og_extraction.tree_cut import optimal_cut, topology
from benchmarks.qfo.species_pair_null import conditioned_null


def species_pair_cut(nodes, membership, edges, species, *, exact=True):
    if exact:
        from benchmarks.og_extraction.exact_species_pair_graph import exact_cut

        return exact_cut(nodes, membership, edges, species)
    # Reuse canonical-edge/membership validation and LCA internal weights only.
    base = graph_cut(nodes, membership, edges, strength_null=True)
    by_id, children, order = topology(nodes)
    leaves = dict(membership)
    strengths = {c: Counter() for c in order}
    weights = Counter()
    for u, v, w in edges:
        s, t = species[u], species[v]
        weights[min(s, t), max(s, t)] += w
        strengths[leaves[u]][s, t] += w
        strengths[leaves[v]][t, s] += w
    for c in reversed(order):
        for child in children[c]:
            strengths[c].update(strengths[child])
    total = math.fsum(weights.values())
    node_rows, scores = [], {}
    for row in base["node_scores"]:
        c = row["cluster_id"]
        k = strengths[c]
        expected = math.fsum(
            (k[s, s] ** 2 / (4 * w) if s == t else k[s, t] * k[t, s] / w)
            for (s, t), w in weights.items()
            if w
        )
        scores[c] = (row["internal_weight"] - expected) / total if total else 0.0
        if abs(scores[c]) <= 1e-12:
            scores[c] = 0.0  # Numerical zero, not a tunable biological threshold.
        node_rows.append(dict(row, keep_score=scores[c], expected_internal_weight=expected))
    root = order[0]
    if not math.isclose(scores[root], 0.0, abs_tol=1e-9):
        raise ValueError("Root objective is not zero")
    # Analytically root=0; remove rounding so zero-gain ties retain parent.
    scores[root] = 0.0
    next(r for r in node_rows if r["cluster_id"] == root)["keep_score"] = 0.0
    cut = optimal_cut(nodes, scores)
    owner = {}
    for c in order:
        owner[c] = c if c in cut.selected else owner.get(by_id[c]["parent_id"])
    pred = {p: owner[leaf] for p, leaf in membership}
    if len(pred) != len(membership) or any(v is None for v in pred.values()):
        raise ValueError("Incomplete conditioned cut")
    independent = conditioned_null(pred, species, edges)
    if not math.isclose(cut.score, independent["gain"], rel_tol=1e-8, abs_tol=1e-9):
        raise ValueError("Independent cut objective mismatch")
    return dict(
        algorithm="fixed-tree-species-pair-strength-null-v1",
        selected=list(cut.selected),
        score=cut.score,
        component_total_weight=total,
        members=[dict(protein_id=p, cluster_id=c) for p, c in sorted(pred.items())],
        node_scores=node_rows,
        reference_labels_used=False,
        tie_policy="keep_parent",
        objective_identity_verified=True,
        diagnostics=dict(groups=len(cut.selected), proteins=len(pred)),
    )
