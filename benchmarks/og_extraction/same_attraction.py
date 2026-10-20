"""One frozen, label-free additive intervention; no production/default changes."""

from collections import Counter
from fractions import Fraction

from benchmarks.og_extraction.exact_species_pair_graph import exact_context
from benchmarks.og_extraction.fixed_tree_graph import graph_cut
from benchmarks.og_extraction.tree_cut import optimal_cut, topology


def dp_state(nodes, scores):
    by_id, children, order = topology(nodes)
    best, keep = {}, {}
    for c in reversed(order):
        split = sum((best[ch] for ch in children[c]), Fraction())
        keep[c] = not children[c] or scores[c] >= split
        best[c] = scores[c] if keep[c] else split
    selected, ancestor, owner = set(), {}, {}
    for c in order:
        p = by_id[c]["parent_id"]
        ancestor[c] = owner.get(p)
        if ancestor[c] is None and keep[c]:
            selected.add(c)
            owner[c] = c
        else:
            owner[c] = ancestor[c]
    return best, keep, selected, ancestor, owner


def independent_gain(prediction, edges, species):
    """Candidate partition minus source merge, directly from raw edge endpoints."""
    weights, strength = Counter(), Counter()
    observed = Fraction()
    groups = set(prediction.values())
    for u, v, w in edges:
        ew = Fraction.from_float(float(w))
        s, t = species[u], species[v]
        weights[min(s, t), max(s, t)] += ew
        strength[prediction[u], s, t] += ew
        strength[prediction[v], t, s] += ew
        if prediction[u] != prediction[v]:
            observed += ew
    expected = Fraction()
    for (s, t), w in weights.items():
        if s == t or not w:
            continue
        a = [strength[g, s, t] for g in groups]
        b = [strength[g, t, s] for g in groups]
        expected += (sum(a) * sum(b) - sum(x * y for x, y in zip(a, b, strict=True))) / w
    total = sum(weights.values(), Fraction())
    return (expected - observed) / total if total else Fraction()


def intervene(nodes, membership, edges, species, baseline_selected):
    # Reuse established validation of graph endpoints, weights and terminal coverage.
    graph_cut(nodes, membership, edges, strength_null=True)
    original, total, strengths, weights = exact_context(nodes, membership, edges, species)
    by_id, children, order = topology(nodes)
    same_expectation = {
        c: sum(
            (strengths[c][s, s] ** 2 / (4 * w) for (s, t), w in weights.items() if s == t and w),
            Fraction(),
        )
        for c in order
    }
    candidate = {
        c: original[c] + same_expectation[c] / total if total else Fraction() for c in order
    }
    bb, bk, bs, ba, bo = dp_state(nodes, original)
    cb, ck, cs, ca, co = dp_state(nodes, candidate)
    if bs != set(baseline_selected):
        raise ValueError("Intervention baseline differs from frozen exact cut")
    optimized = optimal_cut(nodes, candidate)
    if set(optimized.selected) != cs or optimized.score != cb[order[0]]:
        raise ValueError("Candidate DP differs from standard complete-cut solver")
    # Same-null expectation is superadditive on disjoint descendants. Thus every
    # fixed refinement gains no more than before, and an old local keep stays keep.
    if any(bk[c] and not ck[c] for c in order):
        raise ValueError("Candidate unexpectedly refines an original local keep")
    root = order[0]
    expected_root = (
        sum((w for (s, t), w in weights.items() if s == t), Fraction()) / total
        if total
        else Fraction()
    )
    if candidate[root] != expected_root:
        raise ValueError("Candidate nonzero root baseline mismatch")
    predictions = [{p: owner[leaf] for p, leaf in membership} for owner in (bo, co)]
    if any(v is None for pred in predictions for v in pred.values()):
        raise ValueError("Incomplete candidate coverage")
    if any(co[bo[leaf]] != co[leaf] for _, leaf in membership):
        raise ValueError("Candidate is not a coarsening of the original complete cut")
    candidate_baseline = sum((candidate[c] for c in bs), Fraction())
    original_candidate = sum((original[c] for c in cs), Fraction())
    if cb[root] < candidate_baseline or bb[root] < original_candidate:
        raise ValueError("Exact optimization violated feasible-cut ordering")
    for chosen, pred in zip((bs, cs), predictions, strict=True):
        score = sum((candidate[c] for c in chosen), Fraction()) - candidate[root]
        if independent_gain(pred, edges, species) != score:
            raise ValueError("Independent candidate boundary/root identity failed")
        if score > sum((original[c] for c in chosen), Fraction()):
            raise ValueError("Fixed partition candidate gain exceeds original gain")
    rows = []
    for c in order:
        if not children[c]:
            continue
        direct_b = sum((original[ch] for ch in children[c]), Fraction()) - original[c]
        direct_c = sum((candidate[ch] for ch in children[c]), Fraction()) - candidate[c]
        removed = (
            (same_expectation[c] - sum(same_expectation[ch] for ch in children[c])) / total
            if total
            else Fraction()
        )
        if removed < 0 or direct_b - direct_c != removed:
            raise ValueError("Fixed direct-boundary same-null ablation identity failed")
        state_b = "inactive" if ba[c] is not None else "active_keep" if c in bs else "active_split"
        state_c = "inactive" if ca[c] is not None else "active_keep" if c in cs else "active_split"
        path = []
        p = c
        while p is not None:
            path.append(p)
            p = by_id[p]["parent_id"]
        rows.append(
            dict(
                cluster_id=c,
                parent_id=by_id[c]["parent_id"],
                baseline_state=state_b,
                candidate_state=state_c,
                baseline_local_keep=bk[c],
                candidate_local_keep=ck[c],
                baseline_selected_ancestor=ba[c],
                candidate_selected_ancestor=ca[c],
                baseline_direct_gain=float(direct_b),
                candidate_direct_gain=float(direct_c),
                baseline_optimized_margin=float(
                    sum((bb[ch] for ch in children[c]), Fraction()) - original[c]
                ),
                candidate_optimized_margin=float(
                    sum((cb[ch] for ch in children[c]), Fraction()) - candidate[c]
                ),
                removed_same_null_gain=float(removed),
                root_path=list(reversed(path)) if state_b != state_c or bk[c] != ck[c] else [],
                local_decision_changed=bk[c] != ck[c],
                active_state_changed=state_b != state_c,
                ancestor_only_state_change=state_b != state_c and bk[c] == ck[c],
            )
        )
    return dict(
        algorithm="same_attraction_cross_null_v1",
        reference_labels_used=False,
        source_weight=float(total),
        root_baseline=float(candidate[root]),
        root_baseline_exact=str(candidate[root]),
        candidate_score=float(cb[root]),
        candidate_gain_over_merge=float(cb[root] - candidate[root]),
        baseline_gain_original=float(bb[root]),
        baseline_gain_candidate=float(candidate_baseline - candidate[root]),
        candidate_gain_original=float(original_candidate),
        baseline_prediction=predictions[0],
        candidate_prediction=predictions[1],
        baseline_selected=sorted(bs),
        candidate_selected=sorted(cs),
        node_transitions=rows,
        objective_identity_verified=True,
        coarsening_verified=True,
        exact_tie_policy="keep parent",
        production_changed=False,
    )
