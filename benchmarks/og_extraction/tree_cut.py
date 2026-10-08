"""Experimental complete tree cuts with externally supplied, label-free API scores.

No biological score is supplied here. Oracle scores belong only in benchmark code.
A network hierarchy cut is not, on its own, an ancestral OG inference.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass


@dataclass(frozen=True)
class Cut:
    selected: tuple[int, ...]
    score: float | int


def topology(nodes: Sequence[Mapping]):
    by_id = {n["cluster_id"]: n for n in nodes}
    if not nodes or len(by_id) != len(nodes):
        raise ValueError("Empty or duplicate node IDs")
    children = {c: [] for c in by_id}
    roots = []
    for c, n in by_id.items():
        if n.get("split_status") == "UNRESOLVED":
            raise ValueError("UNRESOLVED hierarchy is diagnostic only")
        p = n["parent_id"]
        if p is None:
            roots.append(c)
        elif p not in by_id or p == c:
            raise ValueError("Invalid parent")
        else:
            children[p].append(c)
    if len(roots) != 1:
        raise ValueError("Expected one root")
    order, seen, stack = [], set(), [roots[0]]
    while stack:
        c = stack.pop()
        if c in seen:
            raise ValueError("Cycle")
        seen.add(c)
        order.append(c)
        cs = sorted(children[c])
        if len(cs) == 1:
            raise ValueError("Unary nodes are unsupported")
        if "child_count" in by_id[c] and by_id[c]["child_count"] != len(cs):
            raise ValueError("Child-count mismatch")
        stack.extend(reversed(cs))
    if len(seen) != len(nodes):
        raise ValueError("Disconnected or cyclic hierarchy")
    return by_id, children, order


def optimal_cut(nodes, keep_scores, split_scores=None, eligible=None):
    """Maximize additive keep/split score; ties keep the parent deterministically.

    Every root-to-leaf path intersects exactly one selected node. An ineligible
    structural leaf has no feasible cover and raises instead of losing members.
    """
    by_id, children, order = topology(nodes)
    if set(keep_scores) != set(by_id):
        raise ValueError("Keep scores must cover exactly all nodes")
    split_scores = split_scores or {}
    if set(split_scores) - set(by_id):
        raise ValueError("Unknown split score node")
    if any(not math.isfinite(v) for v in (*keep_scores.values(), *split_scores.values())):
        raise ValueError("Non-finite score")
    allowed = set(by_id) if eligible is None else set(eligible)
    if allowed - set(by_id):
        raise ValueError("Unknown eligible node")
    best, keep = {}, {}
    for c in reversed(order):
        a = keep_scores[c] if c in allowed else None
        b = (
            split_scores.get(c, 0) + sum(best[u] for u in children[c])
            if children[c] and all(best[u] is not None for u in children[c])
            else None
        )
        if a is None and b is None:
            best[c], keep[c] = None, False
        elif a is not None and (b is None or a >= b):
            best[c], keep[c] = a, True
        else:
            best[c], keep[c] = b, False
    if best[order[0]] is None:
        raise ValueError("No complete eligible cut")
    selected, stack = [], [order[0]]
    while stack:
        c = stack.pop()
        if keep[c]:
            selected.append(c)
        else:
            stack.extend(reversed(sorted(children[c])))
    return Cut(tuple(selected), best[order[0]])


def pair_f1_oracle(forest, reference_pairs, max_iterations=100, *, restrict_eligible=False):
    """Exact finite fractional optimization over complete cuts, DIAGNOSTIC ONLY.

    Each tree supplies nodes and node-local tp/pp combinatorial counts. Maximize
    2*sum(tp)/(sum(pp)+reference_pairs) globally, not independent per-tree F1.
    Integer Dinkelbach comparisons avoid floating-point stopping ambiguity.
    This upper bound is specific to fixed-tree cuts and the supplied reference.
    With restrict_eligible, every tree must explicitly supply its eligible nodes;
    this remains a label-driven diagnostic, not the legacy active-view algorithm.
    """
    if reference_pairs <= 0:
        raise ValueError("Pair F1 oracle requires positive reference pair count")
    numerator, denominator = 0, 1
    for iteration in range(1, max_iterations + 1):
        cuts, tp, pp = [], 0, 0
        for tree in forest:
            if restrict_eligible and tree.get("eligible") is None:
                raise ValueError("Restricted oracle requires explicit eligible nodes")
            for c in tree["tp"]:
                if not 0 <= tree["tp"][c] <= tree["pp"][c]:
                    raise ValueError("Invalid pair counts")
            scores = {
                c: 2 * t * denominator - numerator * tree["pp"][c] for c, t in tree["tp"].items()
            }
            cut = optimal_cut(
                tree["nodes"],
                scores,
                eligible=tree["eligible"] if restrict_eligible else None,
            )
            cuts.append(cut.selected)
            tp += sum(tree["tp"][c] for c in cut.selected)
            pp += sum(tree["pp"][c] for c in cut.selected)
        if tp > reference_pairs:
            raise ValueError("Selected true pairs exceed reference total")
        new_num, new_den = 2 * tp, pp + reference_pairs
        gain = new_num * denominator - numerator * new_den
        if gain < 0:
            raise AssertionError("Fractional optimization regressed")
        if gain == 0:
            return dict(
                cuts=cuts,
                true_positive_pairs=tp,
                predicted_pairs=pp,
                reference_pairs=reference_pairs,
                pair_f1=new_num / new_den,
                iterations=iteration,
                diagnostic_only=True,
            )
        numerator, denominator = new_num, new_den
    raise ValueError("Oracle did not converge within the explicit iteration limit")
