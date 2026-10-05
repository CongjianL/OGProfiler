"""Label-free weighted graph cut on a fixed hierarchy; experimental, not ancestral OG.

Objective: sum_v [internal_weight(v) - pair_penalty * n(v)*(n(v)-1)/2]
over a complete, non-overlapping cut. All resolved nodes are candidates. No
reference labels, event qualification, Leiden calls, or production files are used.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from benchmarks.og_extraction.tree_cut import optimal_cut, topology
from ogprofiler.core.manifest import sha256_file

ALGORITHM = "fixed-tree-weighted-pair-penalty-v1"


STRENGTH_ALGORITHM = "fixed-tree-component-strength-null-v1"


def graph_cut(nodes, membership, edges, *, pair_penalty=None, strength_null=False):
    """Score one component. Edges are unique (u, v, weight), with u < v.

    Membership is (protein_id, structural_leaf_id), including isolated proteins.
    Weights and penalty must be finite and nonnegative. Penalty has edge-weight
    units and is explicit: zero/large values may produce root/leaf degeneracy.
    Structural leaves remain indivisible even when their score is negative.
    Runtime O(nodes + edges * tree_height); no internal descendant lists.
    """
    if strength_null:
        if pair_penalty is not None:
            raise ValueError("Strength-null has fixed resolution 1, no pair_penalty")
    elif (
        pair_penalty is None
        or isinstance(pair_penalty, bool)
        or not math.isfinite(pair_penalty)
        or pair_penalty < 0
    ):
        raise ValueError("pair_penalty must be finite and nonnegative")
    by_id, children, order = topology(nodes)
    leaf_by_protein, sizes = {}, Counter()
    for p, leaf in membership:
        if p in leaf_by_protein or leaf not in by_id or children[leaf]:
            raise ValueError("Duplicate protein or invalid terminal membership")
        leaf_by_protein[p] = leaf
        sizes[leaf] += 1
    depth = {}
    for c in order:
        parent = by_id[c]["parent_id"]
        depth[c] = 0 if parent is None else depth[parent] + 1
    for c in reversed(order):
        sizes[c] += sum(sizes[u] for u in children[c])
        if sizes[c] == 0 or ("n_genes" in by_id[c] and sizes[c] != by_id[c]["n_genes"]):
            raise ValueError("Empty node or node membership count mismatch")

    # Each edge contributes once at its LCA, then propagates upward by postorder.
    at_lca, seen_edges = defaultdict(list), set()
    terminal_strength = defaultdict(list)
    for u, v, weight in edges:
        if u not in leaf_by_protein or v not in leaf_by_protein:
            raise ValueError("Edge endpoint outside membership")
        if u >= v or (u, v) in seen_edges:
            raise ValueError("Expected unique canonical edges u < v")
        if isinstance(weight, bool) or not math.isfinite(weight) or weight < 0:
            raise ValueError("Edge weight must be finite and nonnegative")
        seen_edges.add((u, v))
        a, b = leaf_by_protein[u], leaf_by_protein[v]
        terminal_strength[a].append(weight)
        terminal_strength[b].append(weight)
        while depth[a] > depth[b]:
            a = by_id[a]["parent_id"]
        while depth[b] > depth[a]:
            b = by_id[b]["parent_id"]
        while a != b:
            a, b = by_id[a]["parent_id"], by_id[b]["parent_id"]
        at_lca[a].append(weight)
    internal, scores, strength = {}, {}, {}
    for c in reversed(order):
        internal[c] = math.fsum([*at_lca[c], *(internal[u] for u in sorted(children[c]))])
        strength[c] = math.fsum(
            [*terminal_strength[c], *(strength[u] for u in sorted(children[c]))]
        )
    total_weight = internal[order[0]]
    if not math.isfinite(total_weight):
        raise ValueError("Non-finite total weight")
    for c in order:
        if strength_null:
            scores[c] = (
                internal[c] / total_weight - (strength[c] / (2 * total_weight)) ** 2
                if total_weight
                else 0.0
            )
        else:
            scores[c] = internal[c] - pair_penalty * (sizes[c] * (sizes[c] - 1) // 2)
    cut = optimal_cut([by_id[c] for c in sorted(by_id)], scores)
    chosen = set(cut.selected)
    owner = {}
    for c in order:
        owner[c] = c if c in chosen else owner.get(by_id[c]["parent_id"])
    members = [
        dict(protein_id=p, cluster_id=owner[leaf]) for p, leaf in sorted(leaf_by_protein.items())
    ]
    if any(row["cluster_id"] is None for row in members):
        raise ValueError("Incomplete cut")
    return dict(
        algorithm=STRENGTH_ALGORITHM if strength_null else ALGORITHM,
        resolution=1.0 if strength_null else None,
        component_total_weight=total_weight,
        reference_labels_used=False,
        ancestral_inference=False,
        pair_penalty=pair_penalty,
        weight_normalization="component-total-weight objective" if strength_null else "none",
        tie_policy="keep_parent",
        selected=list(cut.selected),
        score=cut.score,
        members=members,
        node_scores=[
            dict(
                cluster_id=c,
                n_genes=sizes[c],
                internal_weight=internal[c],
                node_strength=strength[c],
                keep_score=scores[c],
            )
            for c in sorted(by_id)
        ],
        diagnostics=dict(
            root_selected=order[0] in chosen,
            all_structural_leaves_selected=all(not children[c] for c in chosen),
            groups=len(chosen),
            proteins=len(members),
            edges=len(seen_edges),
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="JSON component: nodes, membership [[protein,leaf]], edges [[u,v,w]]",
    )
    parser.add_argument("--pair-penalty", type=float, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    digest = sha256_file(args.input)
    data = json.loads(args.input.read_text())
    if set(data) != {"nodes", "membership", "edges"}:
        raise ValueError("Expected exactly nodes, membership, edges; no reference labels")
    report = graph_cut(**data, pair_penalty=args.pair_penalty)
    if digest != sha256_file(args.input):
        raise ValueError("Input changed during execution")
    report["input_sha256"] = digest
    report["source_sha256"] = sha256_file(Path(__file__))
    report["tree_cut_source_sha256"] = sha256_file(Path(__file__).with_name("tree_cut.py"))
    with args.out.open("x") as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write("\n")


if __name__ == "__main__":
    main()
