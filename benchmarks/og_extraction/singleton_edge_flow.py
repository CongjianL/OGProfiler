"""Raw singleton edges and legal three-child cuts; read-only diagnostic."""

from collections import Counter
from fractions import Fraction

from benchmarks.og_extraction.tree_cut import topology


def trace_flow(nodes, membership, source_edges, incident_edges, species, node_id):
    by_id, children, order = topology(nodes)
    proteins = {c: set() for c in order}
    for p, leaf in membership:
        proteins[leaf].add(p)
    for c in reversed(order):
        for ch in children[c]:
            proteins[c].update(proteins[ch])
    cs = sorted(children[node_id])
    if len(cs) != 3:
        raise ValueError("Expected three immediate children")
    singletons = {next(iter(proteins[c])): c for c in cs if len(proteins[c]) == 1}
    if len(singletons) != 2:
        raise ValueError("Expected two singleton children")
    source = proteins[order[0]]
    node = proteins[node_id]
    child_owner = {p: c for c in cs for p in proteins[c]}
    seen = set()
    raw = []
    incident_in_source = {}
    sums = Counter()
    for u, v, w in incident_edges:
        if u >= v or (u, v) in seen or (u not in singletons and v not in singletons):
            raise ValueError("Invalid incident edge")
        seen.add((u, v))
        weight = Fraction.from_float(float(w))
        if weight < 0:
            raise ValueError("Negative weight")
        if u in source and v in source:
            incident_in_source[u, v] = weight
        for p, other in ((u, v), (v, u)):
            if p not in singletons:
                continue
            scope = (
                "node_child"
                if other in node
                else "outside_node_inside_source"
                if other in source
                else "outside_source"
            )
            row = dict(
                singleton=p,
                other=other,
                u=u,
                v=v,
                weight=float(weight),
                singleton_species=species[p],
                other_species=species[other],
                other_child=child_owner.get(other),
                scope=scope,
            )
            raw.append(row)
            sums[p, scope, species[other], child_owner.get(other)] += weight
    expected_incident = {
        (u, v): Fraction.from_float(float(w))
        for u, v, w in source_edges
        if u in singletons or v in singletons
    }
    if expected_incident != incident_in_source:
        raise ValueError("Source induced edges and raw incident edges disagree")
    blocks, strengths, actual = Counter(), Counter(), Counter()
    for u, v, w in source_edges:
        w = Fraction.from_float(float(w))
        s, t = species[u], species[v]
        blocks[min(s, t), max(s, t)] += w
        if u in node:
            strengths[child_owner[u], s, t] += w
        if v in node:
            strengths[child_owner[v], t, s] += w
        if u in node and v in node and child_owner[u] != child_owner[v]:
            a, b = sorted((child_owner[u], child_owner[v]))
            actual[a, b, min(s, t), max(s, t)] += w
    total = sum(blocks.values(), Fraction())
    pair_rows = []
    for i, a in enumerate(cs):
        for b in cs[i + 1 :]:
            for (s, t), w in sorted(blocks.items()):
                if not w:
                    continue
                expected = (
                    strengths[a, s, s] * strengths[b, s, s] / (2 * w)
                    if s == t
                    else (
                        strengths[a, s, t] * strengths[b, t, s]
                        + strengths[b, s, t] * strengths[a, t, s]
                    )
                    / w
                )
                observed = actual[a, b, s, t]
                if expected or observed:
                    delta = expected - observed
                    pair_rows.append(
                        dict(
                            a=a,
                            b=b,
                            species_a=s,
                            species_b=t,
                            expected=float(expected),
                            observed=float(observed),
                            weight_deficit=float(delta),
                            deficit_sign=1 if delta > 0 else -1 if delta < 0 else 0,
                            gain=float(delta / total) if total else 0.0,
                        )
                    )
    # All five partitions of three children: only actual node unions are tree nodes.
    partitions = [(cs,), tuple((c,) for c in cs)]
    partitions += [(tuple(c for c in cs if c != separate), (separate,)) for separate in cs]
    candidates = []
    for parts in partitions:
        group_owner = {c: i for i, part in enumerate(parts) for c in part}
        gain = Fraction()
        for i, a in enumerate(cs):
            for b in cs[i + 1 :]:
                if group_owner[a] == group_owner[b]:
                    continue
                for (s, t), w in blocks.items():
                    if w:
                        expected = (
                            strengths[a, s, s] * strengths[b, s, s] / (2 * w)
                            if s == t
                            else (
                                strengths[a, s, t] * strengths[b, t, s]
                                + strengths[b, s, t] * strengths[a, t, s]
                            )
                            / w
                        )
                        gain += expected - actual[a, b, s, t]
        gain = gain / total if total else Fraction()
        unions = [set().union(*(proteins[c] for c in part)) for part in parts]
        union_nodes = [next((c for c in order if proteins[c] == ps), None) for ps in unions]
        candidates.append(
            dict(
                parts=[list(part) for part in parts],
                prediction={p: i for i, ps in enumerate(unions) for p in ps},
                node_ids=union_nodes,
                tree_representable=all(c is not None for c in union_nodes),
                gain_vs_keep=float(gain),
                gain_sign=1 if gain > 0 else -1 if gain < 0 else 0,
            )
        )
    return dict(
        node_id=node_id,
        source_root=order[0],
        direct_children=cs,
        child_sizes={c: len(proteins[c]) for c in cs},
        singletons=singletons,
        incident_source_identity_verified=True,
        incident_edges=raw,
        edge_sums=[
            dict(singleton=p, scope=sc, other_species=s, other_child=c, weight=float(w))
            for (p, sc, s, c), w in sums.items()
        ],
        child_pair_species_blocks=pair_rows,
        candidate_partitions=candidates,
        source_total_weight=float(total),
        reference_labels_used_for_structure=False,
    )
