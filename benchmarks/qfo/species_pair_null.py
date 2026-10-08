"""Species-pair-conditioned null on saved cuts; no selection or reference labels."""

import math
from collections import Counter


def conditioned_null(prediction, species, edges):
    """Preserve each species block's weight and per-group endpoint strengths.

    Cross-species blocks use a bipartite configuration expectation. Same-species
    blocks use an undirected strength null (including expected self loops).
    This is a sensitivity diagnostic, not an optimized partition objective.
    Edge validity is checked by the caller's unconditional decomposition.
    """
    weights = Counter()
    strengths = Counter()
    observed = Counter()
    for u, v, w in edges:
        s, t = species[u], species[v]
        block = min(s, t), max(s, t)
        weights[block] += w
        strengths[prediction[u], s, t] += w
        strengths[prediction[v], t, s] += w
        if prediction[u] != prediction[v]:
            observed[block] += w
    groups = sorted(set(prediction.values()))
    total = sum(weights.values())
    rows = []
    terms = Counter()
    for (s, t), weight in sorted(weights.items()):
        if s == t:
            values = [strengths[g, s, s] for g in groups]
            expected = (
                (sum(values) ** 2 - sum(v * v for v in values)) / (4 * weight) if weight else 0.0
            )
            endpoint_totals = [sum(values)]
            target_totals = [2 * weight]
        else:
            left = [strengths[g, s, t] for g in groups]
            right = [strengths[g, t, s] for g in groups]
            expected = (
                (sum(left) * sum(right) - sum(a * b for a, b in zip(left, right, strict=True)))
                / weight
                if weight
                else 0.0
            )
            endpoint_totals = [sum(left), sum(right)]
            target_totals = [weight, weight]
        if any(
            not math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-9)
            for a, b in zip(endpoint_totals, target_totals, strict=True)
        ):
            raise ValueError("Species block strength conservation failed")
        kind = "same_species" if s == t else "cross_species"
        contribution = (expected - observed[s, t]) / total if total else 0.0
        terms[kind] += contribution
        rows.append(
            dict(
                species_a=s,
                species_b=t,
                total_weight=weight,
                observed_between_groups_weight=observed[s, t],
                expected_between_groups_weight=expected,
                gain_contribution=contribution,
            )
        )
    return dict(
        gain=sum(terms.values()),
        same_species_gain=terms["same_species"],
        cross_species_gain=terms["cross_species"],
        block_conservation_verified=True,
        species_blocks=rows,
    )
