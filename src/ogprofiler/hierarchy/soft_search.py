"""ADR0004 opt-in soft binary search within one shared finite candidate budget.

The evaluator applies original gates. Fallback only reuses original-gated points;
count/stability are not assumed monotonic. No I/O, recursion or OG logic here.
"""

from __future__ import annotations

import math
from dataclasses import replace

from ogprofiler.hierarchy.resolution import ResolutionSearchResult


def search_soft_binary_24(config, evaluate):
    config.validate()
    cap = min(24, config.max_candidate_evaluations)
    tested = []
    originals = {}
    blocked = False

    def accepted():
        return sorted((c for c in tested if c.valid), key=lambda c: c.gamma)

    def sample(gamma, phase):
        nonlocal blocked
        for c in tested:
            if math.isclose(c.gamma, gamma, rel_tol=1e-12, abs_tol=0):
                return c
        if len(tested) >= cap:
            blocked = True
            return None
        c = evaluate(gamma)
        originals[gamma] = c
        violations = c.violations + (() if c.child_count == 2 else ("TARGET_CHILD_COUNT",))
        c = replace(
            c,
            valid=c.valid and c.child_count == 2,
            policy_valid=c.policy_valid and c.child_count == 2,
            violations=violations,
            rejection_reason=c.rejection_reason or (violations[0] if violations else None),
            phase=phase,
            evaluation_index=len(tested) + 1,
            evaluation_budget=cap,
        )
        tested.append(c)
        return c

    def has_target_evidence():
        ordered = sorted(tested, key=lambda c: c.gamma)
        return any(c.child_count == 2 for c in ordered) or any(
            a.child_count < 2 <= b.child_count for a, b in zip(ordered, ordered[1:], strict=False)
        )

    gamma = config.gamma_min
    for _ in range(min(8, config.max_coarse_candidates)):
        if gamma > config.gamma_max:
            break
        c = sample(gamma, "COARSE")
        if c is None or c.child_count >= 2:
            break
        gamma *= config.growth_factor
    # A versioned guard: rejected coarse binary needs a measured upper side.
    rejected_split = (
        tested and tested[-1].child_count >= 2 and not originals[tested[-1].gamma].valid
    )
    if rejected_split:
        gamma = tested[-1].gamma
        for _ in range(cap):
            if accepted() or blocked or len(tested) >= max(1, cap - 3):
                break
            upper = min(gamma * config.growth_factor, config.gamma_max)
            if math.isclose(upper, gamma, rel_tol=1e-12, abs_tol=0):
                break
            c = sample(upper, "UPPER_GUARD")
            if c is None:
                break
            gamma = upper
            if c.child_count > 2 and originals[upper].valid:
                break
    # Global points are conditional: do not spend them after a target bracket is established.
    if not accepted() and not blocked and not has_target_evidence():
        sample(config.gamma_max, "ENDPOINT")
        if not accepted() and not blocked:
            for i in range(1, min(3, config.rescue_grid_points) + 1):
                gamma = config.gamma_min * (config.gamma_max / config.gamma_min) ** (i / 4)
                if sample(gamma, "RESCUE") is None:
                    break
    boundary_turn = 0

    def next_point():
        nonlocal boundary_turn
        lower, upper, crossings, interiors, fallback = [], [], [], [], []
        ordered = sorted(tested, key=lambda c: c.gamma)
        for a, b in zip(ordered, ordered[1:], strict=False):
            midpoint = (a.gamma + b.gamma) / 2
            if any(math.isclose(midpoint, c.gamma, rel_tol=1e-12, abs_tol=0) for c in tested):
                continue
            rank = (-math.log(b.gamma / a.gamma), a.gamma, b.gamma)
            item = (rank, midpoint)
            if a.child_count != 2 and b.child_count == 2:
                lower.append(item)
            elif a.child_count == 2 and b.child_count != 2:
                upper.append(item)
            elif a.child_count < 2 <= b.child_count:
                crossings.append(item)
            elif a.child_count == b.child_count == 2:
                interiors.append(item)
            else:
                midpoint = math.sqrt(a.gamma * b.gamma)
                if not any(
                    math.isclose(midpoint, c.gamma, rel_tol=1e-12, abs_tol=0) for c in tested
                ):
                    fallback.append((rank, midpoint))
        if lower and upper:
            pool = lower if boundary_turn % 2 == 0 else upper
            boundary_turn += 1
        else:
            pool = lower or upper or crossings or interiors or fallback
        phase = "TARGET_PROBE"
        return (min(pool)[1], phase) if pool else None

    while not accepted() and not blocked:
        point = next_point()
        if point is None:
            break
        if len(tested) == cap:
            blocked = cap < 24
            break
        gamma, phase = point
        if len(tested) >= cap - 3:
            phase = "BORROWED_" + phase
        if sample(gamma, phase) is None:
            break
    if accepted():
        upper = accepted()[0].gamma
        lower = max((c.gamma for c in tested if c.gamma < upper), default=config.gamma_min)
        for i in range(1, 4):
            if len(tested) == cap:
                break
            if sample(lower + (upper - lower) * i / 4, "REFINE") is None:
                break
    chosen = accepted()
    if not chosen and not blocked:
        fallback = sorted(
            (c for c in tested if originals[c.gamma].valid and c.child_count > 2),
            key=lambda c: c.gamma,
        )
        if fallback:
            prior = fallback[0]
            original = originals[prior.gamma]
            selected = replace(
                prior,
                valid=True,
                policy_valid=original.policy_valid,
                violations=original.violations,
                rejection_reason=original.rejection_reason,
                phase="FALLBACK_KWAY/" + prior.phase,
                selection_kind="FALLBACK_KWAY",
            )
            tested[tested.index(prior)] = selected
            chosen = [selected]
    if chosen and chosen[0].child_count == 2:
        prior = chosen[0]
        selected = replace(
            prior,
            selection_kind="BINARY",
            refinement_truncated=len(tested) == cap
            and sum(c.phase == "REFINE" for c in tested) < 3,
        )
        tested[tested.index(prior)] = selected
        chosen = [selected]
    status = (
        "ACCEPTED"
        if chosen
        else "EVALUATION_BUDGET_EXHAUSTED"
        if blocked
        else "REJECTED_ALL_TESTED"
    )
    return ResolutionSearchResult(
        chosen[0] if chosen else None, tuple(tested), None if chosen else status, status
    )
