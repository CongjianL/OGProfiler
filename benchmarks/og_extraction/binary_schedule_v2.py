"""Experimental ADR0004 v2: adaptive allocation within one 24-point budget."""

from __future__ import annotations

import math
from dataclasses import replace

from ogprofiler.hierarchy.resolution import ResolutionSearchResult

PROTOCOL = "adr0004-hard-binary-24-v2"


def finite_binary_search_v2(config, evaluate):
    config.validate()
    cap = min(24, config.max_candidate_evaluations)
    tested = []
    blocked = False

    def accepted():
        return sorted((c for c in tested if c.valid), key=lambda c: c.gamma)

    def sample(gamma, phase):
        nonlocal blocked
        if any(math.isclose(gamma, c.gamma, rel_tol=1e-12, abs_tol=0) for c in tested):
            return True
        if len(tested) >= cap:
            blocked = True
            return False
        c = evaluate(gamma)
        violations = c.violations + (() if c.child_count == 2 else ("TARGET_CHILD_COUNT",))
        tested.append(
            replace(
                c,
                valid=c.valid and c.child_count == 2,
                policy_valid=c.policy_valid and c.child_count == 2,
                violations=violations,
                rejection_reason=c.rejection_reason or (violations[0] if violations else None),
                phase=phase,
                evaluation_index=len(tested) + 1,
                evaluation_budget=cap,
            )
        )
        return True

    gamma = config.gamma_min
    for _ in range(8):
        if gamma > config.gamma_max or not sample(gamma, "COARSE") or accepted():
            break
        gamma *= config.growth_factor
    if not accepted() and not blocked:
        sample(config.gamma_max, "ENDPOINT")
    if not accepted() and not blocked:
        for i in range(1, 4):
            gamma = config.gamma_min * (config.gamma_max / config.gamma_min) ** (i / 4)
            if not sample(gamma, "RESCUE"):
                break

    boundary_turn = 0

    def next_point():
        nonlocal boundary_turn
        ordered = sorted(tested, key=lambda c: c.gamma)
        lower_edges, upper_edges, crossings, interiors, fallback = [], [], [], [], []
        for a, b in zip(ordered, ordered[1:], strict=False):
            midpoint = (a.gamma + b.gamma) / 2
            if any(math.isclose(midpoint, c.gamma, rel_tol=1e-12, abs_tol=0) for c in tested):
                continue
            item = ((-math.log(b.gamma / a.gamma), a.gamma, b.gamma), midpoint)
            if a.child_count != 2 and b.child_count == 2:
                lower_edges.append(item)
            elif a.child_count == 2 and b.child_count != 2:
                upper_edges.append(item)
            elif a.child_count < 2 <= b.child_count:
                crossings.append(item)
            elif a.child_count == 2 and b.child_count == 2:
                interiors.append(item)
            else:
                midpoint = math.sqrt(a.gamma * b.gamma)
                if not any(
                    math.isclose(midpoint, c.gamma, rel_tol=1e-12, abs_tol=0) for c in tested
                ):
                    fallback.append((item[0], midpoint))
        # Raw binary with rejected gates needs both boundaries: lower/upper round-robin.
        if lower_edges and upper_edges:
            pool = lower_edges if boundary_turn % 2 == 0 else upper_edges
            boundary_turn += 1
        else:
            pool = lower_edges or upper_edges or crossings or interiors or fallback
        return min(pool)[1] if pool else None

    # 9 nominal target points; after failure, up to all remaining shared slots are borrowed.
    for step in range(12):
        if accepted() or blocked:
            break
        gamma = next_point()
        if gamma is None:
            break
        if len(tested) == cap:
            if cap < 24:
                blocked = True
            # Completed finite schedule at full cap is rejection, not a fictitious 25th request.
            break
        if not sample(gamma, "TARGET_PROBE" if step < 9 else "BORROWED_REFINE"):
            break
    if accepted():
        upper = accepted()[0].gamma
        lower = max((c.gamma for c in tested if c.gamma < upper), default=config.gamma_min)
        for i in range(1, 4):
            if len(tested) == cap:
                break
            if not sample(lower + (upper - lower) * i / 4, "REFINE"):
                break
    chosen = accepted()
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
