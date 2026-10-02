"""Resolution search with split acceptance and multi-seed stability."""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass, replace

import igraph as ig

from ogprofiler.exceptions import HierarchyError
from ogprofiler.hierarchy.leiden import LeidenCallCounter, LeidenResult, run_leiden


@dataclass(frozen=True, slots=True)
class ResolutionSearchConfig:
    strategy: str = "bounded_adaptive_v2"
    admission_policy: str = "nonempty_children_v1"
    max_candidate_evaluations: int = 24
    max_coarse_candidates: int = 10
    rescue_grid_points: int = 8
    gamma_min: float = 0.01
    gamma_max: float = 10.0
    growth_factor: float = 2.0
    local_grid_points: int = 5
    min_child_size: int = 2
    max_child_fraction: float = 0.95
    tiny_fragment_size: int = 2
    max_tiny_fragment_fraction: float = 1.0
    stability_threshold: float = 0.9
    publication_seeds: int = 5
    min_quality: float | None = None

    def validate(self) -> None:
        if self.strategy not in {"adaptive", "log_grid", "bounded_adaptive_v2"}:
            raise HierarchyError(f"Unknown resolution strategy: {self.strategy}")
        if self.admission_policy not in {"legacy_strict", "nonempty_children_v1"}:
            raise HierarchyError("Unknown admission policy")
        if (
            self.strategy == "bounded_adaptive_v2"
            and self.admission_policy != "nonempty_children_v1"
        ):
            raise HierarchyError("bounded_adaptive_v2 requires nonempty_children_v1")
        if (
            self.admission_policy == "nonempty_children_v1"
            and self.strategy != "bounded_adaptive_v2"
        ):
            raise HierarchyError("nonempty_children_v1 requires bounded_adaptive_v2")
        if self.min_quality is not None and not math.isfinite(self.min_quality):
            raise HierarchyError("min_quality must be finite")
        if any(
            type(value) is not int or value < 1
            for value in (
                self.max_candidate_evaluations,
                self.max_coarse_candidates,
                self.rescue_grid_points,
            )
        ):
            raise HierarchyError("Search budgets must be positive integers")
        if not all(
            math.isfinite(value) for value in (self.gamma_min, self.gamma_max, self.growth_factor)
        ):
            raise HierarchyError("Search bounds must be finite")
        if self.gamma_min <= 0 or self.gamma_max < self.gamma_min:
            raise HierarchyError("Resolution bounds must satisfy 0 < gamma_min <= gamma_max")
        if self.growth_factor <= 1 or self.local_grid_points < 2:
            raise HierarchyError("Resolution search configuration is invalid")
        if self.min_child_size < 1 or self.tiny_fragment_size < 1:
            raise HierarchyError("Child sizes must be positive")
        if not 0 < self.max_child_fraction < 1:
            raise HierarchyError("max_child_fraction must be between 0 and 1")
        if not 0 <= self.max_tiny_fragment_fraction <= 1:
            raise HierarchyError("max_tiny_fragment_fraction must be between 0 and 1")
        if not 0 <= self.stability_threshold <= 1:
            raise HierarchyError("stability_threshold must be between 0 and 1")
        if self.publication_seeds < 5:
            raise HierarchyError("publication_seeds must be at least 5")


@dataclass(frozen=True, slots=True)
class SplitCandidate:
    gamma: float
    membership: tuple[int, ...]
    child_count: int
    quality: float
    min_child_size: int
    max_child_fraction: float
    tiny_fragment_fraction: float
    stability: float
    adjusted_rand_index: float
    normalized_mutual_info: float
    inter_edge_fraction: float
    intra_edge_fraction: float
    valid: bool
    rejection_reason: str | None
    violations: tuple[str, ...] = ()
    structural_valid: bool = True
    policy_valid: bool = True
    phase: str = "legacy"
    evaluation_index: int = 0
    evaluation_budget: int = 0
    stability_evaluated: bool = True


@dataclass(frozen=True, slots=True)
class ResolutionSearchResult:
    selected: SplitCandidate | None
    candidates: tuple[SplitCandidate, ...]
    terminal_reason: str | None
    search_status: str = "REJECTED_ALL_TESTED"


def _canonical(membership: tuple[int, ...]) -> tuple[int, ...]:
    first: dict[int, int] = {}
    for index, label in enumerate(membership):
        first.setdefault(label, index)
    mapping = {label: rank for rank, label in enumerate(sorted(first, key=first.__getitem__))}
    return tuple(mapping[label] for label in membership)


def adjusted_rand_index(left: tuple[int, ...], right: tuple[int, ...]) -> float:
    if len(left) != len(right):
        raise HierarchyError("Partition sizes differ")
    if len(left) < 2:
        return 1.0

    def choose2(value: int) -> float:
        return value * (value - 1) / 2

    contingency = Counter(zip(left, right, strict=True))
    left_pairs = sum(choose2(value) for value in Counter(left).values())
    right_pairs = sum(choose2(value) for value in Counter(right).values())
    both = sum(choose2(value) for value in contingency.values())
    total = choose2(len(left))
    expected = left_pairs * right_pairs / total
    denominator = (left_pairs + right_pairs) / 2 - expected
    return 1.0 if math.isclose(denominator, 0.0) else (both - expected) / denominator


def normalized_mutual_info(left: tuple[int, ...], right: tuple[int, ...]) -> float:
    if len(left) != len(right):
        raise HierarchyError("Partition sizes differ")
    if not left:
        return 1.0
    n = len(left)
    left_counts, right_counts = Counter(left), Counter(right)
    mutual = sum(
        (count / n) * math.log((count * n) / (left_counts[left_label] * right_counts[right_label]))
        for (left_label, right_label), count in Counter(zip(left, right, strict=True)).items()
    )
    h_left = -sum((count / n) * math.log(count / n) for count in left_counts.values())
    h_right = -sum((count / n) * math.log(count / n) for count in right_counts.values())
    return 1.0 if math.isclose(h_left + h_right, 0.0) else 2 * mutual / (h_left + h_right)


def _seed_count(mode: str, publication_seeds: int) -> int:
    values = {"fast": 1, "robust": 3, "publication": publication_seeds}
    if mode not in values:
        raise HierarchyError(f"Unknown stability mode: {mode}")
    return values[mode]


def _multi_seed(
    graph: ig.Graph,
    gamma: float,
    method: str,
    weights: str | list[float] | None,
    seed: int,
    mode: str,
    publication_seeds: int,
    counter: LeidenCallCounter,
) -> tuple[LeidenResult, float, float, float]:
    results = [
        run_leiden(graph, gamma, method, weights, seed + index * 104_729, counter)
        for index in range(_seed_count(mode, publication_seeds))
    ]
    for result in results:
        if (
            len(result.membership) != graph.vcount()
            or not all(type(label) is int and label >= 0 for label in result.membership)
            or not math.isfinite(result.quality)
        ):
            raise HierarchyError("COMPUTATION_ERROR: invalid Leiden partition or quality")
    memberships = [_canonical(result.membership) for result in results]
    if len(results) == 1:
        return LeidenResult(memberships[0], results[0].quality), 1.0, 1.0, 1.0
    ari_values: list[float] = []
    nmi_values: list[float] = []
    support = [0.0] * len(results)
    for left in range(len(results)):
        for right in range(left + 1, len(results)):
            ari = adjusted_rand_index(memberships[left], memberships[right])
            nmi = normalized_mutual_info(memberships[left], memberships[right])
            ari_values.append(ari)
            nmi_values.append(nmi)
            support[left] += ari
            support[right] += ari
    best = max(range(len(results)), key=lambda index: (support[index], results[index].quality))
    mean_ari = sum(ari_values) / len(ari_values)
    return (
        LeidenResult(memberships[best], results[best].quality),
        mean_ari,
        mean_ari,
        sum(nmi_values) / len(nmi_values),
    )


def _candidate(
    graph: ig.Graph,
    gamma: float,
    result: LeidenResult,
    stability: float,
    ari: float,
    nmi: float,
    config: ResolutionSearchConfig,
) -> SplitCandidate:
    if len(result.membership) != graph.vcount() or not result.membership:
        raise HierarchyError("Partition membership must cover every vertex exactly once")
    if not all(type(label) is int and label >= 0 for label in result.membership):
        raise HierarchyError("Invalid partition labels")
    if not all(math.isfinite(value) for value in (result.quality, stability, ari, nmi)):
        raise HierarchyError("Non-finite partition quality or stability")
    sizes = tuple(Counter(result.membership).values())
    total = len(result.membership)
    minimum = min(sizes)
    maximum_fraction = max(sizes) / total
    tiny_fraction = sum(size for size in sizes if size <= config.tiny_fragment_size) / total
    inter = sum(
        1
        for left, right in graph.get_edgelist()
        if result.membership[left] != result.membership[right]
    )
    inter_fraction = inter / graph.ecount() if graph.ecount() else 0.0
    violations: list[str] = []
    if len(sizes) <= 1:
        violations.append("NO_SPLIT")
    if config.admission_policy == "legacy_strict" and minimum < config.min_child_size:
        violations.append("MIN_CHILD_SIZE")
    if maximum_fraction > config.max_child_fraction:
        violations.append("MAX_CHILD_FRACTION")
    if tiny_fraction > config.max_tiny_fragment_fraction:
        violations.append("TINY_FRAGMENT_FRACTION")
    if stability < config.stability_threshold:
        violations.append("UNSTABLE")
    if config.min_quality is not None and result.quality < config.min_quality:
        violations.append("LOW_QUALITY")
    reason = violations[0] if violations else None
    if reason in {"MIN_CHILD_SIZE", "MAX_CHILD_FRACTION", "TINY_FRAGMENT_FRACTION"}:
        reason = "GAMMA_LIMIT"
    return SplitCandidate(
        gamma,
        result.membership,
        len(sizes),
        result.quality,
        minimum,
        maximum_fraction,
        tiny_fraction,
        stability,
        ari,
        nmi,
        inter_fraction,
        1.0 - inter_fraction,
        reason is None,
        reason,
        violations=tuple(violations),
        structural_valid=len(sizes) >= 2,
        policy_valid=not any(code != "NO_SPLIT" for code in violations),
    )


def _log_grid(lower: float, upper: float, points: int) -> tuple[float, ...]:
    if math.isclose(lower, upper, rel_tol=1e-12, abs_tol=0):
        return (lower,)
    step = (math.log(upper) - math.log(lower)) / (points - 1)
    return tuple(math.exp(math.log(lower) + index * step) for index in range(points))


def search_resolution(
    graph: ig.Graph,
    config: ResolutionSearchConfig,
    *,
    method: str,
    weights: str | list[float] | None,
    seed: int,
    counter: LeidenCallCounter,
    stability_mode: str = "fast",
) -> ResolutionSearchResult:
    config.validate()
    candidates: dict[float, SplitCandidate] = {}

    def evaluate(gamma: float) -> SplitCandidate:
        result, stability, ari, nmi = _multi_seed(
            graph, gamma, method, weights, seed, stability_mode, config.publication_seeds, counter
        )
        return replace(
            _candidate(graph, gamma, result, stability, ari, nmi, config),
            stability_evaluated=stability_mode != "fast",
        )

    if config.strategy == "bounded_adaptive_v2":
        exhausted = False

        def bounded_evaluate(gamma: float, phase: str) -> SplitCandidate | None:
            nonlocal exhausted
            for tested, candidate in candidates.items():
                if math.isclose(gamma, tested, rel_tol=1e-12, abs_tol=0):
                    return candidate
            if len(candidates) >= config.max_candidate_evaluations:
                exhausted = True
                return None
            candidate = replace(
                evaluate(gamma),
                phase=phase,
                evaluation_index=len(candidates) + 1,
                evaluation_budget=config.max_candidate_evaluations,
            )
            candidates[gamma] = candidate
            return candidate

        gamma = config.gamma_min
        for _ in range(config.max_coarse_candidates):
            # Reserve one endpoint slot if it has not been evaluated.
            if gamma > config.gamma_max or (
                len(candidates) >= config.max_candidate_evaluations - 1
                and not math.isclose(gamma, config.gamma_max, rel_tol=1e-12)
            ):
                break
            candidate = bounded_evaluate(gamma, "coarse")
            if candidate is None or candidate.valid:
                break
            gamma *= config.growth_factor
        if not any(item.valid for item in candidates.values()):
            bounded_evaluate(config.gamma_max, "endpoint")
        if not any(item.valid for item in candidates.values()):
            for gamma in _log_grid(
                config.gamma_min, config.gamma_max, config.rescue_grid_points + 2
            )[1:-1]:
                candidate = bounded_evaluate(gamma, "rescue")
                if candidate is None:
                    break
        valid = [item for item in candidates.values() if item.valid]
        if valid:
            lowest = min(item.gamma for item in valid)
            lower = max((value for value in candidates if value < lowest), default=lowest)
            for gamma in _log_grid(lower, lowest, config.local_grid_points):
                if bounded_evaluate(gamma, "local") is None:
                    break
        ordered = tuple(candidates[value] for value in sorted(candidates))
        selected = min(
            (item for item in ordered if item.valid), key=lambda item: item.gamma, default=None
        )
        status = (
            "ACCEPTED"
            if selected
            else "EVALUATION_BUDGET_EXHAUSTED"
            if exhausted
            else "REJECTED_ALL_TESTED"
        )
        return ResolutionSearchResult(selected, ordered, None if selected else status, status)

    if config.strategy == "log_grid":
        for gamma in _log_grid(config.gamma_min, config.gamma_max, config.local_grid_points):
            candidates[gamma] = evaluate(gamma)
    else:
        gamma = config.gamma_min
        previous = gamma
        accepted: SplitCandidate | None = None
        while gamma <= config.gamma_max * (1 + 1e-12):
            candidates[gamma] = evaluate(gamma)
            if candidates[gamma].valid:
                accepted = candidates[gamma]
                break
            previous = gamma
            gamma *= config.growth_factor
        if accepted is not None and accepted.gamma > config.gamma_min:
            for local_gamma in _log_grid(previous, accepted.gamma, config.local_grid_points):
                if not any(math.isclose(local_gamma, value, rel_tol=1e-12) for value in candidates):
                    candidates[local_gamma] = evaluate(local_gamma)
    ordered = tuple(candidates[value] for value in sorted(candidates))
    selected = min(
        (item for item in ordered if item.valid), key=lambda item: item.gamma, default=None
    )
    reasons = {item.rejection_reason for item in ordered}
    terminal_reason = None
    if selected is None:
        if "UNSTABLE" in reasons:
            terminal_reason = "UNSTABLE"
        elif "LOW_QUALITY" in reasons:
            terminal_reason = "LOW_QUALITY"
        elif reasons == {"NO_SPLIT"}:
            terminal_reason = "NO_SPLIT"
        else:
            terminal_reason = "GAMMA_LIMIT"
    return ResolutionSearchResult(
        selected, ordered, terminal_reason, "ACCEPTED" if selected else "REJECTED_ALL_TESTED"
    )
