"""Recursive component-subtree scheduler with schedule-independent topology IDs."""

from __future__ import annotations

import multiprocessing
import time
from collections.abc import Mapping
from concurrent.futures import FIRST_COMPLETED, Future, ProcessPoolExecutor, wait
from dataclasses import asdict, dataclass

import igraph as ig

from ogprofiler.core.models import HierarchyNode
from ogprofiler.graph.components import Component
from ogprofiler.hierarchy.engine import (
    HierarchyConfig,
    HierarchyMetrics,
    HierarchyResult,
    ResolutionCandidateTrace,
    _peak_rss_bytes,
    _species_count,
    _terminal_reason,
    infer_component_hierarchy,
    node_outcome,
)
from ogprofiler.hierarchy.leiden import LeidenCallCounter
from ogprofiler.hierarchy.resolution import search_resolution


@dataclass(frozen=True, order=True, slots=True)
class SubtreePath:
    component_id: int
    child_ordinals: tuple[int, ...] = ()

    @property
    def task_id(self) -> str:
        suffix = (
            "root"
            if not self.child_ordinals
            else ".".join(f"{ordinal:06d}" for ordinal in self.child_ordinals)
        )
        return f"component-{self.component_id:08d}/{suffix}"

    @property
    def parent(self) -> SubtreePath | None:
        if not self.child_ordinals:
            return None
        return SubtreePath(self.component_id, self.child_ordinals[:-1])

    def child(self, ordinal: int) -> SubtreePath:
        if ordinal < 0:
            raise ValueError("Child ordinal must be non-negative")
        return SubtreePath(self.component_id, (*self.child_ordinals, ordinal))


@dataclass(frozen=True, slots=True)
class _Task:
    path: SubtreePath
    root_indices: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class _Candidate:
    gamma: float
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
    rejection_reason: str | None
    valid: bool
    selected: bool
    violations: tuple[str, ...] = ()
    structural_valid: bool = True
    policy_valid: bool = True
    phase: str = "legacy"
    evaluation_index: int = 0
    evaluation_budget: int = 0
    stability_evaluated: bool = True
    binary_eligible: bool | None = None
    kway_eligible: bool | None = None
    original_violations: tuple[str, ...] = ()
    selection_kind: str | None = None
    refinement_truncated: bool = False


@dataclass(frozen=True, slots=True)
class _Outcome:
    task: _Task
    terminal_reason: str | None
    resolution: float | None
    quality: float | None
    children: tuple[_Task, ...]
    candidates: tuple[_Candidate, ...]
    leiden_calls: int
    peak_rss_bytes: int


_GRAPH: ig.Graph | None = None
_GLOBAL_IDS: tuple[int, ...] = ()
_SPECIES: Mapping[int, int] | None = None
_CONFIG: HierarchyConfig | None = None


def _initialize(
    graph: ig.Graph,
    global_ids: tuple[int, ...],
    species: Mapping[int, int] | None,
    config: HierarchyConfig,
) -> None:
    global _GRAPH, _GLOBAL_IDS, _SPECIES, _CONFIG
    _GRAPH, _GLOBAL_IDS, _SPECIES, _CONFIG = graph, global_ids, species, config


def _evaluate(task: _Task) -> _Outcome:
    if _GRAPH is None or _CONFIG is None:
        raise RuntimeError("Subtree worker was not initialized")
    depth = len(task.path.child_ordinals)
    graph = _GRAPH if depth == 0 else _GRAPH.induced_subgraph(list(task.root_indices))
    global_ids = tuple(_GLOBAL_IDS[index] for index in task.root_indices)
    graph.vs["protein_id"] = list(global_ids)
    reason = _terminal_reason(graph, global_ids, depth, _CONFIG, _SPECIES)
    if reason is not None:
        return _Outcome(task, reason, None, None, (), (), 0, _peak_rss_bytes())

    counter = LeidenCallCounter(n_iterations=_CONFIG.leiden_iterations)
    search = search_resolution(
        graph,
        _CONFIG.resolution,
        method=_CONFIG.method,
        weights="weight",
        seed=_CONFIG.seed,
        counter=counter,
        stability_mode=_CONFIG.stability_mode,
    )
    selected = search.selected
    candidates = tuple(
        _Candidate(
            candidate.gamma,
            candidate.child_count,
            candidate.quality,
            candidate.min_child_size,
            candidate.max_child_fraction,
            candidate.tiny_fragment_fraction,
            candidate.stability,
            candidate.adjusted_rand_index,
            candidate.normalized_mutual_info,
            candidate.inter_edge_fraction,
            candidate.intra_edge_fraction,
            candidate.rejection_reason,
            candidate.valid,
            selected is not None and candidate.gamma == selected.gamma,
            candidate.violations,
            candidate.structural_valid,
            candidate.policy_valid,
            candidate.phase,
            candidate.evaluation_index,
            candidate.evaluation_budget,
            candidate.stability_evaluated,
            candidate.binary_eligible,
            candidate.kway_eligible,
            candidate.original_violations,
            candidate.selection_kind,
            candidate.refinement_truncated,
        )
        for candidate in search.candidates
    )
    if selected is None:
        return _Outcome(
            task,
            search.terminal_reason or "GAMMA_LIMIT",
            None,
            None,
            (),
            candidates,
            counter.count,
            _peak_rss_bytes(),
        )

    groups: dict[int, list[int]] = {}
    for local_index, community in enumerate(selected.membership):
        groups.setdefault(community, []).append(local_index)
    ordered = sorted(
        (tuple(indices) for indices in groups.values()),
        key=lambda indices: min(global_ids[index] for index in indices),
    )
    children = tuple(
        _Task(
            task.path.child(ordinal),
            tuple(task.root_indices[index] for index in indices),
        )
        for ordinal, indices in enumerate(ordered)
    )
    return _Outcome(
        task,
        None,
        selected.gamma,
        selected.quality,
        children,
        candidates,
        counter.count,
        _peak_rss_bytes(),
    )


def infer_component_hierarchy_parallel(
    component: Component,
    config: HierarchyConfig,
    species_by_protein: Mapping[int, int] | None,
    root_graph: ig.Graph,
    *,
    workers: int,
) -> HierarchyResult:
    """Release every accepted split's children as independent process tasks."""
    if workers < 2:
        raise ValueError("Parallel subtree scheduling requires at least two workers")
    # A shared component cap uses deterministic DFS reservation, not per-worker caps.
    if config.component_leiden_call_budget is not None:
        return infer_component_hierarchy(component, config, species_by_protein, root_graph)
    config.validate()
    started = time.perf_counter()
    root_ids = component.vertices
    root = _Task(SubtreePath(component.component_id), tuple(range(len(root_ids))))
    outcomes: dict[SubtreePath, _Outcome] = {}
    context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=context,
        initializer=_initialize,
        initargs=(root_graph, root_ids, species_by_protein, config),
    ) as executor:
        futures: dict[Future[_Outcome], _Task] = {executor.submit(_evaluate, root): root}
        while futures:
            done, _ = wait(futures, return_when=FIRST_COMPLETED)
            for future in done:
                task = futures.pop(future)
                outcome = future.result()
                outcomes[task.path] = outcome
                for child in outcome.children:
                    futures[executor.submit(_evaluate, child)] = child

    ordered_paths = sorted(outcomes)
    cluster_by_path = {path: cluster_id for cluster_id, path in enumerate(ordered_paths)}
    nodes: list[HierarchyNode] = []
    membership: list[tuple[int, int]] = []
    traces: list[ResolutionCandidateTrace] = []
    for path in ordered_paths:
        outcome = outcomes[path]
        cluster_id = cluster_by_path[path]
        global_ids = tuple(root_ids[index] for index in outcome.task.root_indices)
        nodes.append(
            HierarchyNode(
                cluster_id=cluster_id,
                parent_id=None if path.parent is None else cluster_by_path[path.parent],
                component_id=component.component_id,
                depth=len(path.child_ordinals),
                n_genes=len(global_ids),
                n_species=_species_count(global_ids, species_by_protein),
                resolution=outcome.resolution,
                quality=outcome.quality,
                child_count=len(outcome.children),
                selection_kind=next(
                    (c.selection_kind for c in outcome.candidates if c.selected), None
                ),
                refinement_truncated=any(
                    c.refinement_truncated for c in outcome.candidates if c.selected
                ),
                **node_outcome(outcome.terminal_reason, config),
                selection_phase=next(
                    (item.phase for item in outcome.candidates if item.selected), None
                )
                if config.resolution.admission_policy != "legacy_strict"
                else None,
            )
        )
        if outcome.terminal_reason:
            membership.extend((protein_id, cluster_id) for protein_id in global_ids)
        traces.extend(
            ResolutionCandidateTrace(cluster_id=cluster_id, **asdict(candidate))
            for candidate in outcome.candidates
        )
    return HierarchyResult(
        component_id=component.component_id,
        nodes=tuple(nodes),
        terminal_membership=tuple(sorted(membership)),
        resolution_candidates=tuple(traces),
        metrics=HierarchyMetrics(
            runtime_seconds=time.perf_counter() - started,
            peak_rss_bytes=max(outcome.peak_rss_bytes for outcome in outcomes.values()),
            leiden_calls=sum(outcome.leiden_calls for outcome in outcomes.values()),
            subgraph_constructions=max(0, len(outcomes) - 1),
            hierarchy_node_count=len(outcomes),
            terminal_family_count=sum(
                bool(outcome.terminal_reason) for outcome in outcomes.values()
            ),
            resolution_candidate_count=sum(
                len(outcome.candidates) for outcome in outcomes.values()
            ),
        ),
    )
