"""Depth-first, component-local hierarchy prototype with bounded graph lifetime."""

from __future__ import annotations

import resource
import sys
import time
from collections.abc import Mapping
from dataclasses import dataclass, replace
from typing import TypedDict

import igraph as ig
import psutil

from ogprofiler.core.models import HierarchyNode
from ogprofiler.graph.components import Component
from ogprofiler.hierarchy.leiden import LeidenCallCounter
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig, _seed_count, search_resolution


@dataclass(frozen=True, slots=True)
class HierarchyConfig:
    recursion_stop_size: int = 1
    component_leiden_call_budget: int | None = None
    leiden_iterations: int = 10
    method: str = "rber"
    seed: int = 42
    max_depth: int = 20
    stability_mode: str = "robust"
    subtree_workers: int = 1
    subtree_release_size: int = 50_000
    resolution: ResolutionSearchConfig = ResolutionSearchConfig()

    def validate(self) -> None:
        self.resolution.validate()
        if type(self.recursion_stop_size) is not int or self.recursion_stop_size < 1:
            raise ValueError("recursion_stop_size must be a positive integer")
        if type(self.leiden_iterations) is not int or self.leiden_iterations < 1:
            raise ValueError("leiden_iterations must be positive and finite")
        if (
            self.component_leiden_call_budget is not None
            and self.resolution.admission_policy == "legacy_strict"
        ):
            raise ValueError("Component budget requires bounded_adaptive_v2")
        if self.component_leiden_call_budget is not None and (
            type(self.component_leiden_call_budget) is not int
            or self.component_leiden_call_budget < 1
        ):
            raise ValueError("component_leiden_call_budget must be a positive integer or null")


@dataclass(frozen=True, slots=True)
class HierarchyMetrics:
    runtime_seconds: float
    peak_rss_bytes: int
    leiden_calls: int
    subgraph_constructions: int
    hierarchy_node_count: int
    terminal_family_count: int
    resolution_candidate_count: int


@dataclass(frozen=True, slots=True)
class ResolutionCandidateTrace:
    cluster_id: int
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


@dataclass(frozen=True, slots=True)
class HierarchyResult:
    component_id: int
    nodes: tuple[HierarchyNode, ...]
    terminal_membership: tuple[tuple[int, int], ...]
    resolution_candidates: tuple[ResolutionCandidateTrace, ...]
    metrics: HierarchyMetrics


@dataclass(frozen=True, slots=True)
class _WorkItem:
    cluster_id: int
    parent_id: int | None
    depth: int
    root_indices: tuple[int, ...]


def _peak_rss_bytes() -> int:
    maximum = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if sys.platform != "darwin":
        maximum *= 1024
    return max(maximum, int(psutil.Process().memory_info().rss))


def _species_count(
    global_ids: tuple[int, ...], species_by_protein: Mapping[int, int] | None
) -> int:
    if species_by_protein is None:
        return 0
    bitmap = 0
    for protein_id in global_ids:
        bitmap |= 1 << species_by_protein[protein_id]
    return bitmap.bit_count()


def _terminal_reason(
    graph: ig.Graph,
    global_ids: tuple[int, ...],
    depth: int,
    config: HierarchyConfig,
    species_by_protein: Mapping[int, int] | None,
) -> str | None:
    if len(global_ids) == 1:
        return "SINGLETON"
    if config.resolution.admission_policy == "nonempty_children_v1":
        if len(global_ids) <= config.recursion_stop_size:
            return "SIZE_STOP"
        if species_by_protein is not None and _species_count(global_ids, species_by_protein) <= 1:
            return "ONE_SPECIES"
        if depth >= config.max_depth:
            return "DEPTH_LIMIT"
        if graph.ecount() == 0:
            return "NO_EDGE_SUPPORT"
        return None
    if len(global_ids) < 2 * config.resolution.min_child_size:
        return "MIN_SIZE"
    if depth >= config.max_depth:
        return "MAX_DEPTH"
    if species_by_protein is not None and _species_count(global_ids, species_by_protein) <= 1:
        return "ONE_SPECIES"
    if graph.ecount() == 0:
        return "NO_EDGES"
    return None


class NodeOutcome(TypedDict, total=False):
    split_status: str
    terminal_reason: str | None
    search_status: str | None
    termination_kind: str | None
    failure_codes: tuple[str, ...]


def node_outcome(reason: str | None, config: HierarchyConfig) -> NodeOutcome:
    """Shared DFS/subtree state semantics; structural leaves retain membership."""
    if config.resolution.admission_policy == "legacy_strict":
        return {"split_status": "TERMINAL" if reason else "SPLIT", "terminal_reason": reason}
    if reason is None:
        return {"split_status": "SPLIT", "search_status": "ACCEPTED"}
    policy = reason in {"SINGLETON", "SIZE_STOP", "ONE_SPECIES"}
    return {
        "split_status": "TERMINAL" if policy else "UNRESOLVED",
        "terminal_reason": reason,
        "search_status": "POLICY_STOP" if policy else reason,
        "termination_kind": "POLICY" if policy else "SEARCH",
        "failure_codes": () if policy else (reason,),
    }


def infer_component_hierarchy(
    component: Component,
    config: HierarchyConfig,
    species_by_protein: Mapping[int, int] | None = None,
    root_graph: ig.Graph | None = None,
) -> HierarchyResult:
    """Infer one component from root to terminal families using explicit DFS."""

    started = time.perf_counter()
    config.validate()
    root_global_ids = component.vertices
    if root_graph is None:
        root_graph, root_global_ids = component.as_edge_table().to_igraph()
    nodes: dict[int, HierarchyNode] = {
        0: HierarchyNode(
            cluster_id=0,
            parent_id=None,
            component_id=component.component_id,
            depth=0,
            n_genes=len(root_global_ids),
            n_species=_species_count(root_global_ids, species_by_protein),
        )
    }
    stack = [_WorkItem(0, None, 0, tuple(range(len(root_global_ids))))]
    terminal_membership: dict[int, int] = {}
    resolution_candidates: list[ResolutionCandidateTrace] = []
    call_counter = LeidenCallCounter(n_iterations=config.leiden_iterations)
    subgraph_constructions = 0
    next_cluster_id = 1

    while stack:
        work = stack.pop()
        if work.depth == 0:
            graph = root_graph
            global_ids = root_global_ids
        else:
            graph = root_graph.induced_subgraph(list(work.root_indices))
            global_ids = tuple(root_global_ids[index] for index in work.root_indices)
            graph.vs["protein_id"] = list(global_ids)
            subgraph_constructions += 1

        reason = _terminal_reason(graph, global_ids, work.depth, config, species_by_protein)
        if reason is not None:
            nodes[work.cluster_id] = replace(nodes[work.cluster_id], **node_outcome(reason, config))
            terminal_membership.update({protein_id: work.cluster_id for protein_id in global_ids})
            continue

        resolution_config = config.resolution
        if config.component_leiden_call_budget is not None:
            available = (config.component_leiden_call_budget - call_counter.count) // _seed_count(
                config.stability_mode, resolution_config.publication_seeds
            )
            if available < 1:
                nodes[work.cluster_id] = replace(
                    nodes[work.cluster_id], **node_outcome("COMPONENT_BUDGET_EXHAUSTED", config)
                )
                terminal_membership.update(
                    {protein_id: work.cluster_id for protein_id in global_ids}
                )
                continue
            resolution_config = replace(
                resolution_config,
                max_candidate_evaluations=min(
                    available, resolution_config.max_candidate_evaluations
                ),
            )
        search = search_resolution(
            graph,
            resolution_config,
            method=config.method,
            weights="weight",
            seed=config.seed,
            counter=call_counter,
            stability_mode=config.stability_mode,
        )
        selected = search.selected
        resolution_candidates.extend(
            ResolutionCandidateTrace(
                cluster_id=work.cluster_id,
                gamma=candidate.gamma,
                child_count=candidate.child_count,
                quality=candidate.quality,
                min_child_size=candidate.min_child_size,
                max_child_fraction=candidate.max_child_fraction,
                tiny_fragment_fraction=candidate.tiny_fragment_fraction,
                stability=candidate.stability,
                adjusted_rand_index=candidate.adjusted_rand_index,
                normalized_mutual_info=candidate.normalized_mutual_info,
                inter_edge_fraction=candidate.inter_edge_fraction,
                intra_edge_fraction=candidate.intra_edge_fraction,
                rejection_reason=candidate.rejection_reason,
                valid=candidate.valid,
                selected=selected is not None and candidate.gamma == selected.gamma,
                violations=candidate.violations,
                structural_valid=candidate.structural_valid,
                policy_valid=candidate.policy_valid,
                phase=candidate.phase,
                evaluation_index=candidate.evaluation_index,
                evaluation_budget=candidate.evaluation_budget,
                stability_evaluated=candidate.stability_evaluated,
            )
            for candidate in search.candidates
        )
        if selected is None:
            nodes[work.cluster_id] = replace(
                nodes[work.cluster_id],
                **node_outcome(search.terminal_reason or "GAMMA_LIMIT", config),
            )
            terminal_membership.update({protein_id: work.cluster_id for protein_id in global_ids})
            continue

        groups: dict[int, list[int]] = {}
        for local_index, community in enumerate(selected.membership):
            groups.setdefault(community, []).append(local_index)
        ordered_groups = sorted(
            (tuple(indices) for indices in groups.values()),
            key=lambda indices: min(global_ids[index] for index in indices),
        )
        nodes[work.cluster_id] = replace(
            nodes[work.cluster_id],
            resolution=selected.gamma,
            quality=selected.quality,
            child_count=len(ordered_groups),
            **node_outcome(None, config),
            selection_phase=selected.phase
            if config.resolution.admission_policy != "legacy_strict"
            else None,
        )

        children: list[_WorkItem] = []
        for indices in ordered_groups:
            child_global_ids = tuple(global_ids[index] for index in indices)
            child_root_indices = tuple(work.root_indices[index] for index in indices)
            child_id = next_cluster_id
            next_cluster_id += 1
            nodes[child_id] = HierarchyNode(
                cluster_id=child_id,
                parent_id=work.cluster_id,
                component_id=component.component_id,
                depth=work.depth + 1,
                n_genes=len(indices),
                n_species=_species_count(child_global_ids, species_by_protein),
            )
            children.append(
                _WorkItem(
                    cluster_id=child_id,
                    parent_id=work.cluster_id,
                    depth=work.depth + 1,
                    root_indices=child_root_indices,
                )
            )
        stack.extend(reversed(children))

    if config.resolution.admission_policy == "nonempty_children_v1":
        # Match sorted subtree paths even when DFS creates grandchildren before siblings run.
        child_ids: dict[int, list[int]] = {}
        for node in nodes.values():
            if node.parent_id is not None:
                child_ids.setdefault(node.parent_id, []).append(node.cluster_id)
        order: list[int] = []
        pending = [0]
        while pending:
            current = pending.pop()
            order.append(current)
            pending.extend(reversed(sorted(child_ids.get(current, []))))
        mapping = {old: new for new, old in enumerate(order)}
        nodes = {
            mapping[old]: replace(
                node,
                cluster_id=mapping[old],
                parent_id=None if node.parent_id is None else mapping[node.parent_id],
            )
            for old, node in nodes.items()
        }
        terminal_membership = {
            protein: mapping[cluster] for protein, cluster in terminal_membership.items()
        }
        resolution_candidates = [
            replace(candidate, cluster_id=mapping[candidate.cluster_id])
            for candidate in resolution_candidates
        ]

    metrics = HierarchyMetrics(
        runtime_seconds=time.perf_counter() - started,
        peak_rss_bytes=_peak_rss_bytes(),
        leiden_calls=call_counter.count,
        subgraph_constructions=subgraph_constructions,
        hierarchy_node_count=len(nodes),
        terminal_family_count=len(set(terminal_membership.values())),
        resolution_candidate_count=len(resolution_candidates),
    )
    return HierarchyResult(
        component_id=component.component_id,
        nodes=tuple(nodes[cluster_id] for cluster_id in sorted(nodes)),
        terminal_membership=tuple(sorted(terminal_membership.items())),
        resolution_candidates=tuple(resolution_candidates),
        metrics=metrics,
    )
