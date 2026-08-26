from __future__ import annotations

from pathlib import Path

import igraph as ig

from ogprofiler.graph.components import Component, extract_components
from ogprofiler.graph.edges import EdgeTable
from ogprofiler.hierarchy import resolution as resolution_module
from ogprofiler.hierarchy.engine import (
    HierarchyConfig,
    infer_component_hierarchy,
)
from ogprofiler.hierarchy.leiden import LeidenCallCounter, LeidenResult, run_leiden
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig, search_resolution
from ogprofiler.hierarchy.subtree import infer_component_hierarchy_parallel
from ogprofiler.hierarchy.validation import validate_hierarchy
from ogprofiler.storage.hierarchy import write_hierarchy_result


def planted_component() -> Component:
    edges: list[tuple[int, int, float]] = []
    for group in range(3):
        vertices = range(group * 5, (group + 1) * 5)
        for left in vertices:
            for right in vertices:
                if left < right:
                    edges.append((left, right, 10.0))
    edges.extend([(4, 5, 0.01), (9, 10, 0.01)])
    table = EdgeTable.canonicalize(list(range(15)), edges)
    return extract_components(table)[0]


def test_leiden_wrapper_is_seed_reproducible_and_returns_k_way_split() -> None:
    component = planted_component()
    graph, _ = component.as_edge_table().to_igraph()
    first = run_leiden(graph, 0.1, "rber", "weight", 42)
    second = run_leiden(graph, 0.1, "rber", "weight", 42)
    assert first.membership == second.membership
    assert len(set(first.membership)) == 3


def test_resolution_search_refines_to_lowest_valid_candidate(monkeypatch: object) -> None:
    graph = ig.Graph.Full(6)

    def fake_run(
        graph: ig.Graph,
        gamma: float,
        method: str,
        weights: str | list[float] | None,
        seed: int,
        counter: LeidenCallCounter,
    ) -> LeidenResult:
        del graph, method, weights, seed
        counter.count += 1
        membership = (0, 0, 0, 1, 1, 1) if gamma >= 0.3 else (0, 0, 0, 0, 0, 0)
        return LeidenResult(membership, gamma)

    monkeypatch.setattr(resolution_module, "run_leiden", fake_run)  # type: ignore[attr-defined]
    counter = LeidenCallCounter()
    result = search_resolution(
        graph,
        ResolutionSearchConfig(
            gamma_min=0.1,
            gamma_max=1.0,
            growth_factor=2.0,
            local_grid_points=5,
            min_child_size=2,
            max_child_fraction=0.8,
        ),
        method="rber",
        weights=None,
        seed=42,
        counter=counter,
    )
    assert result.selected is not None
    assert 0.3 <= result.selected.gamma < 0.4
    assert result.selected.child_count == 2
    assert counter.count == len(result.candidates)


def test_dfs_hierarchy_preserves_all_invariants_and_k_way_split(tmp_path: Path) -> None:
    component = planted_component()
    config = HierarchyConfig(
        max_depth=3,
        resolution=ResolutionSearchConfig(
            gamma_min=0.1,
            gamma_max=2.0,
            growth_factor=2.0,
            local_grid_points=4,
            min_child_size=2,
            max_child_fraction=0.8,
        ),
    )
    result = infer_component_hierarchy(component, config)
    validate_hierarchy(component, result)
    root_children = [node for node in result.nodes if node.parent_id == 0]
    assert len(root_children) == 3
    assert sorted(node.n_genes for node in root_children) == [5, 5, 5]
    assert len(result.terminal_membership) == 15
    assert result.metrics.terminal_family_count == 3
    assert result.metrics.subgraph_constructions == 3
    assert result.metrics.leiden_calls > 0
    assert result.metrics.resolution_candidate_count == result.metrics.leiden_calls
    assert len(result.resolution_candidates) == result.metrics.leiden_calls
    assert any(candidate.selected for candidate in result.resolution_candidates)
    assert all(not hasattr(node, "gene_ids") for node in result.nodes)

    write_hierarchy_result(tmp_path, result)
    assert (tmp_path / "nodes.parquet").is_file()
    assert (tmp_path / "members.parquet").is_file()
    assert (tmp_path / "candidates.parquet").is_file()
    assert (tmp_path / "metrics.json").is_file()


def test_parallel_subtree_scheduler_matches_frozen_dfs_topology() -> None:
    component = planted_component()
    graph, _ = component.as_edge_table().to_igraph()
    config = HierarchyConfig(
        max_depth=3,
        resolution=ResolutionSearchConfig(
            gamma_min=0.1,
            gamma_max=2.0,
            growth_factor=2.0,
            local_grid_points=4,
            min_child_size=2,
            max_child_fraction=0.8,
        ),
    )
    serial = infer_component_hierarchy(component, config, root_graph=graph)
    parallel = infer_component_hierarchy_parallel(
        component, config, None, graph, workers=2
    )
    validate_hierarchy(component, parallel)
    assert parallel.nodes == serial.nodes
    assert parallel.terminal_membership == serial.terminal_membership
    assert parallel.resolution_candidates == serial.resolution_candidates


def test_singleton_component_bypasses_leiden() -> None:
    component = extract_components(EdgeTable.canonicalize([7], []))[0]
    result = infer_component_hierarchy(component, HierarchyConfig())
    validate_hierarchy(component, result)
    assert result.nodes[0].terminal_reason == "SINGLETON"
    assert result.metrics.leiden_calls == 0
    assert result.terminal_membership == ((7, 0),)


def test_species_and_depth_terminal_rules_bypass_extra_search() -> None:
    component = planted_component()
    one_species = infer_component_hierarchy(
        component,
        HierarchyConfig(),
        species_by_protein={protein_id: 0 for protein_id in component.vertices},
    )
    assert one_species.nodes[0].terminal_reason == "ONE_SPECIES"
    assert one_species.metrics.leiden_calls == 0

    depth_limited = infer_component_hierarchy(
        component,
        HierarchyConfig(
            max_depth=1,
            resolution=ResolutionSearchConfig(
                gamma_min=0.1,
                gamma_max=1.0,
                min_child_size=2,
                max_child_fraction=0.8,
            ),
        ),
    )
    children = [node for node in depth_limited.nodes if node.parent_id == 0]
    assert len(children) == 3
    assert {node.terminal_reason for node in children} == {"MAX_DEPTH"}


def test_robust_search_classifies_unstable_and_low_quality(monkeypatch: object) -> None:
    graph = ig.Graph.Full(6)

    def unstable_run(
        graph: ig.Graph,
        gamma: float,
        method: str,
        weights: str | list[float] | None,
        seed: int,
        counter: LeidenCallCounter,
    ) -> LeidenResult:
        del graph, gamma, method, weights
        counter.count += 1
        membership = (0, 0, 0, 1, 1, 1) if seed % 2 == 0 else (0, 1, 0, 1, 0, 1)
        return LeidenResult(membership, 1.0)

    monkeypatch.setattr(resolution_module, "run_leiden", unstable_run)  # type: ignore[attr-defined]
    robust_counter = LeidenCallCounter()
    unstable = search_resolution(
        graph,
        ResolutionSearchConfig(
            gamma_min=0.1,
            gamma_max=0.1,
            min_child_size=2,
            max_child_fraction=0.8,
            stability_threshold=0.9,
        ),
        method="rber",
        weights=None,
        seed=42,
        counter=robust_counter,
        stability_mode="robust",
    )
    assert unstable.selected is None
    assert unstable.terminal_reason == "UNSTABLE"
    assert unstable.candidates[0].adjusted_rand_index < 0.9
    assert robust_counter.count == 3

    def low_quality_run(
        graph: ig.Graph,
        gamma: float,
        method: str,
        weights: str | list[float] | None,
        seed: int,
        counter: LeidenCallCounter,
    ) -> LeidenResult:
        del graph, gamma, method, weights, seed
        counter.count += 1
        return LeidenResult((0, 0, 0, 1, 1, 1), 0.1)

    monkeypatch.setattr(resolution_module, "run_leiden", low_quality_run)  # type: ignore[attr-defined]
    low_quality = search_resolution(
        graph,
        ResolutionSearchConfig(
            gamma_min=0.1,
            gamma_max=0.1,
            min_child_size=2,
            max_child_fraction=0.8,
            min_quality=0.5,
        ),
        method="rber",
        weights=None,
        seed=42,
        counter=LeidenCallCounter(),
    )
    assert low_quality.selected is None
    assert low_quality.terminal_reason == "LOW_QUALITY"


def test_no_split_and_constraint_failure_have_distinct_terminal_reasons(
    monkeypatch: object,
) -> None:
    graph = ig.Graph.Full(6)

    def fake_run(
        graph: ig.Graph,
        gamma: float,
        method: str,
        weights: str | list[float] | None,
        seed: int,
        counter: LeidenCallCounter,
    ) -> LeidenResult:
        del graph, method, weights, seed
        counter.count += 1
        membership = (0, 0, 0, 0, 0, 0) if gamma < 0.2 else (0, 0, 0, 0, 0, 1)
        return LeidenResult(membership, 1.0)

    monkeypatch.setattr(resolution_module, "run_leiden", fake_run)  # type: ignore[attr-defined]
    no_split = search_resolution(
        graph,
        ResolutionSearchConfig(gamma_min=0.1, gamma_max=0.1),
        method="rber",
        weights=None,
        seed=1,
        counter=LeidenCallCounter(),
    )
    constrained = search_resolution(
        graph,
        ResolutionSearchConfig(gamma_min=0.2, gamma_max=0.2, min_child_size=2),
        method="rber",
        weights=None,
        seed=1,
        counter=LeidenCallCounter(),
    )
    assert no_split.terminal_reason == "NO_SPLIT"
    assert constrained.terminal_reason == "GAMMA_LIMIT"
