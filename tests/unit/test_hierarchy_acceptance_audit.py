"""Characterize current policy; these passing tests are not repair acceptance."""

from __future__ import annotations

import igraph as ig
import pytest

from benchmarks.og_extraction.reference_v1 import load_reference
from ogprofiler.graph.components import extract_components
from ogprofiler.graph.edges import EdgeTable
from ogprofiler.hierarchy import resolution as resolution_module
from ogprofiler.hierarchy.engine import HierarchyConfig, infer_component_hierarchy
from ogprofiler.hierarchy.leiden import LeidenCallCounter, LeidenResult
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig, _candidate, search_resolution
from ogprofiler.hierarchy.validation import validate_hierarchy
from ogprofiler.orthogroups.engine import extract_component_orthogroups


def candidate(membership, stability=1.0):
    graph = ig.Graph.Full(len(membership))
    return _candidate(
        graph,
        0.1,
        LeidenResult(tuple(membership), 1.0),
        stability,
        stability,
        1.0,
        ResolutionSearchConfig(
            strategy="adaptive",
            admission_policy="legacy_strict",
        ),
    )


def test_single_singleton_vetoes_an_otherwise_stable_partition():
    result = candidate([0] + [1] * 49 + [2] * 50)
    assert result.child_count == 3
    assert result.stability == 1.0
    assert result.tiny_fragment_fraction == 0.01
    assert result.max_child_fraction == 0.5
    assert result.rejection_reason == "GAMMA_LIMIT" and not result.valid


def test_primary_rejection_hides_simultaneous_instability():
    result = candidate([0] + [1] * 49 + [2] * 50, 0.1)
    assert result.stability < 0.9
    assert result.rejection_reason == "GAMMA_LIMIT"


def search(monkeypatch, fake, config):
    def multiple(graph, gamma, method, weights, seed, mode, publication_seeds, counter):
        membership, stability = fake(gamma)
        counter.count += 1
        return LeidenResult(membership, 1.0), stability, stability, 1.0

    monkeypatch.setattr(resolution_module, "_multi_seed", multiple)
    return search_resolution(
        ig.Graph.Full(6),
        config,
        method="rber",
        weights=None,
        seed=42,
        counter=LeidenCallCounter(),
        stability_mode="robust",
    )


def test_adaptive_search_does_not_test_configured_upper_endpoint(monkeypatch):
    result = search(
        monkeypatch,
        lambda gamma: ((0, 0, 0, 1, 1, 1) if gamma == 10 else (0,) * 6, 1.0),
        ResolutionSearchConfig(
            strategy="adaptive", admission_policy="legacy_strict", gamma_min=0.01, gamma_max=10
        ),
    )
    assert len(result.candidates) == 10
    assert result.candidates[-1].gamma == 5.12
    assert result.selected is None
    # Independent check: gamma=10 would be valid under identical acceptance.
    assert candidate([0, 0, 0, 1, 1, 1]).valid


def test_no_valid_coarse_point_means_no_local_grid_even_with_feasible_interior(monkeypatch):
    result = search(
        monkeypatch,
        lambda gamma: ((0, 0, 0, 1, 1, 1) if 1.3 < gamma < 1.7 else (0,) * 6, 1.0),
        ResolutionSearchConfig(
            strategy="adaptive",
            admission_policy="legacy_strict",
            gamma_min=1,
            gamma_max=2,
            local_grid_points=5,
        ),
    )
    assert [c.gamma for c in result.candidates] == [1, 2]
    assert result.selected is None


def test_terminal_reason_priority_masks_mixed_rejections(monkeypatch):
    result = search(
        monkeypatch,
        lambda gamma: ((0, 0, 0, 1, 1, 1), 0.8) if gamma == 1 else ((0, 1, 1, 1, 1, 1), 1.0),
        ResolutionSearchConfig(
            strategy="adaptive", admission_policy="legacy_strict", gamma_min=1, gamma_max=2
        ),
    )
    assert [c.rejection_reason for c in result.candidates] == ["UNSTABLE", "GAMMA_LIMIT"]
    assert result.terminal_reason == "UNSTABLE"


def test_failed_acceptance_becomes_a_full_residual_og(monkeypatch):
    def multiple(*args):
        return LeidenResult((0, 1, 1, 1, 1, 1), 1.0), 1.0, 1.0, 1.0

    monkeypatch.setattr(resolution_module, "_multi_seed", multiple)
    component = extract_components(
        EdgeTable.canonicalize(list(range(6)), [(p, p + 1, 1.0) for p in range(5)])
    )[0]
    species = {p: p % 2 for p in range(6)}
    result = infer_component_hierarchy(
        component,
        HierarchyConfig(
            resolution=ResolutionSearchConfig(strategy="adaptive", admission_policy="legacy_strict")
        ),
        species,
    )
    validate_hierarchy(component, result)
    assert len(result.nodes) == 1
    assert result.nodes[0].terminal_reason == "GAMMA_LIMIT"
    assert len(result.terminal_membership) == 6
    from dataclasses import asdict

    og = extract_component_orthogroups(
        [asdict(n) for n in result.nodes],
        [dict(protein_id=p, terminal_cluster_id=c) for p, c in result.terminal_membership],
        species,
        {p: f"p{p}" for p in species},
        total_species=2,
    )
    assert len(og.groups) == 1
    assert og.groups[0].selection_type == "RESIDUAL_NONE"
    assert og.groups[0].protein_ids == tuple(range(6))


@pytest.mark.parametrize("size", [2, 3, 4])
def test_actual_v1_accepts_singleton_children_without_global_size_veto(size):
    env = load_reference()
    graph = ig.Graph.Full(size)
    graph.vs["name"] = [f"{p % 2}|p{p}" for p in range(size)]
    graph.vs["index"] = list(range(size))
    obj = env["HHN"](graph)
    # Stub only the expensive community search; execute actual pinned V1 acceptance.
    obj.BipartiteGraphs = lambda *args: (
        [graph.subgraph([0]), graph.subgraph(list(range(1, size)))],
        1.0,
    )
    records = []
    children = obj.RunCommunityDetection(
        records, [("root", list(range(size)))], "rber", "weight", 1.0
    )
    assert len(children) == 2
    assert sorted(len(members) for _, members in children) == [1, size - 1]
    if size < 4:
        component = extract_components(
            EdgeTable.canonicalize(list(range(size)), [(p, p + 1, 1.0) for p in range(size - 1)])
        )[0]
        result = infer_component_hierarchy(
            component,
            HierarchyConfig(
                resolution=ResolutionSearchConfig(
                    strategy="adaptive", admission_policy="legacy_strict"
                )
            ),
            {p: p % 2 for p in range(size)},
        )
        assert result.nodes[0].terminal_reason == "MIN_SIZE"
        assert result.metrics.leiden_calls == 0
    else:
        assert not candidate([0, 1, 1, 1]).valid
