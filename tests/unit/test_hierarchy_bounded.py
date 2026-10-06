"""ADR 0002 public search contracts."""

import igraph as ig
import pytest

from ogprofiler.hierarchy.leiden import LeidenCallCounter, LeidenResult
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig, search_resolution


def search(monkeypatch, partition, **kwargs):
    calls = []

    def run(graph, gamma, method, weights, seed, counter):
        calls.append(gamma)
        counter.count += 1
        return LeidenResult(partition(gamma), 1.0)

    monkeypatch.setattr("ogprofiler.hierarchy.resolution.run_leiden", run)
    result = search_resolution(
        ig.Graph.Full(6),
        ResolutionSearchConfig(
            strategy="bounded_adaptive_v2", admission_policy="nonempty_children_v1", **kwargs
        ),
        method="rber",
        weights=None,
        seed=42,
        counter=LeidenCallCounter(),
        stability_mode="robust",
    )
    return result, calls


def test_singleton_admission(monkeypatch):
    result, calls = search(monkeypatch, lambda gamma: (0, 1, 1, 1, 1, 1))
    assert result.selected is not None
    assert result.selected.min_child_size == 1
    assert result.search_status == "ACCEPTED"
    assert len(calls) <= 72


def test_endpoint_and_rescue(monkeypatch):
    result, calls = search(
        monkeypatch, lambda gamma: (0, 0, 0, 1, 1, 1) if gamma == 10 else (0,) * 6
    )
    assert result.selected.gamma == 10
    assert result.selected.phase == "endpoint"
    assert len(set(calls)) <= 24


def test_budget_failure_is_not_terminal_success(monkeypatch):
    result, calls = search(monkeypatch, lambda gamma: (0,) * 6, max_candidate_evaluations=3)
    assert result.selected is None
    assert result.search_status == "EVALUATION_BUDGET_EXHAUSTED"
    assert all(candidate.evaluation_budget == 3 for candidate in result.candidates)
    assert len(calls) == 9
    assert 10 in calls


def test_all_violations(monkeypatch):
    result, _ = search(
        monkeypatch,
        lambda gamma: (0, 1, 1, 1, 1, 1),
        max_child_fraction=0.5,
        max_tiny_fragment_fraction=0,
    )
    assert {"MAX_CHILD_FRACTION", "TINY_FRAGMENT_FRACTION"} <= set(result.candidates[0].violations)
    assert result.search_status == "REJECTED_ALL_TESTED"


def test_depth_limit_preserves_members_as_unresolved():
    from ogprofiler.graph.components import Component
    from ogprofiler.hierarchy.engine import HierarchyConfig, infer_component_hierarchy

    graph = ig.Graph.Full(6)
    graph.es["weight"] = [1.0] * graph.ecount()
    component = Component(0, tuple(range(6)), tuple())
    result = infer_component_hierarchy(
        component,
        HierarchyConfig(
            max_depth=0,
            resolution=ResolutionSearchConfig(
                strategy="bounded_adaptive_v2", admission_policy="nonempty_children_v1"
            ),
        ),
        {i: i % 2 for i in range(6)},
        graph,
    )
    assert result.nodes[0].split_status == "UNRESOLVED"
    assert result.nodes[0].search_status == "DEPTH_LIMIT"
    assert len(result.terminal_membership) == 6


def test_rescue_nonmonotonic_and_fast_tag(monkeypatch):
    import math

    rescue = 0.01 * (10 / 0.01) ** (1 / 9)
    result, calls = search(
        monkeypatch,
        lambda gamma: (
            (0, 0, 0, 1, 1, 1) if math.isclose(gamma, rescue, rel_tol=1e-12) else (0,) * 6
        ),
    )
    assert result.selected.phase == "rescue"
    assert result.selected.gamma == pytest.approx(rescue)
    assert len(set(calls)) <= 24


def test_component_budget_and_size_stop():
    from ogprofiler.graph.components import Component
    from ogprofiler.hierarchy.engine import HierarchyConfig, infer_component_hierarchy

    graph = ig.Graph.Full(6)
    graph.es["weight"] = [1.0] * graph.ecount()
    component = Component(0, tuple(range(6)), tuple())
    config = HierarchyConfig(
        component_leiden_call_budget=2,
        stability_mode="robust",
        resolution=ResolutionSearchConfig(
            strategy="bounded_adaptive_v2", admission_policy="nonempty_children_v1"
        ),
    )
    result = infer_component_hierarchy(component, config, {i: i % 2 for i in range(6)}, graph)
    assert result.metrics.leiden_calls == 0
    assert result.nodes[0].failure_codes == ("COMPONENT_BUDGET_EXHAUSTED",)
    from dataclasses import replace

    result = infer_component_hierarchy(
        component, replace(config, recursion_stop_size=6), {i: i % 2 for i in range(6)}, graph
    )
    assert result.nodes[0].terminal_reason == "SIZE_STOP"
    assert result.nodes[0].split_status == "TERMINAL"


def test_new_config_defaults_and_explicit_migration():
    from ogprofiler.config import hierarchy_config, load_config
    from ogprofiler.exceptions import InputError

    config = hierarchy_config(load_config()["hierarchy"])
    assert config.resolution.strategy == "bounded_adaptive_v2"
    assert config.leiden_iterations == 10
    assert config.recursion_stop_size == 1
    assert config.resolution.stability_threshold == 0.9
    assert config.resolution.max_candidate_evaluations == 24
    with pytest.raises(InputError, match="legacy-only"):
        load_config(overrides=["hierarchy.min_family_size=3"])
    legacy = load_config(
        overrides=[
            "hierarchy.admission_policy=legacy_strict",
            "hierarchy.topology_policy=kway_v1",
            "hierarchy.resolution_strategy=adaptive",
            "hierarchy.min_family_size=3",
        ]
    )
    assert hierarchy_config(legacy["hierarchy"]).resolution.min_child_size == 3
    assert hierarchy_config(legacy["hierarchy"]).leiden_iterations == 2


def test_saved_iteration_budget_and_cli_override_are_preserved(tmp_path):
    from ogprofiler.config import hierarchy_config, load_config

    path = tmp_path / "run.yaml"
    path.write_text("hierarchy:\n  leiden_iterations: 2\n")
    assert hierarchy_config(load_config(str(path))["hierarchy"]).leiden_iterations == 2
    config = load_config(str(path), ["hierarchy.leiden_iterations=10"])
    assert hierarchy_config(config["hierarchy"]).leiden_iterations == 10


def test_legacy_yaml_omitted_budget_and_explicit_override(tmp_path):
    from ogprofiler.config import hierarchy_config, load_config

    path = tmp_path / "legacy.yaml"
    path.write_text(
        "hierarchy:\n  admission_policy: legacy_strict\n"
        "  topology_policy: kway_v1\n  resolution_strategy: adaptive\n"
    )
    assert hierarchy_config(load_config(str(path))["hierarchy"]).leiden_iterations == 2
    config = load_config(str(path), ["hierarchy.leiden_iterations=10"])
    assert hierarchy_config(config["hierarchy"]).leiden_iterations == 10
    config = load_config(
        str(path),
        [
            "hierarchy.admission_policy=nonempty_children_v1",
            "hierarchy.resolution_strategy=bounded_adaptive_v2",
        ],
    )
    assert hierarchy_config(config["hierarchy"]).leiden_iterations == 10


def test_new_serial_parallel_equivalence():
    from test_hierarchy import planted_component

    from ogprofiler.hierarchy.engine import HierarchyConfig, infer_component_hierarchy
    from ogprofiler.hierarchy.subtree import infer_component_hierarchy_parallel
    from ogprofiler.hierarchy.validation import validate_hierarchy

    component = planted_component()
    graph, _ = component.as_edge_table().to_igraph()
    config = HierarchyConfig(
        max_depth=3,
        resolution=ResolutionSearchConfig(gamma_min=0.1, gamma_max=2, max_child_fraction=0.8),
    )
    serial = infer_component_hierarchy(component, config, root_graph=graph)
    parallel = infer_component_hierarchy_parallel(component, config, None, graph, workers=2)
    validate_hierarchy(component, serial)
    validate_hierarchy(component, parallel)
    assert serial.nodes == parallel.nodes
    assert serial.terminal_membership == parallel.terminal_membership
    assert serial.resolution_candidates == parallel.resolution_candidates
    assert serial.metrics.leiden_calls == parallel.metrics.leiden_calls


def test_invalid_membership_is_computation_error(monkeypatch):
    from ogprofiler.exceptions import HierarchyError

    with pytest.raises(HierarchyError, match="COMPUTATION_ERROR"):
        search(monkeypatch, lambda gamma: (0, 1))


def test_singleton_real_child_and_fast_measurement_tag(monkeypatch):
    from ogprofiler.graph.components import Component
    from ogprofiler.hierarchy.engine import HierarchyConfig, infer_component_hierarchy

    graph = ig.Graph.Full(6)
    graph.es["weight"] = [1.0] * graph.ecount()

    def run(graph, gamma, method, weights, seed, counter):
        assert counter.n_iterations == 10
        counter.count += 1
        return LeidenResult((0, 1, 1, 1, 1, 1), 1.0)

    monkeypatch.setattr("ogprofiler.hierarchy.resolution.run_leiden", run)
    result = infer_component_hierarchy(
        Component(0, tuple(range(6)), tuple()),
        HierarchyConfig(stability_mode="fast"),
        {0: 0, **{i: 1 for i in range(1, 6)}},
        graph,
    )
    assert result.nodes[0].split_status == "SPLIT"
    assert result.nodes[1].terminal_reason == "SINGLETON"
    assert result.terminal_membership[0] == (0, 1)
    assert all(not item.stability_evaluated for item in result.resolution_candidates)


def test_new_search_keeps_stability_threshold(monkeypatch):
    def run(graph, gamma, method, weights, seed, counter):
        counter.count += 1
        membership = (0, 0, 0, 1, 1, 1) if seed % 2 == 0 else (0, 1, 0, 1, 0, 1)
        return LeidenResult(membership, 1.0)

    monkeypatch.setattr("ogprofiler.hierarchy.resolution.run_leiden", run)
    counter = LeidenCallCounter()
    result = search_resolution(
        ig.Graph.Full(6),
        ResolutionSearchConfig(),
        method="rber",
        weights=None,
        seed=42,
        counter=counter,
        stability_mode="robust",
    )
    assert result.selected is None
    assert result.search_status == "REJECTED_ALL_TESTED"
    assert all("UNSTABLE" in item.violations for item in result.candidates)
    assert counter.count <= 72
