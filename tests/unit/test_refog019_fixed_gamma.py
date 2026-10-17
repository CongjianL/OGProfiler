"""Fixed point must use real V2 multi-seed stability and preserve input config."""

from dataclasses import asdict

import igraph
import pytest

from benchmarks.og_extraction.refog019_fixed_gamma import evaluate_fixed
from ogprofiler.hierarchy import resolution
from ogprofiler.hierarchy.engine import HierarchyConfig
from ogprofiler.hierarchy.leiden import LeidenResult
from ogprofiler.hierarchy.resolution import ResolutionSearchConfig


@pytest.mark.parametrize("unstable", [True, False])
def test_fixed_gamma_keeps_original_gates_and_reports_seed_partitions(monkeypatch, unstable):
    config = HierarchyConfig(
        seed=42,
        stability_mode="robust",
        leiden_iterations=10,
        resolution=ResolutionSearchConfig(topology_policy="soft_binary_24_v2"),
    )
    before = asdict(config)
    graph = igraph.Graph.Full(12)
    graph.es["weight"] = [1.0] * graph.ecount()

    def fake(g, gamma, method, weights, seed, counter):
        assert gamma == 0.953125 and method == "rber" and weights == "weight"
        assert counter.n_iterations == 10
        counter.count += 1
        singleton = 1 if unstable and seed == 209500 else 0
        membership = tuple(1 if p == singleton else 0 for p in range(12))
        return LeidenResult(membership, 1.0)

    monkeypatch.setattr(resolution, "run_leiden", fake)
    candidate, traces, pairs, calls = evaluate_fixed(graph, config)
    assert asdict(config) == before
    assert calls == 3 and [t["seed"] for t in traces] == [42, 104771, 209500]
    assert candidate.child_count == 2 and candidate.max_child_fraction == 11 / 12
    if unstable:
        assert candidate.stability == pytest.approx(3 / 11)
        assert candidate.original_violations == ("UNSTABLE",)
        assert not candidate.binary_eligible
        assert [p["ari"] for p in pairs] == pytest.approx([1, -1 / 11, -1 / 11])
    else:
        assert candidate.stability == 1 and candidate.binary_eligible
        assert candidate.original_violations == ()
