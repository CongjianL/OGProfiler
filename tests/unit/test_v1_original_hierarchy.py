"""Native V1 search/stop behavior differs intentionally from V2 gates."""

import igraph

from benchmarks.og_extraction.v1_original_hierarchy import original_tree, paths_and_clades


class Partition:
    def __init__(self, graph, membership):
        self.graph = graph
        self.membership = membership
        self.q = 1.0

    def subgraphs(self):
        return [
            self.graph.induced_subgraph([i for i, c in enumerate(self.membership) if c == label])
            for label in sorted(set(self.membership))
        ]


def test_v1_two_gene_special_case_has_no_leiden_call():
    graph = igraph.Graph.Full(2)
    nodes, calls = original_tree(graph, (8, 19), {8: 0, 19: 1})
    assert calls == []
    assert nodes[0]["children"] == ["0-0", "0-1"]
    paths, _ = paths_and_clades(nodes)
    assert paths == {8: ["0", "0-0"], 19: ["0", "0-1"]}


def test_original_seed_free_binary_preserves_duplicated_first_evaluation():
    graph = igraph.Graph.Full(3)
    kwargs_seen = []

    def partition(g, *args, **kwargs):
        kwargs_seen.append(kwargs)
        return Partition(g, [0, 1, 1])

    nodes, calls = original_tree(
        graph, (10, 11, 12), {10: 0, 11: 1, 12: 1}, find_partition=partition
    )
    assert len(calls) == 2
    assert [c["gamma"] for c in calls] == [0.5, 0.5]
    assert all(c["n_iterations"] == 10 and "seed" not in c for c in kwargs_seen)
    assert nodes[0]["children"] == ["0-0", "0-1"]
    assert nodes[-1]["terminal_reason"] == "ONE_SPECIES"
    assert set(paths_and_clades(nodes)[0]) == {10, 11, 12}


def test_original_failure_is_terminal_after_1001_calls_not_v2_unresolved():
    graph = igraph.Graph.Full(3)

    def partition(g, *args, **kwargs):
        return Partition(g, [0, 1, 2])

    nodes, calls = original_tree(
        graph, (10, 11, 12), {10: 0, 11: 1, 12: 2}, find_partition=partition
    )
    assert len(calls) == 1001
    assert len(nodes) == 1
    assert nodes[0]["terminal_reason"] == "V1_SEARCH_NOT_EXACTLY_TWO"
    assert nodes[0]["children"] == []
