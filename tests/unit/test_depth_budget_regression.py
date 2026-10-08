from benchmarks.og_extraction.depth_budget_regression import compare_prefix, depth_config
from ogprofiler.config import load_config


def test_depth_config_changes_only_max_depth():
    before = load_config(
        overrides=["hierarchy.topology_policy=soft_binary_24_v2", "hierarchy.max_depth=20"]
    )
    after = depth_config(before)
    assert before["hierarchy"]["max_depth"] == 20
    assert after["hierarchy"]["max_depth"] == 42
    after["hierarchy"]["max_depth"] = 20
    assert before == after


def test_prefix_compares_clades_and_traces_not_numeric_ids():
    old = [
        dict(
            cluster_id=0,
            parent_id=None,
            depth=0,
            n_genes=2,
            child_count=2,
            split_status="SPLIT",
            terminal_reason=None,
        ),
        dict(
            cluster_id=1,
            parent_id=0,
            depth=1,
            n_genes=1,
            child_count=0,
            split_status="TERMINAL",
            terminal_reason="SINGLETON",
        ),
        dict(
            cluster_id=2,
            parent_id=0,
            depth=1,
            n_genes=1,
            child_count=0,
            split_status="TERMINAL",
            terminal_reason="SINGLETON",
        ),
    ]
    new = [
        dict(
            n,
            cluster_id=n["cluster_id"] + 10,
            parent_id=None if n["parent_id"] is None else n["parent_id"] + 10,
        )
        for n in old
    ]
    before = [dict(cluster_id=0, gamma=0.5, selected=True)]
    after = [dict(cluster_id=10, gamma=0.5, selected=True)]
    assert compare_prefix(old, [(0, 1), (1, 2)], before, new, [(0, 11), (1, 12)], after)["passed"]
    after[0]["gamma"] = 0.6
    assert not compare_prefix(old, [(0, 1), (1, 2)], before, new, [(0, 11), (1, 12)], after)[
        "passed"
    ]
