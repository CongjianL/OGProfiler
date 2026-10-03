"""Public frozen-tree depth audit uses protein sets, not cluster ID equality."""

from benchmarks.og_extraction.depth_path_audit import audit_trees, shrink_steps


def test_depth_audit_distinguishes_extra_binary_edges_from_membership_change():
    old = [
        dict(
            cluster_id=0,
            parent_id=None,
            depth=0,
            n_genes=4,
            child_count=3,
            split_status="SPLIT",
            terminal_reason=None,
        ),
        dict(
            cluster_id=1,
            parent_id=0,
            depth=1,
            n_genes=2,
            child_count=2,
            split_status="SPLIT",
            terminal_reason=None,
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
        dict(
            cluster_id=3,
            parent_id=0,
            depth=1,
            n_genes=1,
            child_count=0,
            split_status="TERMINAL",
            terminal_reason="SINGLETON",
        ),
    ]
    old.extend(
        [
            dict(
                cluster_id=4,
                parent_id=1,
                depth=2,
                n_genes=1,
                child_count=0,
                split_status="TERMINAL",
                terminal_reason="SINGLETON",
            ),
            dict(
                cluster_id=5,
                parent_id=1,
                depth=2,
                n_genes=1,
                child_count=0,
                split_status="TERMINAL",
                terminal_reason="SINGLETON",
            ),
        ]
    )
    new = [
        dict(
            cluster_id=10,
            parent_id=None,
            depth=0,
            n_genes=4,
            child_count=2,
            split_status="SPLIT",
            terminal_reason=None,
            selection_kind="BINARY",
        ),
        dict(
            cluster_id=11,
            parent_id=10,
            depth=1,
            n_genes=3,
            child_count=2,
            split_status="SPLIT",
            terminal_reason=None,
            selection_kind="BINARY",
        ),
        dict(
            cluster_id=12,
            parent_id=11,
            depth=2,
            n_genes=2,
            child_count=0,
            split_status="UNRESOLVED",
            terminal_reason="DEPTH_LIMIT",
        ),
        dict(
            cluster_id=13,
            parent_id=11,
            depth=2,
            n_genes=1,
            child_count=0,
            split_status="TERMINAL",
            terminal_reason="SINGLETON",
        ),
        dict(
            cluster_id=14,
            parent_id=10,
            depth=1,
            n_genes=1,
            child_count=0,
            split_status="TERMINAL",
            terminal_reason="SINGLETON",
        ),
    ]
    result = audit_trees(
        new,
        [(0, 12), (1, 12), (2, 13), (3, 14)],
        old,
        [(0, 4), (1, 5), (2, 2), (3, 3)],
        max_depth=2,
    )
    failed = result["failed_nodes"][0]
    assert failed["baseline_exact_cluster"] == 1
    assert failed["exact_clade_depth_delta"] == 1
    assert failed["ancestor_selection_kinds"] == {"BINARY": 2}
    assert result["summary"]["unresolved_proteins"] == 2
    assert result["summary"]["all_depth_limit"]
    assert shrink_steps(226, 0.95) == 59


def test_integer_shrink_bound_includes_natural_terminal_precedence():
    assert shrink_steps(1, 0.95) == 0
    assert shrink_steps(2, 0.95) == 1
    assert shrink_steps(27, 0.95) == 22
    assert 20 + shrink_steps(27, 0.95) == 42
