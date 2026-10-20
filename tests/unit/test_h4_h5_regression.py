"""Compact Slurm result summaries retain unresolved semantics."""

import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from benchmarks.og_extraction.hierarchy_regression import summarize_hierarchy


def test_summary_counts_unresolved_proteins_and_budget(tmp_path: Path):
    root = tmp_path / "hierarchy/components/component=00000000"
    root.mkdir(parents=True)
    pq.write_table(
        pa.Table.from_pylist(
            [
                dict(
                    cluster_id=0,
                    parent_id=None,
                    depth=0,
                    n_genes=6,
                    split_status="SPLIT",
                    terminal_reason=None,
                    search_status="ACCEPTED",
                ),
                dict(
                    cluster_id=1,
                    parent_id=0,
                    depth=1,
                    n_genes=1,
                    split_status="TERMINAL",
                    terminal_reason="SINGLETON",
                    search_status="POLICY_STOP",
                ),
                dict(
                    cluster_id=2,
                    parent_id=0,
                    depth=1,
                    n_genes=5,
                    split_status="UNRESOLVED",
                    terminal_reason="REJECTED_ALL_TESTED",
                    search_status="REJECTED_ALL_TESTED",
                ),
            ]
        ),
        root / "nodes.parquet",
    )
    pq.write_table(
        pa.Table.from_pylist(
            [
                dict(
                    cluster_id=0,
                    gamma=0.1,
                    selected=True,
                    valid=True,
                    violations=[],
                    stability=0.95,
                    stability_evaluated=True,
                    phase="coarse",
                    evaluation_budget=24,
                )
            ]
        ),
        root / "candidates.parquet",
    )
    (root / "metrics.json").write_text(json.dumps(dict(leiden_calls=3)))
    result = summarize_hierarchy(tmp_path, (0,))
    assert result["root_split"]
    assert result["unresolved_proteins"] == 5
    assert result["unresolved_nodes"] == 1
    assert not result["resolved"]
    assert result["candidate_budget_passed"]
    assert result["roots"][0]["selected_candidates"][0]["stability"] == 0.95


def test_config_migration_preserves_scientific_controls():
    from copy import deepcopy

    from benchmarks.og_extraction.hierarchy_regression import migrate_config
    from ogprofiler.config import DEFAULT_CONFIG

    old = deepcopy(DEFAULT_CONFIG)
    old["hierarchy"].update(min_family_size=2, resolution_strategy="adaptive", leiden_iterations=2)
    new = migrate_config(old)
    assert new["hierarchy"]["min_family_size"] is None
    assert new["hierarchy"]["resolution_strategy"] == "bounded_adaptive_v2"
    assert new["hierarchy"]["leiden_iterations"] == 10
    assert new["hierarchy"]["recursion_stop_size"] == 1
    assert new["hierarchy"]["max_candidate_evaluations"] == 24
    for key in ("seed", "method", "stability_threshold", "stability_mode", "max_depth"):
        assert new["hierarchy"][key] == old["hierarchy"][key]
    assert new["edges"] == old["edges"]
