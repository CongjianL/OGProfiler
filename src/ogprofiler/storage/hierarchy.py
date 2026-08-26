"""Component-level hierarchy result persistence."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from ogprofiler.hierarchy.engine import HierarchyResult


def write_hierarchy_result(directory: Path, result: HierarchyResult) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    nodes = result.nodes
    node_table = pa.table(
        {
            "cluster_id": pa.array([node.cluster_id for node in nodes], type=pa.int64()),
            "parent_id": pa.array([node.parent_id for node in nodes], type=pa.int64()),
            "component_id": pa.array([node.component_id for node in nodes], type=pa.int64()),
            "depth": pa.array([node.depth for node in nodes], type=pa.int32()),
            "n_genes": pa.array([node.n_genes for node in nodes], type=pa.int64()),
            "n_species": pa.array([node.n_species for node in nodes], type=pa.int32()),
            "resolution": pa.array([node.resolution for node in nodes], type=pa.float64()),
            "quality": pa.array([node.quality for node in nodes], type=pa.float64()),
            "child_count": pa.array([node.child_count for node in nodes], type=pa.int32()),
            "split_status": pa.array([node.split_status for node in nodes], type=pa.string()),
            "terminal_reason": pa.array([node.terminal_reason for node in nodes], type=pa.string()),
        }
    )
    pq.write_table(node_table, directory / "nodes.parquet", compression="zstd")
    membership_table = pa.table(
        {
            "protein_id": pa.array(
                [protein_id for protein_id, _ in result.terminal_membership], type=pa.int64()
            ),
            "terminal_cluster_id": pa.array(
                [cluster_id for _, cluster_id in result.terminal_membership], type=pa.int64()
            ),
        }
    )
    pq.write_table(membership_table, directory / "members.parquet", compression="zstd")
    candidates = result.resolution_candidates
    candidate_table = pa.table(
        {
            "cluster_id": pa.array(
                [candidate.cluster_id for candidate in candidates], type=pa.int64()
            ),
            "gamma": pa.array([candidate.gamma for candidate in candidates], type=pa.float64()),
            "child_count": pa.array(
                [candidate.child_count for candidate in candidates], type=pa.int32()
            ),
            "quality": pa.array([candidate.quality for candidate in candidates], type=pa.float64()),
            "min_child_size": pa.array(
                [candidate.min_child_size for candidate in candidates], type=pa.int64()
            ),
            "max_child_fraction": pa.array(
                [candidate.max_child_fraction for candidate in candidates], type=pa.float64()
            ),
            "tiny_fragment_fraction": pa.array(
                [candidate.tiny_fragment_fraction for candidate in candidates], type=pa.float64()
            ),
            "stability": pa.array(
                [candidate.stability for candidate in candidates], type=pa.float64()
            ),
            "adjusted_rand_index": pa.array(
                [candidate.adjusted_rand_index for candidate in candidates], type=pa.float64()
            ),
            "normalized_mutual_info": pa.array(
                [candidate.normalized_mutual_info for candidate in candidates], type=pa.float64()
            ),
            "inter_edge_fraction": pa.array(
                [candidate.inter_edge_fraction for candidate in candidates], type=pa.float64()
            ),
            "intra_edge_fraction": pa.array(
                [candidate.intra_edge_fraction for candidate in candidates], type=pa.float64()
            ),
            "rejection_reason": pa.array(
                [candidate.rejection_reason for candidate in candidates], type=pa.string()
            ),
            "valid": pa.array([candidate.valid for candidate in candidates], type=pa.bool_()),
            "selected": pa.array([candidate.selected for candidate in candidates], type=pa.bool_()),
        }
    )
    pq.write_table(candidate_table, directory / "candidates.parquet", compression="zstd")
    (directory / "metrics.json").write_text(
        json.dumps(asdict(result.metrics), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
