#!/usr/bin/env python3
"""Compare normalized frozen V1 hierarchy with a V2 prototype run."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq


def member_hash(members: set[str]) -> str:
    return hashlib.sha256("\n".join(sorted(members)).encode()).hexdigest()[:20]


def read_id_map(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            original_id, recoded_id = line.split("\t")[:2]
            mapping[recoded_id] = original_id
    return mapping


def normalize_v2(run: Path, id_map: dict[str, str]) -> dict[str, Any]:
    all_nodes: list[dict[str, Any]] = []
    all_membership: dict[str, str] = {}
    component_root = run / "hierarchy" / "components"
    for component_dir in sorted(path for path in component_root.iterdir() if path.is_dir()):
        component_id = int(component_dir.name)
        nodes = pq.read_table(component_dir / "nodes.parquet").to_pylist()
        memberships = pq.read_table(component_dir / "members.parquet").to_pylist()
        vertex_path = run / "components" / f"component_{component_id:08d}.vertices.parquet"
        vertices = pq.read_table(vertex_path).to_pylist()
        legacy_by_integer = {
            int(row["protein_id"]): str(row["original_id"])
            for row in vertices
        }
        node_by_id = {int(node["cluster_id"]): node for node in nodes}
        descendants: dict[int, set[str]] = {int(node["cluster_id"]): set() for node in nodes}
        terminal_by_protein = {
            int(row["protein_id"]): int(row["terminal_cluster_id"])
            for row in memberships
        }
        for protein_id, terminal_id in terminal_by_protein.items():
            legacy_id = legacy_by_integer[protein_id]
            original_id = id_map.get(legacy_id, legacy_id)
            current_id: int | None = terminal_id
            while current_id is not None:
                descendants[current_id].add(original_id)
                parent = node_by_id[current_id]["parent_id"]
                current_id = int(parent) if parent is not None else None
        hashes = {cluster_id: member_hash(members) for cluster_id, members in descendants.items()}
        for cluster_id, node in node_by_id.items():
            parent = node["parent_id"]
            all_nodes.append(
                {
                    "component_id": component_id,
                    "cluster_id": cluster_id,
                    "cluster_hash": hashes[cluster_id],
                    "parent_hashes": [hashes[int(parent)]] if parent is not None else [],
                    "depth": int(node["depth"]),
                    "n_genes": len(descendants[cluster_id]),
                    "members": sorted(descendants[cluster_id]),
                    "terminal": node["terminal_reason"] is not None,
                }
            )
        for protein_id, terminal_id in terminal_by_protein.items():
            legacy_id = legacy_by_integer[protein_id]
            original_id = id_map.get(legacy_id, legacy_id)
            all_membership[original_id] = hashes[terminal_id]
    all_nodes.sort(key=lambda item: (item["component_id"], item["depth"], item["cluster_hash"]))
    return {"nodes": all_nodes, "membership": all_membership}


def pair_set(membership: dict[str, str]) -> set[tuple[str, str]]:
    by_cluster: dict[str, list[str]] = {}
    for protein_id, cluster_id in membership.items():
        by_cluster.setdefault(cluster_id, []).append(protein_id)
    pairs: set[tuple[str, str]] = set()
    for members in by_cluster.values():
        members.sort()
        for left_index, left in enumerate(members):
            pairs.update((left, right) for right in members[left_index + 1 :])
    return pairs


def ratio(numerator: float | int, denominator: float | int) -> float | None:
    return float(numerator) / float(denominator) if denominator else None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--v1-normalized", required=True, type=Path)
    parser.add_argument("--v1-hierarchy-metrics", required=True, type=Path)
    parser.add_argument("--v2-run", required=True, type=Path)
    parser.add_argument("--sequence-ids", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    v1 = json.loads(args.v1_normalized.read_text(encoding="utf-8"))
    v1_metrics = json.loads(args.v1_hierarchy_metrics.read_text(encoding="utf-8"))
    id_map = read_id_map(args.sequence_ids)
    v2 = normalize_v2(args.v2_run, id_map)
    (args.out / "v2_normalized_hierarchy.json").write_text(
        json.dumps(v2, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    v2_manifest = json.loads((args.v2_run / "manifest.json").read_text(encoding="utf-8"))
    v2_components = v2_manifest["components"]
    v2_runtime = sum(float(component["runtime_seconds"]) for component in v2_components)
    v2_peak_rss = max(int(component["peak_rss_bytes"]) for component in v2_components)

    v1_membership = {str(key): str(value) for key, value in v1["membership"].items()}
    v2_membership = {str(key): str(value) for key, value in v2["membership"].items()}
    common_proteins = set(v1_membership) & set(v2_membership)
    v1_common = {protein: v1_membership[protein] for protein in common_proteins}
    v2_common = {protein: v2_membership[protein] for protein in common_proteins}
    v1_pairs = pair_set(v1_common)
    v2_pairs = pair_set(v2_common)
    pair_intersection = v1_pairs & v2_pairs
    membership_rows = [
        {
            "protein_id": protein,
            "v1_terminal_cluster": v1_membership.get(protein),
            "v2_terminal_cluster": v2_membership.get(protein),
            "exact_terminal_set_match": v1_membership.get(protein) == v2_membership.get(protein),
        }
        for protein in sorted(set(v1_membership) | set(v2_membership))
    ]
    membership_path = args.out / "membership_comparison.tsv"
    with membership_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(
            "protein_id\tv1_terminal_cluster\tv2_terminal_cluster\texact_terminal_set_match\n"
        )
        for row in membership_rows:
            handle.write(
                f"{row['protein_id']}\t{row['v1_terminal_cluster'] or ''}\t"
                f"{row['v2_terminal_cluster'] or ''}\t"
                f"{str(row['exact_terminal_set_match']).lower()}\n"
            )

    v1_node_sets = {node["cluster_hash"] for node in v1["nodes"]}
    v2_node_sets = {node["cluster_hash"] for node in v2["nodes"]}
    v1_edges = {
        (parent, node["cluster_hash"])
        for node in v1["nodes"]
        for parent in node["parent_hashes"]
    }
    v2_edges = {
        (parent, node["cluster_hash"])
        for node in v2["nodes"]
        for parent in node["parent_hashes"]
    }
    comparison = {
        "proteins": {
            "v1": len(v1_membership),
            "v2": len(v2_membership),
            "common": len(common_proteins),
            "exact_terminal_set_assignments": sum(
                v1_membership.get(protein) == v2_membership.get(protein)
                for protein in set(v1_membership) | set(v2_membership)
            ),
        },
        "terminal_pairing": {
            "v1_pairs": len(v1_pairs),
            "v2_pairs": len(v2_pairs),
            "intersection": len(pair_intersection),
            "precision": ratio(len(pair_intersection), len(v2_pairs)),
            "recall": ratio(len(pair_intersection), len(v1_pairs)),
            "jaccard": ratio(len(pair_intersection), len(v1_pairs | v2_pairs)),
        },
        "topology": {
            "v1_nodes": len(v1_node_sets),
            "v2_nodes": len(v2_node_sets),
            "shared_member_sets": len(v1_node_sets & v2_node_sets),
            "node_jaccard": ratio(
                len(v1_node_sets & v2_node_sets), len(v1_node_sets | v2_node_sets)
            ),
            "v1_edges": len(v1_edges),
            "v2_edges": len(v2_edges),
            "shared_edges": len(v1_edges & v2_edges),
            "edge_jaccard": ratio(len(v1_edges & v2_edges), len(v1_edges | v2_edges)),
        },
        "performance": {
            "v1_runtime_seconds": v1_metrics["runtime_seconds"],
            "v2_runtime_seconds": v2_runtime,
            "runtime_speedup_v1_over_v2": ratio(v1_metrics["runtime_seconds"], v2_runtime),
            "v1_peak_rss_bytes": v1_metrics["peak_rss_bytes"],
            "v2_peak_rss_bytes": v2_peak_rss,
            "rss_ratio_v1_over_v2": ratio(v1_metrics["peak_rss_bytes"], v2_peak_rss),
            "v2_leiden_calls": sum(int(item["leiden_calls"]) for item in v2_components),
            "v2_subgraph_constructions": sum(
                int(item["subgraph_constructions"]) for item in v2_components
            ),
        },
    }
    (args.out / "comparison.json").write_text(
        json.dumps(comparison, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(comparison, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
