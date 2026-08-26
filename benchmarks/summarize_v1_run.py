#!/usr/bin/env python3
"""Normalize frozen V1 SSN/hierarchy outputs and summarize baseline metrics."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import defaultdict, deque
from pathlib import Path
from typing import Any

import igraph


def member_hash(members: set[str]) -> str:
    payload = "\n".join(sorted(members)).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()[:20]


def read_id_map(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        original_id, recoded_id = line.split("\t")[:2]
        mapping[recoded_id] = original_id
    return mapping


def parse_elapsed(value: str) -> float:
    parts = value.strip().split(":")
    if len(parts) == 2:
        minutes, seconds = parts
        return int(minutes) * 60 + float(seconds)
    if len(parts) == 3:
        hours, minutes, seconds = parts
        return int(hours) * 3600 + int(minutes) * 60 + float(seconds)
    return float(value)


def parse_time_report(path: Path) -> dict[str, float | int]:
    values: dict[str, float | int] = {}
    if not path.is_file():
        return values
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "Elapsed (wall clock) time" in line:
            values["runtime_seconds"] = parse_elapsed(line.rsplit(": ", 1)[-1])
        elif "Maximum resident set size" in line:
            values["peak_rss_bytes"] = int(line.rsplit(":", 1)[-1].strip()) * 1024
    return values


def parse_stage_times(path: Path) -> dict[str, float]:
    if not path.is_file():
        return {}
    pattern = re.compile(r"^\[(.+?)\]:\s+(\d+)h\s+(\d+)m\s+([0-9.]+)s")
    result: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = pattern.search(line.strip())
        if match:
            name, hours, minutes, seconds = match.groups()
            result[name] = int(hours) * 3600 + int(minutes) * 60 + float(seconds)
    return result


def normalize_hierarchy(graph: igraph.Graph, id_map: dict[str, str]) -> dict[str, Any]:
    raw_members: list[set[str]] = []
    for vertex in graph.vs:
        recoded = str(vertex["geneIDs"] or "").split()
        raw_members.append({id_map.get(item, item) for item in recoded})

    parents: dict[int, set[int]] = defaultdict(set)
    children: dict[int, set[int]] = defaultdict(set)
    invalid_edges: list[tuple[int, int]] = []
    for edge in graph.es:
        left, right = edge.tuple
        left_members, right_members = raw_members[left], raw_members[right]
        if left_members > right_members:
            parent, child = left, right
        elif right_members > left_members:
            parent, child = right, left
        elif len(left_members) > len(right_members):
            parent, child = left, right
        elif len(right_members) > len(left_members):
            parent, child = right, left
        else:
            invalid_edges.append((left, right))
            continue
        parents[child].add(parent)
        children[parent].add(child)

    roots = [index for index in range(graph.vcount()) if not parents[index]]
    depths: dict[int, int] = {root: 0 for root in roots}
    queue = deque(roots)
    while queue:
        parent = queue.popleft()
        for child in sorted(children[parent]):
            candidate = depths[parent] + 1
            if child not in depths or candidate < depths[child]:
                depths[child] = candidate
                queue.append(child)

    nodes: list[dict[str, Any]] = []
    for index, vertex in enumerate(graph.vs):
        members = raw_members[index]
        parent_hashes = sorted(member_hash(raw_members[parent]) for parent in parents[index])
        nodes.append(
            {
                "cluster_hash": member_hash(members),
                "v1_name": str(vertex["name"]),
                "parent_hashes": parent_hashes,
                "depth": depths.get(index),
                "n_genes": len(members),
                "members": sorted(members),
                "terminal": not children[index],
            }
        )
    nodes.sort(key=lambda item: (item["depth"] is None, item["depth"], item["cluster_hash"]))
    terminals = [node for node in nodes if node["terminal"]]
    membership = {
        member: node["cluster_hash"]
        for node in terminals
        for member in node["members"]
    }
    return {
        "nodes": nodes,
        "membership": membership,
        "root_count": len(roots),
        "terminal_family_count": len(terminals),
        "max_depth": max(depths.values(), default=0),
        "invalid_unoriented_edges": invalid_edges,
        "multi_parent_nodes": sorted(index for index, values in parents.items() if len(values) > 1),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--working-dir", required=True, type=Path)
    parser.add_argument("--time-report", required=True, type=Path)
    parser.add_argument("--stdout", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    ssn = igraph.Graph.Read_GML(str(args.working_dir / "ssn.gml"))
    hierarchy = igraph.Graph.Read_GML(str(args.working_dir / "hm.gml"))
    id_map = read_id_map(args.working_dir / "SequenceIDs.txt")
    normalized = normalize_hierarchy(hierarchy, id_map)
    (args.out / "normalized_hierarchy.json").write_text(
        json.dumps(normalized, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    components = ssn.connected_components()
    metrics: dict[str, Any] = {
        **parse_time_report(args.time_report),
        "stage_runtime_seconds": parse_stage_times(args.stdout),
        "ssn_nodes": ssn.vcount(),
        "ssn_edges": ssn.ecount(),
        "connected_component_count": len(components),
        "largest_component_size": max(components.sizes(), default=0),
        "hierarchy_nodes": hierarchy.vcount(),
        "hierarchy_edges": hierarchy.ecount(),
        "hierarchy_depth": normalized["max_depth"],
        "terminal_family_count": normalized["terminal_family_count"],
        "membership_count": len(normalized["membership"]),
    }
    (args.out / "baseline_metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(metrics, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

