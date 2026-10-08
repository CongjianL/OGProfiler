"""Component-local network-event annotation and verified resume."""

from __future__ import annotations

import json
import os
import uuid
from collections import Counter
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from ogprofiler.core.manifest import sha256_file, sha256_json, write_json
from ogprofiler.evolution.network import EventAnnotation, annotate_network_events
from ogprofiler.exceptions import HierarchyError

NETWORK_EVENT_ALGORITHM_VERSION = "network-overlap-v1"
EVENT_SCHEMA = pa.schema(
    [
        ("component_id", pa.int64()),
        ("cluster_id", pa.int64()),
        ("child_count", pa.int32()),
        ("network_event", pa.string()),
        ("overlap_count", pa.int32()),
        ("overlap_score", pa.float64()),
        ("pairwise_overlap_summary", pa.string()),
        ("confidence", pa.float64()),
        ("legacy_event", pa.string()),
    ]
)


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object in {path}")
    return value


def _write_events(path: Path, annotations: tuple[EventAnnotation, ...]) -> None:
    rows = [
        {
            "component_id": item.component_id,
            "cluster_id": item.cluster_id,
            "child_count": item.child_count,
            "network_event": item.network_event,
            "overlap_count": item.overlap_count,
            "overlap_score": item.overlap_score,
            "pairwise_overlap_summary": item.pairwise_overlap_summary,
            "confidence": item.confidence,
            "legacy_event": item.legacy_event,
        }
        for item in annotations
    ]
    pq.write_table(pa.Table.from_pylist(rows, schema=EVENT_SCHEMA), path, compression="zstd")


def annotate_component(
    run_root: Path, component_id: int, overlap_threshold: float, command: list[str]
) -> tuple[Path, bool, Counter[str]]:
    hierarchy = run_root / "hierarchy" / "components" / f"component={component_id:08d}"
    nodes_path = hierarchy / "nodes.parquet"
    members_path = hierarchy / "members.parquet"
    proteins_path = run_root / "input" / "proteins.parquet"
    if not nodes_path.is_file() or not members_path.is_file() or not proteins_path.is_file():
        raise HierarchyError(f"Missing hierarchy inputs for component {component_id}")
    output = run_root / "evolution" / "components" / f"component={component_id:08d}"
    events_path = output / "events.parquet"
    manifest_path = output / "network-event-manifest.json"
    parameters = {"component_id": component_id, "overlap_threshold": overlap_threshold}
    inputs = {
        path.relative_to(run_root).as_posix(): sha256_file(path)
        for path in (nodes_path, members_path, proteins_path)
    }
    if events_path.is_file() and manifest_path.is_file():
        try:
            previous = _load_json(manifest_path)
            if bool(
                previous["algorithm_version"] == NETWORK_EVENT_ALGORITHM_VERSION
                and previous["parameters"] == parameters
                and previous["input_checksums"] == inputs
                and previous["output_checksums"]["events.parquet"] == sha256_file(events_path)
            ):
                counts = Counter(
                    str(value) for value in pq.read_table(events_path)["network_event"].to_pylist()
                )
                return events_path, True, counts
        except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
            pass
    nodes = pq.read_table(nodes_path).to_pylist()
    members = pq.read_table(members_path).to_pylist()
    proteins = pq.read_table(proteins_path, columns=["protein_id", "species_id"]).to_pylist()
    species = {int(row["protein_id"]): int(row["species_id"]) for row in proteins}
    annotations = annotate_network_events(nodes, members, species, overlap_threshold)
    output.mkdir(parents=True, exist_ok=True)
    temporary = output / f".events.{uuid.uuid4().hex}.tmp.parquet"
    _write_events(temporary, annotations)
    os.replace(temporary, events_path)
    counts = Counter(item.network_event for item in annotations)
    write_json(
        manifest_path,
        {
            "algorithm_version": NETWORK_EVENT_ALGORITHM_VERSION,
            "command": command,
            "parameters": parameters,
            "parameters_sha256": sha256_json(parameters),
            "input_checksums": inputs,
            "output_checksums": {"events.parquet": sha256_file(events_path)},
            "counts": {"nodes": len(annotations), "events": dict(sorted(counts.items()))},
        },
    )
    return events_path, False, counts


def run_network_annotation_stage(
    run_root: Path, overlap_threshold: float, command: list[str]
) -> tuple[Path, int, int]:
    hierarchy_root = run_root / "hierarchy" / "components"
    component_dirs = sorted(hierarchy_root.glob("component=*"))
    if not component_dirs:
        # An edge-free dataset has no non-singleton hierarchy to annotate.
        # Accept only a resolved scheduler and exactly matching singleton inputs.
        try:
            scheduler = _load_json(run_root / "hierarchy/scheduler-manifest.json")
            index = pq.ParquetFile(run_root / "components/index.parquet").read().to_pylist()
            singletons = (
                pq.ParquetFile(run_root / "components/singleton_terminal_families.parquet")
                .read()
                .to_pylist()
            )
            proteins = pq.ParquetFile(run_root / "input/proteins.parquet").read().to_pylist()
            pairs = [(int(r["protein_id"]), int(r["component_id"])) for r in index]
            singleton_pairs = [(int(r["protein_id"]), int(r["component_id"])) for r in singletons]
            valid = (
                scheduler["hierarchy_status"] == "RESOLVED"
                and scheduler["counts"]["failed"] == 0
                and scheduler["counts"]["unresolved"] == 0
                and scheduler["counts"]["singleton_components"] == len(pairs)
                and bool(pairs)
                and len({p for p, _ in pairs}) == len(pairs)
                and len({c for _, c in pairs}) == len(pairs)
                and sorted(pairs) == sorted(singleton_pairs)
                and sorted(p for p, _ in pairs) == sorted(int(r["protein_id"]) for r in proteins)
            )
        except (OSError, KeyError, TypeError, ValueError):
            valid = False
        if not valid:
            raise HierarchyError("Network annotation requires component hierarchy outputs")
    total = Counter[str]()
    reused = 0
    component_ids: list[int] = []
    for directory in component_dirs:
        component_id = int(directory.name.split("=", 1)[1])
        component_ids.append(component_id)
        _, was_reused, counts = annotate_component(
            run_root, component_id, overlap_threshold, command
        )
        reused += int(was_reused)
        total.update(counts)
    manifest_path = run_root / "evolution" / "network-event-manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    write_json(
        manifest_path,
        {
            "algorithm_version": NETWORK_EVENT_ALGORITHM_VERSION,
            "command": command,
            "overlap_threshold": overlap_threshold,
            "component_ids": component_ids,
            "component_count": len(component_ids),
            "reused_components": reused,
            "event_counts": dict(sorted(total.items())),
        },
    )
    return manifest_path, len(component_ids), reused
