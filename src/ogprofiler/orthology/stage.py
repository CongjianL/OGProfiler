"""Checksummed streaming Zstandard ortholog-candidate output stage."""

from __future__ import annotations

import json
import os
import uuid
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any, cast

import pyarrow as pa
import pyarrow.parquet as pq

from ogprofiler.core.manifest import sha256_file, sha256_json, write_json
from ogprofiler.exceptions import HierarchyError
from ogprofiler.orthology.engine import (
    ORTHOLOGY_SUPPORTING_EVENTS,
    OrthologCandidate,
    generate_component_candidates,
)

ORTHOLOGY_ALGORITHM_VERSION = "hierarchy-cross-child-v1"
OUTPUT_NAME = "ortholog_pairs.tsv.zst"
HEADER = (
    "protein_a_id\tprotein_b_id\tspecies_a_id\tspecies_b_id\tcomponent_id\t"
    "supporting_cluster_id\trelationship\n"
)


def _rows(path: Path) -> list[dict[str, Any]]:
    try:
        return cast(list[dict[str, Any]], pq.read_table(path).to_pylist())
    except Exception as error:
        raise HierarchyError(f"Failed to read orthology input {path}: {error}") from error


def _component_paths(run_root: Path) -> list[tuple[int, Path, Path, Path]]:
    paths: list[tuple[int, Path, Path, Path]] = []
    for directory in sorted((run_root / "hierarchy" / "components").glob("component=*")):
        component_id = int(directory.name.split("=", 1)[1])
        nodes = directory / "nodes.parquet"
        members = directory / "members.parquet"
        events = (
            run_root
            / "evolution"
            / "components"
            / directory.name
            / "events.parquet"
        )
        for path in (nodes, members, events):
            if not path.is_file():
                raise HierarchyError(f"Missing orthology input: {path}")
        paths.append((component_id, nodes, members, events))
    return paths


def _candidate_stream(
    component_paths: list[tuple[int, Path, Path, Path]],
    species_by_protein: dict[int, int],
) -> Iterator[OrthologCandidate]:
    for component_id, nodes, members, events in component_paths:
        yield from generate_component_candidates(
            component_id=component_id,
            nodes=_rows(nodes),
            terminal_memberships=_rows(members),
            events=_rows(events),
            species_by_protein=species_by_protein,
        )


def _encode(candidate: OrthologCandidate) -> bytes:
    return (
        f"{candidate.protein_a_id}\t{candidate.protein_b_id}\t"
        f"{candidate.species_a_id}\t{candidate.species_b_id}\t"
        f"{candidate.component_id}\t{candidate.supporting_cluster_id}\t"
        f"{candidate.relationship}\n"
    ).encode()


def write_candidate_stream(
    path: Path,
    candidates: Iterator[OrthologCandidate],
    *,
    chunk_size: int,
    progress: Callable[[int, int], None] | None = None,
) -> tuple[int, int]:
    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    count = 0
    chunks = 0
    buffer = bytearray(HEADER.encode("utf-8"))
    try:
        with pa.output_stream(str(temporary), compression="zstd") as output:
            for candidate in candidates:
                buffer.extend(_encode(candidate))
                count += 1
                if count % chunk_size == 0:
                    output.write(buffer)
                    buffer.clear()
                    chunks += 1
                    if progress is not None:
                        progress(count, chunks)
            if buffer:
                output.write(buffer)
                chunks += 1
                if progress is not None:
                    progress(count, chunks)
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return count, chunks


def _verified(
    manifest_path: Path,
    output_path: Path,
    inputs: dict[str, str],
    parameters: dict[str, Any],
) -> bool:
    if not manifest_path.is_file() or not output_path.is_file():
        return False
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        return bool(
            manifest["algorithm_version"] == ORTHOLOGY_ALGORITHM_VERSION
            and manifest["input_checksums"] == inputs
            and manifest["parameters"] == parameters
            and manifest["output_checksums"][OUTPUT_NAME] == sha256_file(output_path)
        )
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
        return False


def run_orthology_stage(
    run_root: Path,
    command: list[str],
    *,
    chunk_size: int,
    progress: Callable[[int, int], None] | None = None,
) -> tuple[Path, bool, int]:
    component_paths = _component_paths(run_root)
    proteins_path = run_root / "input" / "proteins.parquet"
    if not proteins_path.is_file():
        raise HierarchyError(f"Missing orthology input: {proteins_path}")
    all_inputs = [proteins_path]
    singleton_path = run_root / "components" / "singleton_terminal_families.parquet"
    if singleton_path.is_file():
        all_inputs.append(singleton_path)
    if not component_paths and not singleton_path.is_file():
        raise HierarchyError("Orthology generation requires hierarchy or singleton outputs")
    for _, nodes, members, events in component_paths:
        all_inputs.extend((nodes, members, events))
    inputs = {path.relative_to(run_root).as_posix(): sha256_file(path) for path in all_inputs}
    parameters = {
        "chunk_size": chunk_size,
        "events": sorted(ORTHOLOGY_SUPPORTING_EVENTS),
        "same_species_filter": True,
        "relationship": "CO_ORTHOLOG_CANDIDATE",
    }
    output_root = run_root / "results"
    output_root.mkdir(parents=True, exist_ok=True)
    output_path = output_root / OUTPUT_NAME
    manifest_path = output_root / "ortholog-manifest.json"
    if _verified(manifest_path, output_path, inputs, parameters):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        return manifest_path, True, int(manifest["counts"]["pairs"])

    protein_rows = _rows(proteins_path)
    species = {int(row["protein_id"]): int(row["species_id"]) for row in protein_rows}
    pair_count, chunks = write_candidate_stream(
        output_path,
        _candidate_stream(component_paths, species),
        chunk_size=chunk_size,
        progress=progress,
    )
    supporting_event_counts = {event: 0 for event in sorted(ORTHOLOGY_SUPPORTING_EVENTS)}
    for _, _, _, events in component_paths:
        for row in _rows(events):
            event = str(row["network_event"])
            if event in supporting_event_counts:
                supporting_event_counts[event] += 1
    write_json(
        manifest_path,
        {
            "algorithm_version": ORTHOLOGY_ALGORITHM_VERSION,
            "command": command,
            "parameters": parameters,
            "parameters_sha256": sha256_json(parameters),
            "input_checksums": inputs,
            "output_checksums": {OUTPUT_NAME: sha256_file(output_path)},
            "counts": {
                "pairs": pair_count,
                "chunks": chunks,
                "components": len(component_paths),
                "supporting_nodes": sum(supporting_event_counts.values()),
                "supporting_event_counts": supporting_event_counts,
            },
        },
    )
    return manifest_path, False, pair_count
