"""Deterministic terminal-family assembly and exchange-table rendering."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import uuid
from collections.abc import Collection, Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol, cast

import pyarrow.parquet as pq

from ogprofiler.exceptions import ExportError

FAMILY_ID_ALGORITHM = "canonical-terminal-membership-rank-v1"


@dataclass(frozen=True, slots=True)
class ProteinMetadata:
    protein_id: int
    species_id: int
    original_id: str


@dataclass(frozen=True, slots=True)
class TerminalFamily:
    family_id: str
    component_id: int
    cluster_id: int
    members: tuple[ProteinMetadata, ...]
    terminal_reason: str
    network_event: str


@dataclass(frozen=True, slots=True)
class ExportTables:
    families: tuple[dict[str, Any], ...]
    members: tuple[dict[str, Any], ...]
    hierarchy: tuple[dict[str, Any], ...]
    events: tuple[dict[str, Any], ...]
    terminal_families: tuple[TerminalFamily, ...]


def _rows(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        raise ExportError(f"Missing export input: {path}")
    try:
        return cast(list[dict[str, Any]], pq.ParquetFile(path).read().to_pylist())
    except Exception as error:
        raise ExportError(f"Failed to read export input {path}: {error}") from error


def _family_sort_key(members: Iterable[ProteinMetadata]) -> tuple[tuple[int, bytes], ...]:
    return tuple(
        (member.species_id, member.original_id.encode("utf-8"))
        for member in sorted(
            members,
            key=lambda item: (item.species_id, item.original_id.encode("utf-8"), item.protein_id),
        )
    )


def _event_index(run_root: Path, component_id: int) -> dict[int, dict[str, Any]]:
    path = (
        run_root / "evolution" / "components" / f"component={component_id:08d}" / "events.parquet"
    )
    rows = _rows(path)
    result = {int(row["cluster_id"]): row for row in rows}
    if len(result) != len(rows):
        raise ExportError(f"Duplicate event cluster IDs for component {component_id}")
    return result


def assemble_export_tables(
    run_root: Path,
    *,
    component_ids: Collection[int] | None = None,
    id_prefix: str = "OG",
) -> ExportTables:
    index_path = run_root / "components/index.parquet"
    if component_ids is None and index_path.is_file():
        component_ids = {int(row["component_id"]) for row in _rows(index_path)}
    proteins = {
        int(row["protein_id"]): ProteinMetadata(
            int(row["protein_id"]), int(row["species_id"]), str(row["original_id"])
        )
        for row in _rows(run_root / "input" / "proteins.parquet")
    }
    if not proteins:
        raise ExportError("Protein metadata is empty")

    family_drafts: list[tuple[int, int, tuple[ProteinMetadata, ...], str, str]] = []
    hierarchy_rows: list[dict[str, Any]] = []
    event_rows: list[dict[str, Any]] = []
    covered: set[int] = set()

    directories = sorted((run_root / "hierarchy" / "components").glob("component=*"))
    for directory in directories:
        component_id = int(directory.name.split("=", 1)[1])
        if component_ids is not None and component_id not in component_ids:
            continue
        nodes = _rows(directory / "nodes.parquet")
        memberships = _rows(directory / "members.parquet")
        events = _event_index(run_root, component_id)
        nodes_by_id = {int(row["cluster_id"]): row for row in nodes}
        if len(nodes_by_id) != len(nodes):
            raise ExportError(f"Duplicate hierarchy cluster IDs for component {component_id}")
        if set(nodes_by_id) != set(events):
            raise ExportError(f"Hierarchy/event node mismatch for component {component_id}")
        members_by_cluster: dict[int, list[ProteinMetadata]] = {}
        for row in memberships:
            protein_id = int(row["protein_id"])
            cluster_id = int(row["terminal_cluster_id"])
            if protein_id in covered:
                raise ExportError(f"Protein {protein_id} belongs to multiple terminal families")
            if protein_id not in proteins:
                raise ExportError(f"Unknown protein ID in hierarchy membership: {protein_id}")
            covered.add(protein_id)
            members_by_cluster.setdefault(cluster_id, []).append(proteins[protein_id])
        terminal_ids = {
            cluster_id
            for cluster_id, node in nodes_by_id.items()
            if node.get("terminal_reason") is not None
        }
        if terminal_ids != set(members_by_cluster):
            raise ExportError(
                f"Terminal hierarchy/membership mismatch for component {component_id}"
            )
        for cluster_id, node in sorted(nodes_by_id.items()):
            event_row = events[cluster_id]
            hierarchy_rows.append(
                {
                    "cluster_id": cluster_id,
                    "parent_id": node.get("parent_id"),
                    "component_id": component_id,
                    "depth": int(node["depth"]),
                    "n_genes": int(node["n_genes"]),
                    "n_species": int(node["n_species"]),
                    "resolution": node.get("resolution"),
                    "quality": node.get("quality"),
                    "child_count": int(node["child_count"]),
                    "terminal_reason": node.get("terminal_reason"),
                    "split_status": node.get("split_status"),
                    "search_status": node.get("search_status"),
                    "termination_kind": node.get("termination_kind"),
                    "failure_codes": ";".join(node.get("failure_codes") or ()),
                    "selection_phase": node.get("selection_phase"),
                }
            )
            event_rows.append(
                {
                    "component_id": component_id,
                    "cluster_id": cluster_id,
                    "network_event": str(event_row["network_event"]),
                    "overlap_score": float(event_row["overlap_score"]),
                    "confidence": float(event_row["confidence"]),
                }
            )
        for cluster_id in sorted(terminal_ids):
            node = nodes_by_id[cluster_id]
            family_drafts.append(
                (
                    component_id,
                    cluster_id,
                    tuple(members_by_cluster[cluster_id]),
                    str(node["terminal_reason"]),
                    str(events[cluster_id]["network_event"]),
                )
            )

    singleton_path = run_root / "components" / "singleton_terminal_families.parquet"
    if singleton_path.is_file():
        for row in _rows(singleton_path):
            protein_id = int(row["protein_id"])
            component_id = int(row["component_id"])
            if protein_id in covered:
                matching = [
                    draft
                    for draft in family_drafts
                    if draft[0] == component_id
                    and len(draft[2]) == 1
                    and draft[2][0].protein_id == protein_id
                ]
                if len(matching) == 1:
                    continue
                raise ExportError(f"Protein {protein_id} is both singleton and hierarchical")
            if protein_id not in proteins:
                raise ExportError(f"Unknown singleton protein ID: {protein_id}")
            covered.add(protein_id)
            member = proteins[protein_id]
            reason = str(row.get("terminal_reason", "SINGLETON"))
            family_drafts.append((component_id, 0, (member,), reason, "SPECIES_SPECIFIC"))
            hierarchy_rows.append(
                {
                    "cluster_id": 0,
                    "parent_id": None,
                    "component_id": component_id,
                    "depth": 0,
                    "n_genes": 1,
                    "n_species": 1,
                    "resolution": None,
                    "quality": None,
                    "child_count": 0,
                    "terminal_reason": reason,
                }
            )
            event_rows.append(
                {
                    "component_id": component_id,
                    "cluster_id": 0,
                    "network_event": "SPECIES_SPECIFIC",
                    "overlap_score": 0.0,
                    "confidence": 1.0,
                }
            )

    missing = sorted(set(proteins) - covered)
    if missing:
        preview = ", ".join(str(value) for value in missing[:5])
        raise ExportError(f"Proteins without terminal-family assignment: {preview}")
    if len(covered) != len(proteins):
        raise ExportError("Terminal-family membership is not a protein partition")

    ordered = sorted(family_drafts, key=lambda item: _family_sort_key(item[2]))
    terminal_families: list[TerminalFamily] = []
    family_rows: list[dict[str, Any]] = []
    member_rows: list[dict[str, Any]] = []
    for rank, (component_id, cluster_id, members, reason, network_event) in enumerate(ordered):
        family_id = f"{id_prefix}{rank:09d}"
        sorted_members = tuple(
            sorted(
                members,
                key=lambda item: (
                    item.species_id,
                    item.original_id.encode("utf-8"),
                    item.protein_id,
                ),
            )
        )
        family = TerminalFamily(
            family_id, component_id, cluster_id, sorted_members, reason, network_event
        )
        terminal_families.append(family)
        family_rows.append(
            {
                "family_id": family_id,
                "component_id": component_id,
                "cluster_id": cluster_id,
                "n_genes": len(sorted_members),
                "n_species": len({member.species_id for member in sorted_members}),
                "terminal_reason": reason,
                "network_event": network_event,
            }
        )
        member_rows.extend(
            {
                "family_id": family_id,
                "protein_id": member.protein_id,
                "species_id": member.species_id,
                "original_id": member.original_id,
            }
            for member in sorted_members
        )

    hierarchy_rows.sort(key=lambda row: (int(row["component_id"]), int(row["cluster_id"])))
    event_rows.sort(key=lambda row: (int(row["component_id"]), int(row["cluster_id"])))
    return ExportTables(
        tuple(family_rows),
        tuple(member_rows),
        tuple(hierarchy_rows),
        tuple(event_rows),
        tuple(terminal_families),
    )


def write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


class MemberCollection(Protocol):
    @property
    def members(self) -> tuple[ProteinMetadata, ...]: ...


def family_membership_sha256(family: MemberCollection) -> str:
    value = [[member.species_id, member.original_id] for member in family.members]
    encoded = json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def terminal_table_paths(run_root: Path) -> tuple[Path, Path]:
    """Exchange input selection for hierarchy diagnostics, never final OG membership."""
    root = run_root / "results"
    if (root / "terminal_families.tsv").is_file():
        return root / "terminal_families.tsv", root / "terminal_members.tsv"
    diagnostic = root / "terminal-families"
    if (diagnostic / "export-manifest.json").is_file():
        return diagnostic / "families.tsv", diagnostic / "members.tsv"
    return root / "families.tsv", root / "members.tsv"
