"""Final OG assembly distinct from the terminal-family diagnostic partition."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ogprofiler.exceptions import ExportError, HierarchyError
from ogprofiler.orthogroups.models import OrthogroupConfig
from ogprofiler.orthogroups.stage import verified_orthogroup_inputs
from ogprofiler.output.results import (
    ExportTables,
    ProteinMetadata,
    _family_sort_key,
    _rows,
    assemble_export_tables,
    family_membership_sha256,
)

OG_ID_ALGORITHM = "canonical-orthogroup-membership-rank-v1"


@dataclass(frozen=True, slots=True)
class ExportOrthogroup:
    family_id: str
    component_id: int
    local_group_id: int
    source_cluster_id: int | None
    members: tuple[ProteinMetadata, ...]
    selection_type: str
    v1_event: str | None
    processing_level: int
    membership_hash: str


@dataclass(frozen=True, slots=True)
class OrthogroupExportTables:
    families: tuple[dict[str, Any], ...]
    members: tuple[dict[str, Any], ...]
    groups: tuple[ExportOrthogroup, ...]
    unassigned: tuple[dict[str, Any], ...]
    diagnostics: ExportTables
    statistics: tuple[dict[str, Any], ...]


def assemble_orthogroup_export_tables(
    run_root: Path,
    config: OrthogroupConfig | None = None,
) -> OrthogroupExportTables:
    try:
        components, _ = verified_orthogroup_inputs(run_root, config)
    except HierarchyError as error:
        raise ExportError(str(error)) from error
    diagnostic = assemble_export_tables(run_root, component_ids=components, id_prefix="TF")
    protein_rows = _rows(run_root / "input/proteins.parquet")
    proteins = {
        int(row["protein_id"]): ProteinMetadata(
            int(row["protein_id"]), int(row["species_id"]), str(row["original_id"])
        )
        for row in protein_rows
    }
    if len(proteins) != len(protein_rows):
        raise ExportError("Duplicate protein metadata")
    index = {
        int(row["protein_id"]): int(row["component_id"])
        for row in _rows(run_root / "components/index.parquet")
    }
    assigned: set[int] = set()
    unassigned_ids: set[int] = set()
    unassigned = []
    drafts = []
    for component in components:
        directory = run_root / "orthogroups/components" / f"component={component:08d}"
        rows = _rows(directory / "groups.parquet")
        groups_by_id = {int(row["local_group_id"]): row for row in rows}
        if len(rows) != len(groups_by_id):
            raise ExportError(f"Duplicate local OG IDs in component {component}")
        members_by_group: dict[int, list[ProteinMetadata]] = defaultdict(list)
        for row in _rows(directory / "members.parquet"):
            protein = int(row["protein_id"])
            local = int(row["local_group_id"])
            if (
                int(row["component_id"]) != component
                or index.get(protein) != component
                or local not in groups_by_id
                or protein in assigned
            ):
                raise ExportError(f"Invalid or overlapping OG member in component {component}")
            assigned.add(protein)
            members_by_group[local].append(proteins[protein])
        for local, row in groups_by_id.items():
            members = tuple(
                sorted(
                    members_by_group[local],
                    key=lambda p: (p.species_id, p.original_id.encode("utf-8")),
                )
            )
            if (
                int(row["component_id"]) != component
                or not members
                or int(row["n_genes"]) != len(members)
                or int(row["n_species"]) != len({p.species_id for p in members})
            ):
                raise ExportError(f"OG counts/membership mismatch in component {component}")
            draft = ExportOrthogroup(
                "",
                component,
                local,
                row["source_cluster_id"],
                members,
                str(row["selection_type"]),
                row["v1_event"],
                int(row["processing_level"]),
                str(row["membership_hash"]),
            )
            if family_membership_sha256(draft) != draft.membership_hash:
                raise ExportError(f"OG membership hash mismatch in component {component}")
            drafts.append(draft)
        for row in _rows(directory / "unassigned.parquet"):
            protein = int(row["protein_id"])
            if (
                int(row["component_id"]) != component
                or index.get(protein) != component
                or protein in unassigned_ids
                or protein in assigned
            ):
                raise ExportError(f"Invalid OG unassigned member in component {component}")
            unassigned_ids.add(protein)
            metadata = proteins[protein]
            unassigned.append(
                {**row, "species_id": metadata.species_id, "original_id": metadata.original_id}
            )
    if assigned & unassigned_ids or assigned | unassigned_ids != set(proteins):
        raise ExportError("OG assigned/unassigned sets do not partition prepared proteins")
    groups: list[ExportOrthogroup] = []
    family_rows: list[dict[str, Any]] = []
    member_rows: list[dict[str, Any]] = []
    for rank, draft in enumerate(sorted(drafts, key=lambda group: _family_sort_key(group.members))):
        family_id = f"OG{rank:09d}"
        group = ExportOrthogroup(
            family_id,
            draft.component_id,
            draft.local_group_id,
            draft.source_cluster_id,
            draft.members,
            draft.selection_type,
            draft.v1_event,
            draft.processing_level,
            draft.membership_hash,
        )
        groups.append(group)
        family_rows.append(
            dict(
                family_id=family_id,
                component_id=group.component_id,
                local_group_id=group.local_group_id,
                source_cluster_id=group.source_cluster_id,
                selection_type=group.selection_type,
                v1_event=group.v1_event,
                processing_level=group.processing_level,
                n_genes=len(group.members),
                n_species=len({p.species_id for p in group.members}),
                membership_hash=group.membership_hash,
            )
        )
        member_rows.extend(
            dict(
                family_id=family_id,
                protein_id=p.protein_id,
                species_id=p.species_id,
                original_id=p.original_id,
            )
            for p in group.members
        )
    counts = dict(
        orthogroups=len(groups),
        terminal_families=len(diagnostic.families),
        input_proteins=len(proteins),
        assigned_proteins=len(assigned),
        unassigned_proteins=len(unassigned_ids),
        duplicate_assignments=0,
    )
    statistics = tuple(dict(metric=key, value=value) for key, value in sorted(counts.items()))
    return OrthogroupExportTables(
        tuple(family_rows),
        tuple(member_rows),
        tuple(groups),
        tuple(sorted(unassigned, key=lambda row: row["protein_id"])),
        diagnostic,
        statistics,
    )
