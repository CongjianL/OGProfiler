"""Final TSV and opt-in per-family FASTA export stage."""

from __future__ import annotations

import json
import os
import uuid
from pathlib import Path
from typing import Any

from ogprofiler.core.manifest import sha256_file, sha256_json, write_json
from ogprofiler.exceptions import ExportError, HierarchyError
from ogprofiler.input.fasta import parse_fasta
from ogprofiler.orthogroups.models import OrthogroupConfig
from ogprofiler.orthogroups.stage import verified_orthogroup_inputs
from ogprofiler.output.orthogroups import (
    OG_ID_ALGORITHM,
    ExportOrthogroup,
    assemble_orthogroup_export_tables,
)
from ogprofiler.output.results import (
    FAMILY_ID_ALGORITHM,
    TerminalFamily,
    _rows,
    assemble_export_tables,
    family_membership_sha256,
    write_tsv,
)

EXPORT_ALGORITHM_VERSION = "orthogroup-and-terminal-diagnostic-export-v2"
TERMINAL_TABLE_FIELDS = {
    "families.tsv": [
        "family_id",
        "component_id",
        "cluster_id",
        "n_genes",
        "n_species",
        "terminal_reason",
        "network_event",
    ],
    "members.tsv": ["family_id", "protein_id", "species_id", "original_id"],
    "hierarchy.tsv": [
        "cluster_id",
        "parent_id",
        "component_id",
        "depth",
        "n_genes",
        "n_species",
        "resolution",
        "quality",
        "child_count",
        "terminal_reason",
    ],
    "events.tsv": [
        "component_id",
        "cluster_id",
        "network_event",
        "overlap_score",
        "confidence",
    ],
}


TABLE_FIELDS = {
    **TERMINAL_TABLE_FIELDS,
    "families.tsv": [
        "family_id",
        "component_id",
        "local_group_id",
        "source_cluster_id",
        "selection_type",
        "v1_event",
        "processing_level",
        "n_genes",
        "n_species",
        "membership_hash",
    ],
    "terminal_families.tsv": TERMINAL_TABLE_FIELDS["families.tsv"],
    "terminal_members.tsv": TERMINAL_TABLE_FIELDS["members.tsv"],
    "unassigned.tsv": [
        "component_id",
        "protein_id",
        "species_id",
        "original_id",
        "terminal_cluster_id",
        "reason",
    ],
    "statistics.tsv": ["metric", "value"],
}


def _terminal_input_paths(run_root: Path) -> list[Path]:
    paths = [
        run_root / "input" / "proteins.parquet",
        run_root / "input" / "proteins.faa",
        run_root / "components" / "singleton_terminal_families.parquet",
    ]
    index_path = run_root / "components/index.parquet"
    current = None
    if index_path.is_file():
        paths.append(index_path)
        current = {int(row["component_id"]) for row in _rows(index_path)}
    for directory in sorted((run_root / "hierarchy" / "components").glob("component=*")):
        if current is not None and int(directory.name.split("=", 1)[1]) not in current:
            continue
        paths.extend((directory / "nodes.parquet", directory / "members.parquet"))
        paths.append(run_root / "evolution" / "components" / directory.name / "events.parquet")
    missing = [path for path in paths if not path.is_file()]
    if missing:
        raise ExportError(f"Missing export input: {missing[0]}")
    return paths


def _manifest_verified(
    manifest_path: Path,
    output_root: Path,
    inputs: dict[str, str],
    parameters: dict[str, Any],
    table_fields: dict[str, list[str]],
) -> bool:
    if not manifest_path.is_file():
        return False
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if (
            manifest["algorithm_version"] != EXPORT_ALGORITHM_VERSION
            or manifest["input_checksums"] != inputs
            or manifest["parameters"] != parameters
        ):
            return False
        outputs = manifest["output_checksums"]
        if not set(table_fields).issubset(outputs):
            return False
        if len(outputs) != len(table_fields) + int(manifest["counts"]["fasta_files"]):
            return False
        return all(
            (output_root / name).is_file() and sha256_file(output_root / name) == checksum
            for name, checksum in outputs.items()
        )
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
        return False


def _read_sequences(path: Path) -> dict[int, str]:
    result: dict[int, str] = {}
    for record in parse_fasta(path, "error"):
        if not record.identifier.startswith("OGP2P"):
            raise ExportError(f"Unexpected prepared FASTA identifier: {record.identifier}")
        try:
            protein_id = int(record.identifier[5:])
        except ValueError as error:
            raise ExportError(f"Invalid prepared FASTA identifier: {record.identifier}") from error
        result[protein_id] = record.sequence
    return result


def _write_family_fasta(
    run_root: Path,
    output_root: Path,
    groups: tuple[TerminalFamily, ...] | tuple[ExportOrthogroup, ...],
    selected: set[str],
) -> list[Path]:
    known = {family.family_id for family in groups}
    unknown = sorted(selected - known)
    if unknown:
        raise ExportError(f"Unknown family ID requested for FASTA export: {unknown[0]}")
    sequences = _read_sequences(run_root / "input" / "proteins.faa")
    output_root = output_root / "fasta"
    output_root.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    for family in groups:
        if family.family_id not in selected:
            continue
        path = output_root / f"{family.family_id}.faa"
        temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
        try:
            with temporary.open("w", encoding="utf-8", newline="\n") as handle:
                for member in family.members:
                    sequence = sequences.get(member.protein_id)
                    if sequence is None:
                        raise ExportError(
                            f"Missing prepared sequence for protein {member.protein_id}"
                        )
                    handle.write(
                        f">OGP2P{member.protein_id:012d} original_id={member.original_id} "
                        f"protein_id={member.protein_id} "
                        f"species_id={member.species_id}\n"
                    )
                    for start in range(0, len(sequence), 80):
                        handle.write(sequence[start : start + 80] + "\n")
            os.replace(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)
        outputs.append(path)
    return outputs


def run_export_stage(
    run_root: Path,
    command: list[str],
    *,
    fasta_families: tuple[str, ...] = (),
    all_family_fasta: bool = False,
    strategy: str = "v1_compatible",
    config: OrthogroupConfig | None = None,
) -> tuple[Path, bool, int]:
    """Default OG export; explicit terminal diagnostics use a separate namespace."""
    output_root = run_root / "results"
    if strategy == "v1_compatible":
        try:
            _, paths = verified_orthogroup_inputs(run_root, config)
        except HierarchyError as error:
            raise ExportError(str(error)) from error
        input_paths = [*paths, run_root / "input/proteins.faa"]
        fields = TABLE_FIELDS
        id_algorithm = OG_ID_ALGORITHM
    elif strategy == "terminal":
        input_paths = _terminal_input_paths(run_root)
        output_root = output_root / "terminal-families"
        fields = TERMINAL_TABLE_FIELDS
        id_algorithm = FAMILY_ID_ALGORITHM
    else:
        raise ExportError("export strategy must be v1_compatible or terminal")
    if any(not path.is_file() for path in input_paths):
        raise ExportError("Missing export input, including prepared proteins.faa")
    inputs = {path.relative_to(run_root).as_posix(): sha256_file(path) for path in input_paths}
    parameters = {
        "strategy": strategy,
        "family_id_algorithm": id_algorithm,
        "fasta_families": sorted(set(fasta_families)),
        "all_family_fasta": all_family_fasta,
    }
    output_root.mkdir(parents=True, exist_ok=True)
    manifest_path = output_root / "export-manifest.json"
    if _manifest_verified(manifest_path, output_root, inputs, parameters, fields):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        return manifest_path, True, int(manifest["counts"]["families"])
    previous_fasta: set[str] = set()
    if manifest_path.is_file():
        try:
            previous = json.loads(manifest_path.read_text())
            previous_fasta = {
                name
                for name in previous["output_checksums"]
                if Path(name).parent == Path("fasta") and Path(name).suffix == ".faa"
            }
        except (OSError, KeyError, ValueError, TypeError):
            pass
    groups: tuple[TerminalFamily, ...] | tuple[ExportOrthogroup, ...]
    if strategy == "v1_compatible":
        og_tables = assemble_orthogroup_export_tables(run_root, config)
        groups = og_tables.groups
        rows_by_name = {
            "families.tsv": og_tables.families,
            "members.tsv": og_tables.members,
            "hierarchy.tsv": og_tables.diagnostics.hierarchy,
            "events.tsv": og_tables.diagnostics.events,
            "terminal_families.tsv": og_tables.diagnostics.families,
            "terminal_members.tsv": og_tables.diagnostics.members,
            "unassigned.tsv": og_tables.unassigned,
            "statistics.tsv": og_tables.statistics,
        }
        extra_counts = {row["metric"]: row["value"] for row in og_tables.statistics}
    else:
        tables = assemble_export_tables(run_root)
        groups = tables.terminal_families
        rows_by_name = {
            "families.tsv": tables.families,
            "members.tsv": tables.members,
            "hierarchy.tsv": tables.hierarchy,
            "events.tsv": tables.events,
        }
        extra_counts = {}
    selected = set(fasta_families)
    if all_family_fasta:
        selected = {group.family_id for group in groups}
    if selected - {group.family_id for group in groups}:
        raise ExportError("Unknown family ID requested for FASTA export")
    # Completion is the final publication step; failed writes never leave a valid marker.
    manifest_path.unlink(missing_ok=True)
    for name, table_fields in fields.items():
        write_tsv(output_root / name, table_fields, rows_by_name[name])
    fasta_outputs = _write_family_fasta(run_root, output_root, groups, selected) if selected else []
    current_fasta = {path.relative_to(output_root).as_posix() for path in fasta_outputs}
    stale_fasta = previous_fasta - current_fasta
    if stale_fasta:
        archive = output_root / "previous-fasta" / uuid.uuid4().hex
        archive.mkdir(parents=True)
        for name in sorted(stale_fasta):
            if (output_root / name).is_file():
                os.replace(output_root / name, archive / Path(name).name)
    outputs = [output_root / name for name in fields] + fasta_outputs
    manifest = {
        "algorithm_version": EXPORT_ALGORITHM_VERSION,
        "command": command,
        "parameters": parameters,
        "parameters_sha256": sha256_json(parameters),
        "input_checksums": inputs,
        "output_checksums": {
            path.relative_to(output_root).as_posix(): sha256_file(path) for path in outputs
        },
        "counts": {
            **extra_counts,
            "families": len(rows_by_name["families.tsv"]),
            "members": len(rows_by_name["members.tsv"]),
            "hierarchy_nodes": len(rows_by_name["hierarchy.tsv"]),
            "events": len(rows_by_name["events.tsv"]),
            "fasta_files": len(fasta_outputs),
        },
        "family_membership_sha256": {
            group.family_id: family_membership_sha256(group) for group in groups
        },
    }
    temporary = manifest_path.with_name(f".export-manifest.{uuid.uuid4().hex}.tmp")
    try:
        write_json(temporary, manifest)
        os.replace(temporary, manifest_path)
    finally:
        temporary.unlink(missing_ok=True)
    return manifest_path, False, len(rows_by_name["families.tsv"])
