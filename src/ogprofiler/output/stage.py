"""Final TSV and opt-in per-family FASTA export stage."""

from __future__ import annotations

import json
import os
import uuid
from pathlib import Path
from typing import Any

from ogprofiler.core.manifest import sha256_file, sha256_json, write_json
from ogprofiler.exceptions import ExportError
from ogprofiler.input.fasta import parse_fasta
from ogprofiler.output.results import (
    FAMILY_ID_ALGORITHM,
    ExportTables,
    assemble_export_tables,
    family_membership_sha256,
    write_tsv,
)

EXPORT_ALGORITHM_VERSION = "terminal-family-export-v1"
TABLE_FIELDS = {
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


def _input_paths(run_root: Path) -> list[Path]:
    paths = [
        run_root / "input" / "proteins.parquet",
        run_root / "input" / "proteins.faa",
        run_root / "components" / "singleton_terminal_families.parquet",
    ]
    for directory in sorted((run_root / "hierarchy" / "components").glob("component=*")):
        paths.extend((directory / "nodes.parquet", directory / "members.parquet"))
        paths.append(
            run_root / "evolution" / "components" / directory.name / "events.parquet"
        )
    missing = [path for path in paths if not path.is_file()]
    if missing:
        raise ExportError(f"Missing export input: {missing[0]}")
    return paths


def _manifest_verified(
    manifest_path: Path,
    output_root: Path,
    inputs: dict[str, str],
    parameters: dict[str, Any],
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
        if not set(TABLE_FIELDS).issubset(outputs):
            return False
        if len(outputs) != len(TABLE_FIELDS) + int(manifest["counts"]["fasta_files"]):
            return False
        return all(
            (output_root / name).is_file()
            and sha256_file(output_root / name) == checksum
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


def _write_family_fasta(run_root: Path, tables: ExportTables, selected: set[str]) -> list[Path]:
    known = {family.family_id for family in tables.terminal_families}
    unknown = sorted(selected - known)
    if unknown:
        raise ExportError(f"Unknown family ID requested for FASTA export: {unknown[0]}")
    sequences = _read_sequences(run_root / "input" / "proteins.faa")
    output_root = run_root / "results" / "fasta"
    output_root.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    for family in tables.terminal_families:
        if family.family_id not in selected:
            continue
        path = output_root / f"{family.family_id}.faa"
        temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            for member in family.members:
                sequence = sequences.get(member.protein_id)
                if sequence is None:
                    raise ExportError(f"Missing prepared sequence for protein {member.protein_id}")
                handle.write(
                    f">{member.original_id} protein_id={member.protein_id} "
                    f"species_id={member.species_id}\n"
                )
                for start in range(0, len(sequence), 80):
                    handle.write(sequence[start : start + 80] + "\n")
        os.replace(temporary, path)
        outputs.append(path)
    return outputs


def run_export_stage(
    run_root: Path,
    command: list[str],
    *,
    fasta_families: tuple[str, ...] = (),
    all_family_fasta: bool = False,
) -> tuple[Path, bool, int]:
    input_paths = _input_paths(run_root)
    inputs = {path.relative_to(run_root).as_posix(): sha256_file(path) for path in input_paths}
    parameters = {
        "family_id_algorithm": FAMILY_ID_ALGORITHM,
        "fasta_families": sorted(set(fasta_families)),
        "all_family_fasta": all_family_fasta,
    }
    output_root = run_root / "results"
    output_root.mkdir(parents=True, exist_ok=True)
    manifest_path = output_root / "export-manifest.json"
    if _manifest_verified(manifest_path, output_root, inputs, parameters):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        return manifest_path, True, int(manifest["counts"]["families"])

    tables = assemble_export_tables(run_root)
    rows_by_name = {
        "families.tsv": tables.families,
        "members.tsv": tables.members,
        "hierarchy.tsv": tables.hierarchy,
        "events.tsv": tables.events,
    }
    for name, fields in TABLE_FIELDS.items():
        write_tsv(output_root / name, fields, rows_by_name[name])

    selected = set(fasta_families)
    if all_family_fasta:
        selected = {family.family_id for family in tables.terminal_families}
    fasta_outputs = _write_family_fasta(run_root, tables, selected) if selected else []
    outputs = [output_root / name for name in TABLE_FIELDS] + fasta_outputs
    write_json(
        manifest_path,
        {
            "algorithm_version": EXPORT_ALGORITHM_VERSION,
            "command": command,
            "parameters": parameters,
            "parameters_sha256": sha256_json(parameters),
            "input_checksums": inputs,
            "output_checksums": {
                path.relative_to(output_root).as_posix(): sha256_file(path) for path in outputs
            },
            "counts": {
                "families": len(tables.families),
                "members": len(tables.members),
                "hierarchy_nodes": len(tables.hierarchy),
                "events": len(tables.events),
                "fasta_files": len(fasta_outputs),
            },
            "family_membership_sha256": {
                family.family_id: family_membership_sha256(family)
                for family in tables.terminal_families
            },
        },
    )
    return manifest_path, False, len(tables.families)
