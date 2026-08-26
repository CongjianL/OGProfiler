"""Deterministic proteome discovery, ID assignment, and metadata output."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from ogprofiler.core.ids import sorted_input_paths, sorted_original_ids
from ogprofiler.core.manifest import checksum_manifest, combined_checksum
from ogprofiler.core.models import DatasetManifest, Protein, Species
from ogprofiler.exceptions import InputError
from ogprofiler.input.fasta import FastaRecord, parse_fasta


@dataclass(frozen=True, slots=True)
class PrepareResult:
    species: tuple[Species, ...]
    proteins: tuple[Protein, ...]
    dataset_manifest: DatasetManifest


def discover_proteomes(directory: Path, extensions: list[str]) -> list[Path]:
    directory = directory.expanduser().resolve()
    if not directory.is_dir():
        raise InputError(f"Proteome input is not a directory: {directory}")
    normalized_extensions = {extension.casefold() for extension in extensions}
    paths = [
        path
        for path in directory.iterdir()
        if path.is_file() and path.suffix.casefold() in normalized_extensions
    ]
    if not paths:
        expected = ", ".join(sorted(normalized_extensions))
        raise InputError(f"No proteome FASTA files in {directory}; expected: {expected}")
    return sorted_input_paths(paths)


def _species_name(path: Path) -> str:
    name = path.stem
    if not name:
        raise InputError(f"Could not derive species name from {path}")
    return name


def _write_parquet(output_dir: Path, species: list[Species], proteins: list[Protein]) -> None:
    species_table = pa.table(
        {
            "species_id": pa.array([item.species_id for item in species], type=pa.int32()),
            "species_name": pa.array([item.species_name for item in species], type=pa.string()),
            "source_file": pa.array([item.source_file for item in species], type=pa.string()),
        }
    )
    proteins_table = pa.table(
        {
            "protein_id": pa.array([item.protein_id for item in proteins], type=pa.int64()),
            "species_id": pa.array([item.species_id for item in proteins], type=pa.int32()),
            "original_id": pa.array([item.original_id for item in proteins], type=pa.string()),
            "length": pa.array([item.length for item in proteins], type=pa.int32()),
        }
    )
    pq.write_table(species_table, output_dir / "species.parquet", compression="zstd")
    pq.write_table(proteins_table, output_dir / "proteins.parquet", compression="zstd")


def prepare_proteomes(
    proteome_directory: Path,
    output_directory: Path,
    config: dict[str, Any],
) -> PrepareResult:
    paths = discover_proteomes(proteome_directory, config["input"]["extensions"])
    root = proteome_directory.expanduser().resolve()
    checksums = checksum_manifest(paths, root)
    policy = config["input"]["illegal_character_policy"]

    species: list[Species] = []
    records_by_species: list[dict[str, FastaRecord]] = []
    species_names: set[str] = set()
    for species_id, path in enumerate(paths):
        name = _species_name(path)
        if name in species_names:
            raise InputError(f"Duplicate species name derived from input files: {name}")
        species_names.add(name)
        records = {record.identifier: record for record in parse_fasta(path, policy)}
        species.append(Species(species_id, name, path.relative_to(root).as_posix()))
        records_by_species.append(records)

    proteins: list[Protein] = []
    ordered_records: list[tuple[Protein, FastaRecord]] = []
    for species_item, records in zip(species, records_by_species, strict=True):
        for original_id in sorted_original_ids(records):
            record = records[original_id]
            protein = Protein(
                protein_id=len(proteins),
                species_id=species_item.species_id,
                original_id=record.identifier,
                length=len(record.sequence),
            )
            proteins.append(protein)
            ordered_records.append((protein, record))

    output_directory.mkdir(parents=True, exist_ok=True)
    _write_parquet(output_directory, species, proteins)
    with (output_directory / "proteins.faa").open("w", encoding="utf-8", newline="\n") as handle:
        for protein, record in ordered_records:
            handle.write(
                f">OGP2P{protein.protein_id:012d} species_id={protein.species_id} "
                f"original_id={protein.original_id}\n"
            )
            for start in range(0, len(record.sequence), 80):
                handle.write(record.sequence[start : start + 80] + "\n")

    manifest = DatasetManifest(
        dataset_sha256=combined_checksum(checksums),
        input_checksums=checksums,
        n_species=len(species),
        n_proteins=len(proteins),
    )
    return PrepareResult(tuple(species), tuple(proteins), manifest)
