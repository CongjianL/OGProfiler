"""Search-stage orchestration, provenance, and verified resume.

The search stage runs per-species-pair searches (query species FASTA against
each species database), matching the V1/OrthoFinder architecture. A single
all-vs-all search against one concatenated database produces different hit
counts because DIAMOND's e-value depends on the database size.
"""

from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path
from typing import Any

from ogprofiler.core.manifest import sha256_file, sha256_json, write_json
from ogprofiler.exceptions import InputError
from ogprofiler.search.base import SearchBackend, SearchParameters, SearchStageResult
from ogprofiler.search.hits import parse_tabular_hits

SEARCH_ALGORITHM_VERSION = "directional-search-v4"


def create_backend(config: dict[str, Any]) -> SearchBackend:
    backend_name = str(config["search"]["backend"])
    if backend_name == "diamond":
        from ogprofiler.search.diamond import DiamondBackend

        return DiamondBackend(executable=str(config["search"]["executable"]))
    if backend_name == "mmseqs":
        from ogprofiler.search.mmseqs import MmseqsBackend

        return MmseqsBackend(
            executable=str(config["search"]["executable"]),
            sensitivity=str(config["search"]["mmseqs_sensitivity"]),
        )
    if backend_name == "blastp":
        from ogprofiler.search.blast import BlastBackend

        executable = str(config["search"]["executable"])
        makeblastdb = str(Path(executable).with_name("makeblastdb"))
        return BlastBackend(executable=executable, makeblastdb_executable=makeblastdb)
    raise InputError(f"Unknown search backend: {backend_name}")


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def _resume_matches(
    manifest_path: Path,
    hits_path: Path,
    *,
    backend: SearchBackend,
    backend_version: str,
    input_checksums: dict[str, str],
    parameters: SearchParameters,
) -> bool:
    if not manifest_path.is_file() or not hits_path.is_file():
        return False
    try:
        manifest = _load_json(manifest_path)
        return bool(
            manifest["algorithm_version"] == SEARCH_ALGORITHM_VERSION
            and manifest["backend"] == backend.name
            and manifest["backend_version"] == backend_version
            and manifest["parameters"] == parameters.to_dict()
            and manifest["input_checksums"] == input_checksums
            and manifest["output_checksums"]["hits.parquet"] == sha256_file(hits_path)
        )
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
        return False


def _write_species_fasta(input_fasta: Path, species_dir: Path) -> tuple[int, ...]:
    """Split the prepared proteins.faa into one FASTA file per species.

    Each written header keeps only the ``OGP2P...`` identifier so DIAMOND emits
    the same query/target IDs the hit parser expects.
    """
    species_dir.mkdir(parents=True, exist_ok=True)
    writers: dict[int, Any] = {}
    species_ids: list[int] = []
    current_species: int | None = None
    with input_fasta.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                species_id = int(line.split("species_id=")[1].split()[0])
                current_species = species_id
                if species_id not in writers:
                    path = species_dir / f"species_{species_id}.faa"
                    writers[species_id] = path.open("w", encoding="utf-8")
                    species_ids.append(species_id)
                writers[species_id].write(line.split()[0] + "\n")
            elif current_species is not None:
                writers[current_species].write(line)
    for writer in writers.values():
        writer.close()
    return tuple(sorted(species_ids))


def _merge_raw_files(paths: list[Path], output: Path) -> None:
    with output.open("w", encoding="utf-8") as handle:
        for path in paths:
            if not path.is_file():
                continue
            with path.open(encoding="utf-8") as source:
                handle.write(source.read())


def run_search_stage(
    run_root: Path,
    backend: SearchBackend,
    parameters: SearchParameters,
    command: list[str],
) -> SearchStageResult:
    input_fasta = run_root / "input" / "proteins.faa"
    proteins_path = run_root / "input" / "proteins.parquet"
    if not input_fasta.is_file() or not proteins_path.is_file():
        raise InputError(f"Prepared input is missing under run workspace: {run_root / 'input'}")

    search_root = run_root / "search"
    backend_root = search_root / backend.name
    hits_path = search_root / "hits.parquet"
    manifest_path = search_root / "search-manifest.json"
    input_checksums = {
        "proteins.faa": sha256_file(input_fasta),
        "proteins.parquet": sha256_file(proteins_path),
    }
    version = backend.version()
    if _resume_matches(
        manifest_path,
        hits_path,
        backend=backend,
        backend_version=version,
        input_checksums=input_checksums,
        parameters=parameters,
    ):
        return SearchStageResult(hits_path, manifest_path, reused=True)

    backend_root.mkdir(parents=True, exist_ok=True)
    species_dir = backend_root / "species_fasta"
    database_dir = backend_root / "databases"
    raw_dir = backend_root / "raw"
    species_ids = _write_species_fasta(input_fasta, species_dir)

    # Build one database per species.
    database_commands: list[tuple[str, ...]] = []
    database_files: list[Path] = []
    for species_id in species_ids:
        species_fasta = species_dir / f"species_{species_id}.faa"
        database_path = database_dir / f"species_{species_id}"
        database_commands.append(backend.build_database(species_fasta, database_path))
        database_files.extend(backend.database_files(database_path))

    # Run one search per (query species, target species) pair, in parallel.
    # Each per-pair search is single-threaded; concurrency comes from the pool.
    single_thread = replace(parameters, threads=1)
    search_commands: list[tuple[str, ...]] = []
    raw_paths: list[Path] = []

    def _search_pair(query_species: int, target_species: int) -> tuple[str, ...]:
        query_path = species_dir / f"species_{query_species}.faa"
        database_path = database_dir / f"species_{target_species}"
        raw_path = raw_dir / f"raw_{query_species}_{target_species}.tsv"
        return backend.search(query_path, database_path, raw_path, single_thread)

    pairs = [
        (query_species, target_species)
        for query_species in species_ids
        for target_species in species_ids
    ]
    with ThreadPoolExecutor(max_workers=max(1, parameters.threads)) as pool:
        futures = [
            pool.submit(_search_pair, query_species, target_species)
            for query_species, target_species in pairs
        ]
        for (query_species, target_species), future in zip(pairs, futures, strict=True):
            search_commands.append(future.result())
            raw_paths.append(raw_dir / f"raw_{query_species}_{target_species}.tsv")

    raw_path = backend_root / "hits.tsv"
    _merge_raw_files(raw_paths, raw_path)
    hit_count = parse_tabular_hits(
        raw_path, proteins_path, hits_path, backend_name=backend.name
    )
    if not database_files:
        raise InputError(f"Search backend did not publish database files: {backend.name}")

    manifest = {
        "algorithm_version": SEARCH_ALGORITHM_VERSION,
        "backend": backend.name,
        "backend_version": version,
        "command": command,
        "search_mode": "per_species_pair",
        "species_ids": list(species_ids),
        "n_pairs": len(pairs),
        "database_commands": [list(item) for item in database_commands],
        "search_command": list(search_commands[0]) if search_commands else [],
        "parameters": parameters.to_dict(),
        "parameters_sha256": sha256_json(parameters.to_dict()),
        "input_checksums": input_checksums,
        "output_checksums": {
            **{
                f"database/{path.name}": sha256_file(path)
                for path in database_files
            },
            "hits.tsv": sha256_file(raw_path),
            "hits.parquet": sha256_file(hits_path),
        },
        "hit_count": hit_count,
        "directional": True,
        "hit_schema": [
            "query_id",
            "target_id",
            "query_species",
            "target_species",
            "bitscore",
            "identity",
            "query_coverage",
            "target_coverage",
            "evalue",
        ],
    }
    write_json(manifest_path, manifest)
    return SearchStageResult(hits_path, manifest_path, reused=False)
