"""Search-stage orchestration, provenance, and verified resume."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ogprofiler.core.manifest import sha256_file, sha256_json, write_json
from ogprofiler.exceptions import InputError
from ogprofiler.search.base import SearchBackend, SearchParameters, SearchStageResult
from ogprofiler.search.hits import parse_diamond_hits

SEARCH_ALGORITHM_VERSION = "directional-search-v1"


def create_backend(config: dict[str, Any]) -> SearchBackend:
    backend_name = str(config["search"]["backend"])
    if backend_name == "diamond":
        from ogprofiler.search.diamond import DiamondBackend

        return DiamondBackend(executable=str(config["search"]["executable"]))
    raise InputError(f"Search backend is configured but not implemented in Phase 3: {backend_name}")


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
    database_path = backend_root / "proteins"
    raw_path = backend_root / "hits.tsv"
    database_command = backend.build_database(input_fasta, database_path)
    search_command = backend.search(input_fasta, database_path, raw_path, parameters)
    hit_count = parse_diamond_hits(raw_path, proteins_path, hits_path)
    database_file = database_path.with_suffix(".dmnd")
    manifest = {
        "algorithm_version": SEARCH_ALGORITHM_VERSION,
        "backend": backend.name,
        "backend_version": version,
        "command": command,
        "database_command": list(database_command),
        "search_command": list(search_command),
        "parameters": parameters.to_dict(),
        "parameters_sha256": sha256_json(parameters.to_dict()),
        "input_checksums": input_checksums,
        "output_checksums": {
            "database.dmnd": sha256_file(database_file),
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
