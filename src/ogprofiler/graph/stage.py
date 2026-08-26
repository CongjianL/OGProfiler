"""Phase 5 connected-component stage with verified artifact reuse."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ogprofiler.core.manifest import checksum_manifest, sha256_file, sha256_json, write_json
from ogprofiler.exceptions import ComponentError
from ogprofiler.graph.partition import build_component_artifacts

COMPONENT_ALGORITHM_VERSION = "union-find-partition-v1"


@dataclass(frozen=True, slots=True)
class ComponentStageConfig:
    edge_batch_size: int = 65_536
    max_open_files: int = 64


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object in {path}")
    return value


def _artifact_paths(output_directory: Path) -> list[Path]:
    paths = [
        output_directory / "index.parquet",
        output_directory / "statistics.parquet",
        output_directory / "singleton_terminal_families.parquet",
    ]
    edge_root = output_directory / "edges"
    if edge_root.is_dir():
        paths.extend(sorted(path for path in edge_root.rglob("*.parquet") if path.is_file()))
    return paths


def run_component_stage(
    run_root: Path, config: ComponentStageConfig, command: list[str]
) -> tuple[Path, bool]:
    proteins_path = run_root / "input" / "proteins.parquet"
    edges_path = run_root / "edges" / "retained_edges.parquet"
    if not proteins_path.is_file() or not edges_path.is_file():
        raise ComponentError("Component construction requires proteins and retained edges")
    output_directory = run_root / "components"
    output_directory.mkdir(parents=True, exist_ok=True)
    manifest_path = output_directory / "component-manifest.json"
    parameters = {
        "edge_batch_size": config.edge_batch_size,
        "max_open_files": config.max_open_files,
        "component_order": "n_vertices_desc_then_min_protein_id",
    }
    inputs = {
        "proteins.parquet": sha256_file(proteins_path),
        "retained_edges.parquet": sha256_file(edges_path),
    }
    if manifest_path.is_file():
        try:
            previous = _load_json(manifest_path)
            artifacts = _artifact_paths(output_directory)
            checksums = checksum_manifest(artifacts, output_directory)
            if bool(
                previous["algorithm_version"] == COMPONENT_ALGORITHM_VERSION
                and previous["parameters"] == parameters
                and previous["input_checksums"] == inputs
                and previous["output_checksums"] == checksums
            ):
                return output_directory / "index.parquet", True
        except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
            pass

    result = build_component_artifacts(
        proteins_path,
        edges_path,
        output_directory,
        batch_size=config.edge_batch_size,
        max_open_files=config.max_open_files,
    )
    artifacts = _artifact_paths(output_directory)
    write_json(
        manifest_path,
        {
            "algorithm_version": COMPONENT_ALGORITHM_VERSION,
            "command": command,
            "parameters": parameters,
            "parameters_sha256": sha256_json(parameters),
            "input_checksums": inputs,
            "output_checksums": checksum_manifest(artifacts, output_directory),
            "counts": {
                "proteins": result.protein_count,
                "retained_edges": result.edge_count,
                "components": result.component_count,
                "singletons": result.singleton_count,
                "partition_files": sum(
                    1 for path in (output_directory / "edges").rglob("*.parquet")
                ),
            },
        },
    )
    return output_directory / "index.parquet", False
