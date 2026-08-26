"""Production component-local hierarchy stage with verified resume."""

from __future__ import annotations

import json
import os
import shutil
import uuid
from dataclasses import asdict
from pathlib import Path
from typing import Any

from ogprofiler.core.manifest import sha256_file, sha256_json, write_json
from ogprofiler.hierarchy.engine import HierarchyConfig, infer_component_hierarchy
from ogprofiler.hierarchy.loader import ComponentGraphLoader
from ogprofiler.hierarchy.validation import validate_hierarchy
from ogprofiler.storage.hierarchy import write_hierarchy_result

HIERARCHY_ALGORITHM_VERSION = "hierarchical-leiden-v1"
ARTIFACT_NAMES = ("nodes.parquet", "members.parquet", "candidates.parquet", "metrics.json")


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected object in {path}")
    return value


def hierarchy_component_identity(
    run_root: Path, component_id: int, config: HierarchyConfig
) -> tuple[str, dict[str, str], dict[str, Any]]:
    component_root = run_root / "components"
    partition_root = component_root / "edges" / f"component={component_id:08d}"
    input_paths = [component_root / "index.parquet", run_root / "input" / "proteins.parquet"]
    if partition_root.is_dir():
        input_paths.extend(sorted(partition_root.glob("*.parquet")))
    inputs = {path.relative_to(run_root).as_posix(): sha256_file(path) for path in input_paths}
    parameters = {"component_id": component_id, "hierarchy": asdict(config)}
    identity = sha256_json(
        {
            "algorithm_version": HIERARCHY_ALGORITHM_VERSION,
            "parameters": parameters,
            "inputs": inputs,
        }
    )
    return identity, inputs, parameters


def hierarchy_component_is_verified(
    run_root: Path, component_id: int, config: HierarchyConfig
) -> bool:
    _, inputs, parameters = hierarchy_component_identity(run_root, component_id, config)
    output = run_root / "hierarchy" / "components" / f"component={component_id:08d}"
    manifest_path = output / "hierarchy-manifest.json"
    if not manifest_path.is_file() or not all((output / name).is_file() for name in ARTIFACT_NAMES):
        return False
    try:
        previous = _load_json(manifest_path)
        checksums = {name: sha256_file(output / name) for name in ARTIFACT_NAMES}
        return bool(
            previous["algorithm_version"] == HIERARCHY_ALGORITHM_VERSION
            and previous["parameters"] == parameters
            and previous["input_checksums"] == inputs
            and previous["output_checksums"] == checksums
        )
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
        return False


def run_hierarchy_component_stage(
    run_root: Path,
    component_id: int,
    config: HierarchyConfig,
    command: list[str],
) -> tuple[Path, bool]:
    _, inputs, parameters = hierarchy_component_identity(run_root, component_id, config)
    output = run_root / "hierarchy" / "components" / f"component={component_id:08d}"
    if hierarchy_component_is_verified(run_root, component_id, config):
        return output, True
    loaded = ComponentGraphLoader(run_root).load(component_id)
    if config.subtree_workers > 1 and len(loaded.component.vertices) >= config.subtree_release_size:
        from ogprofiler.hierarchy.subtree import infer_component_hierarchy_parallel

        result = infer_component_hierarchy_parallel(
            loaded.component,
            config,
            loaded.species_by_protein,
            loaded.graph,
            workers=config.subtree_workers,
        )
    else:
        result = infer_component_hierarchy(
            loaded.component,
            config,
            loaded.species_by_protein,
            root_graph=loaded.graph,
        )
    validate_hierarchy(loaded.component, result)
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = output.parent / f".{output.name}.{uuid.uuid4().hex}.tmp"
    try:
        write_hierarchy_result(staging, result)
        checksums = {name: sha256_file(staging / name) for name in ARTIFACT_NAMES}
        write_json(
            staging / "hierarchy-manifest.json",
            {
                "algorithm_version": HIERARCHY_ALGORITHM_VERSION,
                "command": command,
                "parameters": parameters,
                "parameters_sha256": sha256_json(parameters),
                "input_checksums": inputs,
                "output_checksums": checksums,
                "metrics": asdict(result.metrics),
            },
        )
        output.mkdir(parents=True, exist_ok=True)
        for name in (*ARTIFACT_NAMES, "hierarchy-manifest.json"):
            os.replace(staging / name, output / name)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return output, False
