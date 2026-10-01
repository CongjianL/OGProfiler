"""Component OG artifacts, manifest-last publication and checksum-verified resume.

A single orchestrator owns run.db. Workers operate on distinct components;
concurrent orchestrators in the same workspace are not supported.
"""

from __future__ import annotations

import json
import multiprocessing
import os
import tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from ogprofiler.core.checkpoint import CheckpointStore
from ogprofiler.core.manifest import sha256_file, sha256_json, write_json
from ogprofiler.core.workspace import Workspace
from ogprofiler.exceptions import HierarchyError, InputError
from ogprofiler.orthogroups.engine import (
    OG_EXTRACTION_ALGORITHM_VERSION,
    OrthogroupConflictError,
    extract_component_orthogroups,
)
from ogprofiler.orthogroups.legacy_events import V1_EVENT_ALGORITHM_VERSION, annotate_v1_events
from ogprofiler.orthogroups.models import OrthogroupConfig as OrthogroupConfig

OG_STAGE_VERSION = "component-orthogroups-parquet-v1"
SCHEMA_VERSION = 1
SCHEMAS = {
    "groups.parquet": pa.schema(
        [
            ("component_id", pa.int64()),
            ("local_group_id", pa.int64()),
            ("source_cluster_id", pa.int64()),
            ("selection_type", pa.string()),
            ("v1_event", pa.string()),
            ("processing_level", pa.int32()),
            ("n_genes", pa.int64()),
            ("n_species", pa.int32()),
            ("membership_hash", pa.string()),
        ]
    ),
    "members.parquet": pa.schema(
        [
            ("component_id", pa.int64()),
            ("local_group_id", pa.int64()),
            ("protein_id", pa.int64()),
        ]
    ),
    "selection_trace.parquet": pa.schema(
        [
            ("component_id", pa.int64()),
            ("trace_order", pa.int64()),
            ("cluster_id", pa.int64()),
            ("processing_level", pa.int32()),
            ("selection_event", pa.string()),
            ("status", pa.string()),
            ("consumed_by", pa.int64()),
        ]
    ),
    "unassigned.parquet": pa.schema(
        [
            ("component_id", pa.int64()),
            ("protein_id", pa.int64()),
            ("terminal_cluster_id", pa.int64()),
            ("reason", pa.string()),
        ]
    ),
    "v1_events.parquet": pa.schema(
        [
            ("component_id", pa.int64()),
            ("cluster_id", pa.int64()),
            ("reference_order", pa.int64()),
            ("undirected_degree", pa.int32()),
            ("eligible_child_count", pa.int32()),
            ("v1_event", pa.string()),
        ]
    ),
}


def _atomic_json(path: Path, value: Any) -> None:
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        write_json(temporary, value)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _output(run_root: Path, component_id: int) -> Path:
    return run_root / "orthogroups/components" / f"component={component_id:08d}"


def _identity(
    run_root: Path, component_id: int, config: OrthogroupConfig, shared: dict[str, str]
) -> dict[str, Any]:
    hierarchy = run_root / "hierarchy/components" / f"component={component_id:08d}"
    inputs = dict(shared)
    if hierarchy.is_dir():
        paths = [hierarchy / name for name in ("nodes.parquet", "members.parquet")]
        paths.append(run_root / "evolution/components" / hierarchy.name / "events.parquet")
        for path in paths:
            inputs[path.relative_to(run_root).as_posix()] = (
                sha256_file(path) if path.is_file() else "MISSING"
            )
    return {
        "algorithm_version": OG_STAGE_VERSION,
        "schema_version": SCHEMA_VERSION,
        "engine_version": OG_EXTRACTION_ALGORITHM_VERSION,
        "event_algorithm_version": V1_EVENT_ALGORITHM_VERSION,
        "parameters": {"component_id": component_id, **asdict(config)},
        "input_checksums": inputs,
    }


def _verified(output: Path, identity: dict[str, Any]) -> bool:
    try:
        manifest = json.loads((output / "og-manifest.json").read_text())
        if manifest["status"] != "DONE" or any(
            manifest[key] != value for key, value in identity.items()
        ):
            return False
        if set(manifest["output_checksums"]) != set(SCHEMAS):
            return False
        for name, schema in SCHEMAS.items():
            path = output / name
            if sha256_file(path) != manifest["output_checksums"][name]:
                return False
            if not pq.ParquetFile(path).schema_arrow.equals(schema):
                return False
        return True
    except (OSError, KeyError, ValueError, TypeError, pa.ArrowException):
        return False


def _write_rows(path: Path, schema: pa.Schema, rows: Any) -> None:
    """Bound serialization memory for large membership/trace outputs."""
    with pq.ParquetWriter(path, schema, compression="zstd") as writer:
        batch = []
        for row in rows:
            batch.append(row)
            if len(batch) == 65_536:
                writer.write_table(pa.Table.from_pylist(batch, schema=schema))
                batch.clear()
        if batch:
            writer.write_table(pa.Table.from_pylist(batch, schema=schema))


def _build_component(
    run_root: Path,
    component_id: int,
    config: OrthogroupConfig,
    command: list[str],
    identity: dict[str, Any],
    total_species: int,
) -> Path:
    output = _output(run_root, component_id)
    output.mkdir(parents=True, exist_ok=True)
    # Remove the completion marker before any possible failed rebuild.
    (output / "og-manifest.json").unlink(missing_ok=True)
    for name, checksum in identity["input_checksums"].items():
        if checksum == "MISSING":
            raise HierarchyError(f"Missing OG input: {run_root / name}")
    index = pq.read_table(
        run_root / "components/index.parquet", filters=[("component_id", "=", component_id)]
    ).to_pylist()
    protein_ids = [
        int(row["protein_id"]) for row in index if int(row["component_id"]) == component_id
    ]
    if not protein_ids or len(set(protein_ids)) != len(protein_ids):
        raise HierarchyError(f"Invalid component index for {component_id}")
    proteins = pq.read_table(
        run_root / "input/proteins.parquet",
        columns=["protein_id", "species_id", "original_id"],
        filters=[("protein_id", "in", protein_ids)],
    ).to_pylist()
    species = {int(row["protein_id"]): int(row["species_id"]) for row in proteins}
    original = {int(row["protein_id"]): str(row["original_id"]) for row in proteins}
    if len(proteins) != len(species) or set(species) != set(protein_ids):
        raise HierarchyError("OG protein metadata differs from component index")
    isolates = pq.read_table(
        run_root / "components/singleton_terminal_families.parquet",
        filters=[("component_id", "=", component_id)],
    ).to_pylist()
    isolate_ids = [int(row["protein_id"]) for row in isolates]
    if isolate_ids and (len(protein_ids) != 1 or isolate_ids != protein_ids):
        raise HierarchyError("SSN isolate differs from component index")
    hierarchy = run_root / "hierarchy/components" / output.name
    if hierarchy.is_dir():
        nodes = pq.ParquetFile(hierarchy / "nodes.parquet").read().to_pylist()
        members = pq.ParquetFile(hierarchy / "members.parquet").read().to_pylist()
    elif isolate_ids:
        # Ephemeral view only: the hierarchy store is never modified.
        nodes = [
            dict(
                component_id=component_id,
                cluster_id=0,
                parent_id=None,
                depth=0,
                n_genes=1,
                n_species=1,
            )
        ]
        members = [dict(protein_id=protein_ids[0], terminal_cluster_id=0)]
    else:
        raise HierarchyError(f"Missing hierarchy for non-isolate component {component_id}")
    if {int(row["protein_id"]) for row in members} != set(protein_ids):
        raise HierarchyError("Terminal membership differs from component index")
    result = extract_component_orthogroups(
        nodes,
        members,
        species,
        original,
        total_species=total_species,
        overlap_count=config.species_overlap_count,
        ssn_isolates=isolate_ids,
    )
    if result.component_id != component_id:
        raise HierarchyError("Hierarchy component ID differs from component index")
    annotations = annotate_v1_events(nodes, members, species, config.species_overlap_count)
    rows = {
        "groups.parquet": (
            {key: getattr(group, key) for key in SCHEMAS["groups.parquet"].names}
            for group in result.groups
        ),
        "members.parquet": (
            dict(component_id=component_id, local_group_id=group.local_group_id, protein_id=protein)
            for group in result.groups
            for protein in group.protein_ids
        ),
        "selection_trace.parquet": (
            dict(component_id=component_id, trace_order=i, **asdict(trace))
            for i, trace in enumerate(result.trace)
        ),
        "unassigned.parquet": (
            dict(component_id=component_id, **asdict(item)) for item in result.unassigned
        ),
        "v1_events.parquet": (asdict(item) for item in annotations),
    }
    with tempfile.TemporaryDirectory(prefix=".og-build-", dir=output) as temporary:
        build = Path(temporary)
        for name, schema in SCHEMAS.items():
            _write_rows(build / name, schema, rows[name])
        # Manifest-last: a crash midway leaves no reusable completion marker.
        for name in SCHEMAS:
            os.replace(build / name, output / name)
    _atomic_json(
        output / "og-manifest.json",
        {
            **identity,
            "status": "DONE",
            "command": command,
            "output_checksums": {name: sha256_file(output / name) for name in SCHEMAS},
            "counts": {
                "selected": len(result.groups),
                "assigned": sum(g.n_genes for g in result.groups),
                "unassigned": len(result.unassigned),
                "duplicate": 0,
            },
            "remaining_cluster_ids": result.remaining_cluster_ids,
        },
    )
    (output / "og-failure.json").unlink(missing_ok=True)
    return output


def _execute_component(
    run_root: Path,
    component_id: int,
    config: OrthogroupConfig,
    command: list[str],
    identity: dict[str, Any],
    total_species: int,
) -> Path:
    try:
        return _build_component(run_root, component_id, config, command, identity, total_species)
    except Exception as error:  # publication/engine failure boundary
        output = _output(run_root, component_id)
        output.mkdir(parents=True, exist_ok=True)
        (output / "og-manifest.json").unlink(missing_ok=True)
        failure: dict[str, Any] = {
            **identity,
            "status": "FAILED",
            "error_type": type(error).__name__,
            "error": str(error),
        }
        if isinstance(error, OrthogroupConflictError):
            failure["duplicate_proteins"] = error.duplicate_members
            failure["counts"] = {
                "duplicate": len(error.duplicate_members),
                "selected_candidates": len(error.result.groups),
                "unassigned": len(error.result.unassigned),
            }
        _atomic_json(output / "og-failure.json", failure)
        raise


def _run_orthogroup_stage(
    run_root: Path,
    config: OrthogroupConfig,
    command: list[str],
    *,
    workers: int = 1,
    retries: int = 1,
) -> tuple[Path, int, int]:
    """Run all indexed components; recover stale RUNNING and retry failed components.

    Resume is governed by both SQLite identity and verified complete artifacts.
    Runtime concurrency/retry changes do not invalidate scientific artifacts.
    """
    for key, value, minimum in (("workers", workers, 1), ("retries", retries, 0)):
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise InputError(f"orthogroup {key} must be an integer >= {minimum}")
    workspace = Workspace.create(run_root)
    store = CheckpointStore(workspace.database_path)
    root = run_root / "orthogroups"
    root.mkdir(exist_ok=True)
    manifest_path = root / "og-manifest.json"
    manifest_path.unlink(missing_ok=True)
    shared_paths = [
        run_root / name
        for name in (
            "input/proteins.parquet",
            "input/species.parquet",
            "components/index.parquet",
            "components/singleton_terminal_families.parquet",
        )
    ]
    for path in shared_paths:
        if not path.is_file():
            raise HierarchyError(f"Missing OG input: {path}")
    shared = {path.relative_to(run_root).as_posix(): sha256_file(path) for path in shared_paths}
    species_rows = (
        pq.ParquetFile(run_root / "input/species.parquet").read(columns=["species_id"]).to_pylist()
    )
    total_species = len({int(row["species_id"]) for row in species_rows})
    if not total_species or total_species != len(species_rows):
        raise HierarchyError("Dataset species table must contain distinct species")
    index = pq.ParquetFile(run_root / "components/index.parquet").read().to_pylist()
    proteins = (
        pq.ParquetFile(run_root / "input/proteins.parquet")
        .read(columns=["protein_id", "species_id", "original_id"])
        .to_pylist()
    )
    if (
        len({int(row["protein_id"]) for row in proteins}) != len(proteins)
        or {int(row["protein_id"]) for row in proteins} != {int(row["protein_id"]) for row in index}
        or not {int(row["species_id"]) for row in proteins}.issubset(
            {int(row["species_id"]) for row in species_rows}
        )
    ):
        raise HierarchyError("Protein metadata, species table and component index disagree")
    if any(
        not isinstance(row["original_id"], str) or not row["original_id"] for row in proteins
    ) or len({(row["species_id"], row["original_id"]) for row in proteins}) != len(proteins):
        raise HierarchyError(
            "Protein original identities must be non-empty and unique within species"
        )
    if len({int(row["protein_id"]) for row in index}) != len(index):
        raise HierarchyError("Duplicate proteins in component index")
    component_ids = sorted({int(row["component_id"]) for row in index})
    if not component_ids:
        raise HierarchyError("OG extraction requires indexed components")
    component_by_protein = {int(row["protein_id"]): int(row["component_id"]) for row in index}
    sizes = Counter(component_by_protein.values())
    isolate_rows = (
        pq.ParquetFile(run_root / "components/singleton_terminal_families.parquet")
        .read(columns=["protein_id", "component_id"])
        .to_pylist()
    )
    if len({row["protein_id"] for row in isolate_rows}) != len(isolate_rows) or any(
        component_by_protein.get(int(row["protein_id"])) != int(row["component_id"])
        or sizes[int(row["component_id"])] != 1
        for row in isolate_rows
    ):
        raise HierarchyError("Singleton table differs from component index")
    pending = []
    identities = {}
    reused = 0
    failures = {}
    for component in component_ids:
        try:
            identity = _identity(run_root, component, config, shared)
            identities[component] = identity
            record = store.register(
                "orthogroups", str(component), sha256_json(identity), OG_STAGE_VERSION
            )
            if _verified(_output(run_root, component), identity):
                if record.status != "DONE":
                    store.finish("orthogroups", str(component), str(_output(run_root, component)))
                reused += 1
            else:
                if record.status == "DONE":
                    store.invalidate("orthogroups", str(component), "artifact verification failed")
                pending.append(component)
        except (OSError, HierarchyError, pa.ArrowException) as error:
            failures[component] = str(error)
    pending.sort(key=lambda component: (-sizes[component], component))
    remaining = pending
    for attempt in range(retries + 1):
        if not remaining:
            break
        current, remaining = remaining, []
        for component in current:
            store.start("orthogroups", str(component))
        outcomes: list[tuple[int, str | None]] = []
        if workers == 1:
            for component in current:
                try:
                    _execute_component(
                        run_root, component, config, command, identities[component], total_species
                    )
                    outcomes.append((component, None))
                except Exception as error:  # component failure boundary
                    outcomes.append((component, str(error)))
        else:
            with ProcessPoolExecutor(
                max_workers=workers, mp_context=multiprocessing.get_context("spawn")
            ) as executor:
                futures = {
                    executor.submit(
                        _execute_component,
                        run_root,
                        component,
                        config,
                        command,
                        identities[component],
                        total_species,
                    ): component
                    for component in current
                }
                for future in as_completed(futures):
                    component = futures[future]
                    try:
                        future.result()
                        outcomes.append((component, None))
                    except Exception as error:  # worker failure boundary
                        outcomes.append((component, str(error)))
        for component, failure_detail in outcomes:
            if failure_detail is None:
                store.finish("orthogroups", str(component), str(_output(run_root, component)))
            else:
                store.fail("orthogroups", str(component), failure_detail)
                if attempt < retries:
                    remaining.append(component)
                else:
                    failures[component] = failure_detail
    component_manifests = {
        str(component): sha256_file(_output(run_root, component) / "og-manifest.json")
        for component in component_ids
        if component not in failures
    }
    _atomic_json(
        manifest_path,
        {
            "algorithm_version": OG_STAGE_VERSION,
            "status": "FAILED" if failures else "DONE",
            "command": command,
            "parameters": asdict(config),
            "parameters_sha256": sha256_json(asdict(config)),
            "runtime": {"workers": workers, "retries": retries},
            "component_ids": component_ids,
            "component_manifests": component_manifests,
            "reused_components": reused,
            "failed_components": failures,
        },
    )
    if failures:
        raise HierarchyError(
            f"OG extraction failed for components {sorted(failures)}; inspect {manifest_path}"
        )
    return manifest_path, len(component_ids), reused


def run_orthogroup_stage(
    run_root: Path,
    config: OrthogroupConfig,
    command: list[str],
    *,
    workers: int = 1,
    retries: int = 1,
) -> tuple[Path, int, int]:
    """Public stage boundary: malformed storage inputs become actionable CLI errors."""
    try:
        return _run_orthogroup_stage(run_root, config, command, workers=workers, retries=retries)
    except (OSError, ValueError, KeyError, TypeError, pa.ArrowException) as error:
        raise HierarchyError(f"OG stage storage/input error: {error}") from error


def verified_orthogroup_inputs(
    run_root: Path,
    expected_config: OrthogroupConfig | None = None,
) -> tuple[tuple[int, ...], tuple[Path, ...]]:
    """Read-only consumer gate: check current source identity and complete artifacts."""
    try:
        root_manifest = run_root / "orthogroups/og-manifest.json"
        manifest = json.loads(root_manifest.read_text())
        if manifest["status"] != "DONE" or manifest["algorithm_version"] != OG_STAGE_VERSION:
            raise HierarchyError("OG stage is not complete; run orthogroups before export")
        config = OrthogroupConfig(**manifest["parameters"])
        if expected_config is not None and config != expected_config:
            raise HierarchyError("OG configuration changed; rerun orthogroups before export")
        shared_names = (
            "input/proteins.parquet",
            "input/species.parquet",
            "components/index.parquet",
            "components/singleton_terminal_families.parquet",
        )
        shared = {name: sha256_file(run_root / name) for name in shared_names}
        rows = pq.ParquetFile(run_root / "components/index.parquet").read().to_pylist()
        components = tuple(sorted({int(row["component_id"]) for row in rows}))
        if (
            not components
            or list(components) != manifest["component_ids"]
            or set(manifest["component_manifests"]) != {str(c) for c in components}
        ):
            raise HierarchyError("OG component inventory changed; rerun orthogroups")
        paths = {root_manifest, *(run_root / name for name in shared_names)}
        for component in components:
            output = _output(run_root, component)
            identity = _identity(run_root, component, config, shared)
            component_manifest = output / "og-manifest.json"
            if sha256_file(component_manifest) != manifest["component_manifests"][
                str(component)
            ] or not _verified(output, identity):
                raise HierarchyError(
                    f"OG component {component} is stale or corrupt; rerun orthogroups"
                )
            paths.add(component_manifest)
            paths.update(output / name for name in SCHEMAS)
            paths.update(run_root / name for name in identity["input_checksums"])
        return components, tuple(sorted(paths))
    except (OSError, ValueError, TypeError, KeyError, pa.ArrowException) as error:
        raise HierarchyError(f"Invalid OG artifacts; rerun orthogroups: {error}") from error
