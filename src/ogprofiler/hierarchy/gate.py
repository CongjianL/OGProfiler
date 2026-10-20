"""Shared upstream completion gate for OG generation and final consumers.

Structural leaves with search failures remain useful diagnostic artifacts, not OGs.
Legacy fixtures without new state columns are accepted only as legacy input.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pyarrow.parquet as pq

from ogprofiler.core.manifest import sha256_file
from ogprofiler.exceptions import HierarchyError


def require_resolved_hierarchy(run_root: Path) -> None:
    from ogprofiler.config import hierarchy_config, load_config
    from ogprofiler.hierarchy.stage import hierarchy_component_is_verified

    config = None
    config_path = run_root / "run.yaml"
    if config_path.is_file():
        config = hierarchy_config(load_config(str(config_path))["hierarchy"])
    rows = pq.read_table(
        run_root / "components/index.parquet", columns=["component_id"]
    ).to_pylist()
    sizes = Counter(int(row["component_id"]) for row in rows)
    for component_id in sorted(sizes):
        root = run_root / "hierarchy/components" / f"component={component_id:08d}"
        nodes_path = root / "nodes.parquet"
        if not nodes_path.is_file():
            # Singleton families are persisted separately and do not require a hierarchy.
            if sizes[component_id] == 1:
                continue
            raise HierarchyError(f"Missing hierarchy component {component_id}")
        nodes = pq.ParquetFile(nodes_path).read().to_pylist()
        if any(row.get("split_status") == "UNRESOLVED" for row in nodes):
            raise HierarchyError(
                f"Hierarchy component {component_id} is UNRESOLVED; OG publication stopped"
            )
        manifest_path = root / "hierarchy-manifest.json"
        new_state = any(row.get("search_status") is not None for row in nodes)
        requires_manifest = new_state or (
            config is not None and config.resolution.admission_policy == "nonempty_children_v1"
        )
        if not manifest_path.is_file():
            if requires_manifest:
                raise HierarchyError(f"Missing verified hierarchy manifest for {component_id}")
            continue
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("hierarchy_status") != "RESOLVED" or not manifest.get(
            "structural_validation_passed"
        ):
            raise HierarchyError(f"Hierarchy component {component_id} lacks resolved verification")
        for name, checksum in manifest["output_checksums"].items():
            if sha256_file(root / name) != checksum:
                raise HierarchyError(
                    f"Hierarchy component {component_id} artifact checksum mismatch"
                )
        for name, checksum in manifest["input_checksums"].items():
            if sha256_file(run_root / name) != checksum:
                raise HierarchyError(f"Hierarchy component {component_id} input checksum mismatch")
        # The manifest records effective stage overrides; run.yaml is the prepare baseline.
        from ogprofiler.hierarchy.engine import HierarchyConfig
        from ogprofiler.hierarchy.resolution import ResolutionSearchConfig

        recorded = dict(manifest["parameters"]["hierarchy"])
        recorded["resolution"] = ResolutionSearchConfig(**recorded["resolution"])
        effective = HierarchyConfig(**recorded)
        effective.resolution.validate()
        if not hierarchy_component_is_verified(run_root, component_id, effective):
            raise HierarchyError(f"Hierarchy component {component_id} is stale or corrupt")
