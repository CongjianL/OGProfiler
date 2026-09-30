"""Phase 4 edge-stage orchestration and verified resume."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq

from ogprofiler.core.manifest import sha256_file, sha256_json, write_json
from ogprofiler.exceptions import InputError
from ogprofiler.similarity.engine import (
    EdgeBuildConfig,
    build_retained_edges,
    filter_coverage,
)
from ogprofiler.similarity.io import (
    iter_directional_hits,
    write_normalized_hits,
    write_retained_edges,
)
from ogprofiler.similarity.normalization import normalize_hits

EDGE_ALGORITHM_VERSION = "of315-nbs-lrb-directional-mean-v4"


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object in {path}")
    return value


def run_edge_stage(
    run_root: Path, config: EdgeBuildConfig, command: list[str]
) -> tuple[Path, bool]:
    hits_path = run_root / "search" / "hits.parquet"
    proteins_path = run_root / "input" / "proteins.parquet"
    if not hits_path.is_file() or not proteins_path.is_file():
        raise InputError("Edge construction requires prepared proteins and search/hits.parquet")
    output_path = run_root / "edges" / "retained_edges.parquet"
    normalized_path = run_root / "edges" / "normalized_hits.parquet"
    manifest_path = run_root / "edges" / "edge-manifest.json"
    parameters = {
        "method": config.method,
        "normalization": config.normalization,
        "nbs_fallback": config.nbs_fallback,
        "apply_coverage_filter": config.apply_coverage_filter,
        "min_query_coverage": config.min_query_coverage,
        "min_target_coverage": config.min_target_coverage,
        "min_bidirectional_coverage": config.min_bidirectional_coverage,
        "best_hit_tolerance": config.best_hit_tolerance,
        "symmetrization": config.symmetrization,
    }
    inputs = {
        "hits.parquet": sha256_file(hits_path),
        "proteins.parquet": sha256_file(proteins_path),
    }
    if manifest_path.is_file() and output_path.is_file() and normalized_path.is_file():
        try:
            previous = _load_json(manifest_path)
            if bool(
                previous["algorithm_version"] == EDGE_ALGORITHM_VERSION
                and previous["parameters"] == parameters
                and previous["input_checksums"] == inputs
                and previous["output_checksums"]["retained_edges.parquet"]
                == sha256_file(output_path)
                and previous["output_checksums"]["normalized_hits.parquet"]
                == sha256_file(normalized_path)
            ):
                return output_path, True
        except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
            pass

    protein_rows = pq.read_table(proteins_path, columns=["protein_id", "length"]).to_pylist()
    lengths = {int(row["protein_id"]): int(row["length"]) for row in protein_rows}
    raw_hits = list(iter_directional_hits(hits_path))
    normalized = normalize_hits(raw_hits, lengths, config.normalization, config.nbs_fallback)
    coverage_filtered = (
        filter_coverage(normalized, config) if config.apply_coverage_filter else normalized
    )
    selected, edges = build_retained_edges(normalized, config)
    write_normalized_hits(normalized_path, normalized)
    write_retained_edges(output_path, edges)
    write_json(
        manifest_path,
        {
            "algorithm_version": EDGE_ALGORITHM_VERSION,
            "command": command,
            "parameters": parameters,
            "parameters_sha256": sha256_json(parameters),
            "input_checksums": inputs,
            "output_checksums": {
                "normalized_hits.parquet": sha256_file(normalized_path),
                "retained_edges.parquet": sha256_file(output_path),
            },
            "counts": {
                "directional_hits": len(raw_hits),
                "normalized_hits": len(normalized),
                "coverage_filtered_hits": len(coverage_filtered),
                "retained_directional_hits": len(selected),
                "retained_edges": len(edges),
            },
        },
    )
    return output_path, False
