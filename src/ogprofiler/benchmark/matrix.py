"""Deterministic one-factor-at-a-time scientific parameter matrix."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ogprofiler.config import DEFAULT_CONFIG
from ogprofiler.core.manifest import sha256_json, write_json


@dataclass(frozen=True, slots=True)
class MatrixRun:
    run_id: str
    axis: str
    value: Any
    overrides: tuple[str, ...]
    config_sha256: str


DEFAULT_AXES: dict[str, list[dict[str, Any]]] = {
    "normalization": [
        {"similarity.normalization": "legacy_nbs"},
        {"similarity.normalization": "raw_bitscore"},
        {"similarity.normalization": "length_scaled_bitscore"},
    ],
    "coverage": [
        {"edges.min_query_coverage": 0.0, "edges.min_target_coverage": 0.0},
        {"edges.min_query_coverage": 50.0, "edges.min_target_coverage": 50.0},
        {"edges.min_query_coverage": 70.0, "edges.min_target_coverage": 70.0},
    ],
    "symmetrization": [
        {"edges.symmetrization": value}
        for value in ("max", "min", "mean", "geometric_mean")
    ],
    "leiden_method": [
        {"hierarchy.method": value} for value in ("rber", "rbcv", "cpm", "modularity")
    ],
    "gamma_strategy": [
        {"hierarchy.resolution_strategy": "adaptive"},
        {"hierarchy.resolution_strategy": "log_grid"},
    ],
    "split_acceptance": [
        {"hierarchy.min_family_size": 2, "hierarchy.max_child_fraction": 0.95},
        {"hierarchy.min_family_size": 3, "hierarchy.max_child_fraction": 0.90},
        {"hierarchy.min_family_size": 5, "hierarchy.max_child_fraction": 0.80},
    ],
    "seed": [{"hierarchy.seed": value} for value in (7, 42, 104729)],
}


def _default_value(dotted: str) -> Any:
    section, key = dotted.split(".", 1)
    return DEFAULT_CONFIG[section][key]


def _is_baseline(overrides: dict[str, Any]) -> bool:
    return all(_default_value(key) == value for key, value in overrides.items())


def generate_ofat_matrix(
    axes: dict[str, list[dict[str, Any]]] | None = None,
) -> tuple[MatrixRun, ...]:
    selected_axes = axes or DEFAULT_AXES
    runs: list[MatrixRun] = []
    baseline_hash = sha256_json(DEFAULT_CONFIG)
    runs.append(MatrixRun("baseline", "baseline", "default", (), baseline_hash))
    seen: set[tuple[str, ...]] = {()}
    for axis, values in selected_axes.items():
        for value in values:
            if _is_baseline(value):
                continue
            overrides = tuple(
                f"{key}={str(item).lower() if isinstance(item, bool) else item}"
                for key, item in sorted(value.items())
            )
            if overrides in seen:
                continue
            seen.add(overrides)
            identity = sha256_json({"axis": axis, "overrides": overrides})[:12]
            config = {section: dict(values) for section, values in DEFAULT_CONFIG.items()}
            for key, item in value.items():
                section, name = key.split(".", 1)
                config[section][name] = item
            runs.append(
                MatrixRun(
                    f"{axis}-{identity}",
                    axis,
                    value,
                    overrides,
                    sha256_json(config),
                )
            )
    return tuple(runs)


def write_matrix(output: Path, runs: tuple[MatrixRun, ...]) -> tuple[Path, Path]:
    output.mkdir(parents=True, exist_ok=True)
    json_path = output / "parameter-matrix.json"
    tsv_path = output / "parameter-matrix.tsv"
    write_json(
        json_path,
        {
            "schema_version": "parameter-matrix-v1",
            "design": "one-factor-at-a-time",
            "run_count": len(runs),
            "runs": [
                {
                    "run_id": run.run_id,
                    "axis": run.axis,
                    "value": run.value,
                    "overrides": list(run.overrides),
                    "config_sha256": run.config_sha256,
                }
                for run in runs
            ],
        },
    )
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["run_id", "axis", "value", "overrides", "config_sha256"])
        for run in runs:
            writer.writerow(
                [run.run_id, run.axis, repr(run.value), ";".join(run.overrides), run.config_sha256]
            )
    return json_path, tsv_path
