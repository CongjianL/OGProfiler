"""Cross-method benchmark schema and compact comparison aggregation."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, cast

from ogprofiler.core.manifest import write_json
from ogprofiler.exceptions import InputError

METHOD_REGISTRY = {
    "ogprofiler1": {"status": "implemented", "input": "normalized V1 artifacts"},
    "ogprofiler2": {"status": "implemented", "input": "scientific-benchmark-v1"},
    "orthofinder": {"status": "planned-import", "input": "Orthogroups.tsv/Orthologues"},
    "sonicparanoid": {"status": "planned-import", "input": "ortholog_groups.tsv"},
}


def _load(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise InputError(f"Failed to read benchmark metrics {path}: {error}") from error
    if value.get("schema_version") != "scientific-benchmark-v1":
        raise InputError(f"Unsupported benchmark schema in {path}")
    return cast(dict[str, Any], value)


def compare_methods(paths: list[Path], output: Path) -> tuple[Path, Path]:
    metrics = [_load(path) for path in paths]
    rows: list[dict[str, Any]] = []
    for item in metrics:
        family = item["family"]
        orthology = item.get("orthology") or {}
        values = [
            family["pairwise_clustering"].get("f1"),
            family["exact_family_recovery"].get("rate"),
            family["single_copy_family_recovery"].get("rate"),
            item["evolution"].get("species_overlap_consistency"),
            orthology.get("f1"),
        ]
        available = [float(value) for value in values if value is not None]
        rows.append(
            {
                "method": item["method"],
                "dataset": item["dataset"],
                "family_pair_f1": family["pairwise_clustering"].get("f1"),
                "exact_family_recovery": family["exact_family_recovery"].get("rate"),
                "single_copy_recovery": family["single_copy_family_recovery"].get("rate"),
                "species_overlap_consistency": item["evolution"].get(
                    "species_overlap_consistency"
                ),
                "orthology_f1": orthology.get("f1"),
                "composite_score": sum(available) / len(available) if available else None,
            }
        )
    rows.sort(
        key=lambda row: (
            row["dataset"],
            -(row["composite_score"] if row["composite_score"] is not None else -1.0),
            row["method"],
        )
    )
    output.mkdir(parents=True, exist_ok=True)
    json_path = output / "method-comparison.json"
    tsv_path = output / "method-comparison.tsv"
    write_json(
        json_path,
        {
            "schema_version": "method-comparison-v1",
            "method_registry": METHOD_REGISTRY,
            "rows": rows,
        },
    )
    fields = list(rows[0]) if rows else ["method", "dataset"]
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return json_path, tsv_path
