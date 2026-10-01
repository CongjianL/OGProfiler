"""Recovery metrics and applicability maps for synthetic evolution fixtures."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from collections.abc import Iterator
from pathlib import Path
from statistics import mean
from typing import Any, cast

import pyarrow as pa

from ogprofiler.benchmark.metrics import family_metrics
from ogprofiler.core.manifest import sha256_file, write_json
from ogprofiler.exceptions import InputError
from ogprofiler.output.results import terminal_table_paths


def _rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise InputError(f"Missing synthetic benchmark input: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _ratio(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def _f1(precision: float | None, recall: float | None) -> float | None:
    if precision is None or recall is None:
        return None
    return 0.0 if precision + recall == 0 else 2 * precision * recall / (precision + recall)


def _decoded_lines(stream: pa.NativeFile, *, chunk_size: int = 1 << 20) -> Iterator[str]:
    pending = b""
    while chunk := stream.read(chunk_size):
        parts = (pending + chunk).splitlines(keepends=True)
        pending = b""
        if parts and not parts[-1].endswith((b"\n", b"\r")):
            pending = parts.pop()
        for line in parts:
            yield line.decode("utf-8")
    if pending:
        yield pending.decode("utf-8")


def _predicted_clades(
    run_root: Path,
) -> tuple[set[frozenset[str]], dict[frozenset[str], str]]:
    hierarchy = _rows(run_root / "results" / "hierarchy.tsv")
    terminal_families, terminal_members = terminal_table_paths(run_root)
    families = _rows(terminal_families)
    members = _rows(terminal_members)
    events = _rows(run_root / "results" / "events.tsv")
    family_members: dict[str, set[str]] = defaultdict(set)
    for row in members:
        family_members[row["family_id"]].add(row["original_id"])
    terminal = {
        (int(row["component_id"]), int(row["cluster_id"])): family_members[row["family_id"]]
        for row in families
    }
    children: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    nodes: set[tuple[int, int]] = set()
    for row in hierarchy:
        key = (int(row["component_id"]), int(row["cluster_id"]))
        nodes.add(key)
        if row["parent_id"]:
            children[(key[0], int(row["parent_id"]))].append(key)
    cache: dict[tuple[int, int], set[str]] = {}

    def descendants(key: tuple[int, int]) -> set[str]:
        if key not in cache:
            cache[key] = (
                set(terminal[key])
                if key in terminal
                else set().union(*(descendants(child) for child in children[key]))
            )
        return cache[key]

    clades = {frozenset(descendants(key)) for key in nodes if len(descendants(key)) >= 2}
    event_by_clade: dict[frozenset[str], str] = {}
    for row in events:
        key = (int(row["component_id"]), int(row["cluster_id"]))
        values = frozenset(descendants(key))
        if len(values) >= 2:
            event_by_clade[values] = row["network_event"]
    return clades, event_by_clade


def _true_clades(
    dataset_root: Path,
) -> tuple[set[frozenset[str]], list[tuple[frozenset[str], str]]]:
    genealogy = _rows(dataset_root / "genealogy.tsv")
    children: dict[str, list[str]] = defaultdict(list)
    proteins: dict[str, str] = {}
    events: dict[str, str] = {}
    for row in genealogy:
        if row["parent_id"]:
            children[row["parent_id"]].append(row["node_id"])
        if row["protein_id"]:
            proteins[row["node_id"]] = row["protein_id"]
        events[row["node_id"]] = row["event"]
    cache: dict[str, set[str]] = {}

    def descendants(node_id: str) -> set[str]:
        if node_id not in cache:
            cache[node_id] = (
                {proteins[node_id]}
                if node_id in proteins
                else set().union(*(descendants(child) for child in children[node_id]))
            )
        return cache[node_id]

    clades = {
        frozenset(descendants(node_id)) for node_id in events if len(descendants(node_id)) >= 2
    }
    event_clades = [
        (frozenset(descendants(node_id)), event)
        for node_id, event in events.items()
        if event in {"SPECIATION", "DUPLICATION"} and len(descendants(node_id)) >= 2
    ]
    return clades, event_clades


def hierarchy_and_event_metrics(run_root: Path, dataset_root: Path) -> dict[str, Any]:
    predicted, predicted_events = _predicted_clades(run_root)
    truth, truth_events = _true_clades(dataset_root)
    exact = len(predicted & truth)
    precision = _ratio(exact, len(predicted))
    recall = _ratio(exact, len(truth))
    best_jaccards = [
        max(
            (len(clade & candidate) / len(clade | candidate) for candidate in predicted),
            default=0.0,
        )
        for clade in truth
    ]
    matched_events = 0
    correct_events = 0
    confusion: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for clade, truth_event in truth_events:
        predicted_event = predicted_events.get(clade)
        if predicted_event is None:
            continue
        expected = "SPECIATION_LIKE" if truth_event == "SPECIATION" else "DUPLICATION_LIKE"
        matched_events += 1
        correct_events += int(predicted_event == expected)
        confusion[expected][predicted_event] += 1
    return {
        "hierarchy": {
            "truth_clades": len(truth),
            "predicted_clades": len(predicted),
            "exact_clades": exact,
            "precision": precision,
            "recall": recall,
            "f1": _f1(precision, recall),
            "mean_best_truth_jaccard": mean(best_jaccards) if best_jaccards else None,
        },
        "events": {
            "truth_events": len(truth_events),
            "matched_exact_clades": matched_events,
            "correct": correct_events,
            "coverage": _ratio(matched_events, len(truth_events)),
            "accuracy_on_matched": _ratio(correct_events, matched_events),
            "end_to_end_accuracy": _ratio(correct_events, len(truth_events)),
            "confusion": {key: dict(value) for key, value in sorted(confusion.items())},
        },
    }


def synthetic_orthology_metrics(run_root: Path, dataset_root: Path) -> dict[str, Any] | None:
    predicted_path = run_root / "results" / "ortholog_pairs.tsv.zst"
    if not predicted_path.is_file():
        return None
    truth = {
        tuple(sorted((row["protein_a"], row["protein_b"])))
        for row in _rows(dataset_root / "true_orthologs.tsv")
    }
    original_by_integer = {
        int(row["protein_id"]): row["original_id"]
        for row in _rows(terminal_table_paths(run_root)[1])
    }
    predicted: set[tuple[str, str]] = set()
    with pa.input_stream(str(predicted_path), compression="zstd") as stream:
        for row in csv.DictReader(_decoded_lines(stream), delimiter="\t"):
            left = original_by_integer[int(row["protein_a_id"])]
            right = original_by_integer[int(row["protein_b_id"])]
            predicted.add((left, right) if left < right else (right, left))
    intersection = len(predicted & truth)
    precision = intersection / len(predicted) if predicted else float(not truth)
    recall = intersection / len(truth) if truth else float(not predicted)
    return {
        "truth_pairs": len(truth),
        "predicted_pairs": len(predicted),
        "intersection": intersection,
        "precision": precision,
        "recall": recall,
        "f1": _f1(precision, recall),
    }


def evaluate_synthetic_run(
    run_root: Path, dataset_root: Path, *, method: str = "ogprofiler2"
) -> dict[str, Any]:
    try:
        manifest = json.loads((dataset_root / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise InputError(f"Failed to read synthetic manifest: {error}") from error
    truth = _rows(dataset_root / "ground_truth.tsv")
    families = _rows(run_root / "results" / "families.tsv")
    members = _rows(run_root / "results" / "members.tsv")
    structure = hierarchy_and_event_metrics(run_root, dataset_root)
    return {
        "schema_version": "synthetic-recovery-v1",
        "method": method,
        "scenario_id": manifest["scenario_id"],
        "scenario": manifest["scenario"],
        "input_checksums": {
            **{
                path.relative_to(run_root).as_posix(): sha256_file(path)
                for path in terminal_table_paths(run_root)
            },
            "dataset_manifest": sha256_file(dataset_root / "manifest.json"),
            "families.tsv": sha256_file(run_root / "results" / "families.tsv"),
            "members.tsv": sha256_file(run_root / "results" / "members.tsv"),
            "hierarchy.tsv": sha256_file(run_root / "results" / "hierarchy.tsv"),
            "events.tsv": sha256_file(run_root / "results" / "events.tsv"),
        },
        "family": family_metrics(families, members, truth, large_family_size=20),
        **structure,
        "orthology": synthetic_orthology_metrics(run_root, dataset_root),
    }


def write_synthetic_evaluation(path: Path, result: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_json(path, result)


def aggregate_applicability(
    metric_paths: list[Path],
    output: Path,
    *,
    family_f1_threshold: float = 0.90,
    hierarchy_f1_threshold: float = 0.80,
    event_accuracy_threshold: float = 0.80,
    orthology_f1_threshold: float = 0.90,
) -> tuple[Path, Path]:
    rows: list[dict[str, Any]] = []
    for path in metric_paths:
        try:
            item = cast(dict[str, Any], json.loads(path.read_text(encoding="utf-8")))
        except (OSError, json.JSONDecodeError) as error:
            raise InputError(f"Failed to read synthetic metrics {path}: {error}") from error
        if item.get("schema_version") != "synthetic-recovery-v1":
            raise InputError(f"Unsupported synthetic metrics schema in {path}")
        scenario = item["scenario"]
        family_f1 = item["family"]["pairwise_clustering"]["f1"]
        hierarchy_f1 = item["hierarchy"]["f1"]
        event_accuracy = item["events"]["end_to_end_accuracy"]
        orthology_f1 = (item.get("orthology") or {}).get("f1")
        family_applicable = family_f1 is not None and family_f1 >= family_f1_threshold
        hierarchy_applicable = hierarchy_f1 is not None and hierarchy_f1 >= hierarchy_f1_threshold
        event_applicable = event_accuracy is not None and event_accuracy >= event_accuracy_threshold
        orthology_applicable = orthology_f1 is not None and orthology_f1 >= orthology_f1_threshold
        rows.append(
            {
                "scenario_id": item["scenario_id"],
                **scenario,
                "family_f1": family_f1,
                "hierarchy_f1": hierarchy_f1,
                "event_accuracy": event_accuracy,
                "orthology_f1": orthology_f1,
                "family_applicable": family_applicable,
                "hierarchy_applicable": hierarchy_applicable,
                "event_applicable": event_applicable,
                "orthology_applicable": orthology_applicable,
                "overall_applicable": all(
                    (
                        family_applicable,
                        hierarchy_applicable,
                        event_applicable,
                        orthology_applicable,
                    )
                ),
            }
        )
    rows.sort(
        key=lambda row: (
            row["divergence"],
            row["duplication_rate"],
            row["loss_rate"],
            row["expansion"],
            row["fusion_rate"],
            row["replicate"],
        )
    )
    output.mkdir(parents=True, exist_ok=True)
    json_path = output / "leiden-applicability.json"
    tsv_path = output / "leiden-applicability.tsv"
    thresholds = {
        "family_f1": family_f1_threshold,
        "hierarchy_f1": hierarchy_f1_threshold,
        "event_accuracy": event_accuracy_threshold,
        "orthology_f1": orthology_f1_threshold,
    }
    write_json(
        json_path,
        {
            "schema_version": "leiden-applicability-v1",
            "thresholds": thresholds,
            "evaluated_scenarios": len(rows),
            "family_applicable_scenarios": sum(bool(row["family_applicable"]) for row in rows),
            "hierarchy_applicable_scenarios": sum(
                bool(row["hierarchy_applicable"]) for row in rows
            ),
            "event_applicable_scenarios": sum(bool(row["event_applicable"]) for row in rows),
            "orthology_applicable_scenarios": sum(
                bool(row["orthology_applicable"]) for row in rows
            ),
            "overall_applicable_scenarios": sum(bool(row["overall_applicable"]) for row in rows),
            "rows": rows,
        },
    )
    fields = list(rows[0]) if rows else ["scenario_id"]
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    axis_rows: list[dict[str, Any]] = []
    for axis in ("divergence", "duplication_rate", "loss_rate", "expansion", "fusion_rate"):
        grouped: dict[float | int, list[dict[str, Any]]] = defaultdict(list)
        for row in rows:
            grouped[row[axis]].append(row)
        for value, group in sorted(grouped.items()):
            axis_row: dict[str, Any] = {
                "axis": axis,
                "value": value,
                "scenarios": len(group),
            }
            for metric in ("family_f1", "hierarchy_f1", "event_accuracy", "orthology_f1"):
                values = [float(row[metric]) for row in group if row[metric] is not None]
                axis_row[f"mean_{metric}"] = mean(values) if values else None
            for decision in (
                "family_applicable",
                "hierarchy_applicable",
                "event_applicable",
                "orthology_applicable",
                "overall_applicable",
            ):
                axis_row[f"{decision}_count"] = sum(bool(row[decision]) for row in group)
            axis_rows.append(axis_row)
    write_json(
        output / "axis-summary.json",
        {"schema_version": "synthetic-axis-summary-v1", "rows": axis_rows},
    )
    axis_fields = list(axis_rows[0]) if axis_rows else ["axis", "value"]
    with (output / "axis-summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, axis_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(axis_rows)
    return json_path, tsv_path
