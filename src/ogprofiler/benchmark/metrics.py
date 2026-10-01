"""Family, evolutionary-consistency, and reference orthology metrics."""

from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from collections.abc import Iterator
from pathlib import Path
from statistics import mean, median
from typing import Any

import pyarrow as pa

from ogprofiler.core.manifest import sha256_file, write_json
from ogprofiler.exceptions import InputError
from ogprofiler.output.results import terminal_table_paths


def _tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise InputError(f"Missing benchmark input: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _choose2(value: int) -> int:
    return value * (value - 1) // 2


def _canonical_pair(left: str, right: str) -> tuple[str, str]:
    return (left, right) if left < right else (right, left)


def _decoded_lines(stream: pa.NativeFile, *, chunk_size: int = 1 << 20) -> Iterator[str]:
    """Yield UTF-8 lines from an Arrow stream without materializing the file."""
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


def _ratio(numerator: int | float, denominator: int | float) -> float | None:
    return float(numerator) / float(denominator) if denominator else None


def _f1(precision: float | None, recall: float | None) -> float | None:
    if precision is None or recall is None or precision + recall == 0:
        return None if precision is None or recall is None else 0.0
    return 2 * precision * recall / (precision + recall)


def _quantiles(values: list[int]) -> dict[str, float | int | None]:
    if not values:
        return {"min": None, "q25": None, "median": None, "q75": None, "max": None, "mean": None}
    ordered = sorted(values)

    def percentile(fraction: float) -> float:
        position = fraction * (len(ordered) - 1)
        lower, upper = math.floor(position), math.ceil(position)
        return ordered[lower] + (position - lower) * (ordered[upper] - ordered[lower])

    return {
        "min": ordered[0],
        "q25": percentile(0.25),
        "median": median(ordered),
        "q75": percentile(0.75),
        "max": ordered[-1],
        "mean": mean(ordered),
    }


def _group_members(rows: list[dict[str, str]], key: str, member: str) -> dict[str, set[str]]:
    result: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        result[row[key]].add(row[member])
    return dict(result)


def family_metrics(
    families: list[dict[str, str]],
    members: list[dict[str, str]],
    truth: list[dict[str, str]],
    *,
    large_family_size: int,
) -> dict[str, Any]:
    predicted = _group_members(members, "family_id", "original_id")
    true = _group_members(truth, "family", "protein_id")
    prediction_by_member = {
        member: family_id for family_id, values in predicted.items() for member in values
    }
    truth_by_member = {member: family_id for family_id, values in true.items() for member in values}
    shared = set(prediction_by_member) & set(truth_by_member)
    contingency = Counter(
        (prediction_by_member[member], truth_by_member[member]) for member in shared
    )
    predicted_pairs = sum(_choose2(len(values & shared)) for values in predicted.values())
    true_pairs = sum(_choose2(len(values & shared)) for values in true.values())
    intersection = sum(_choose2(value) for value in contingency.values())
    precision = _ratio(intersection, predicted_pairs)
    recall = _ratio(intersection, true_pairs)
    predicted_sets = set(map(frozenset, predicted.values()))
    exact_truth = sum(frozenset(values) in predicted_sets for values in true.values())
    species_count = len({row["species"] for row in truth})
    truth_species: dict[str, Counter[str]] = defaultdict(Counter)
    for row in truth:
        truth_species[row["family"]][row["species"]] += 1
    single_copy = {
        family_id
        for family_id, counts in truth_species.items()
        if len(counts) == species_count and set(counts.values()) == {1}
    }
    exact_single = sum(frozenset(true[family_id]) in predicted_sets for family_id in single_copy)
    large = {key: value for key, value in true.items() if len(value) >= large_family_size}
    fragments = [
        len({prediction_by_member[member] for member in values if member in prediction_by_member})
        for values in large.values()
    ]
    coverages = [int(row["n_species"]) / species_count for row in families]
    sizes = [int(row["n_genes"]) for row in families]
    return {
        "family_number": len(families),
        "family_size_distribution": _quantiles(sizes),
        "species_coverage": {
            "mean": mean(coverages) if coverages else None,
            "complete_fraction": _ratio(sum(value == 1.0 for value in coverages), len(coverages)),
        },
        "membership_coverage": _ratio(len(shared), len(truth_by_member)),
        "pairwise_clustering": {
            "predicted_pairs": predicted_pairs,
            "true_pairs": true_pairs,
            "intersection": intersection,
            "precision": precision,
            "recall": recall,
            "f1": _f1(precision, recall),
        },
        "exact_family_recovery": {
            "count": exact_truth,
            "total": len(true),
            "rate": _ratio(exact_truth, len(true)),
        },
        "single_copy_family_recovery": {
            "count": exact_single,
            "total": len(single_copy),
            "rate": _ratio(exact_single, len(single_copy)),
        },
        "large_family_fragmentation": {
            "threshold": large_family_size,
            "families": len(large),
            "mean_fragments": mean(fragments) if fragments else None,
            "max_fragments": max(fragments) if fragments else None,
            "unfragmented_fraction": _ratio(sum(value == 1 for value in fragments), len(fragments)),
        },
    }


def evolution_metrics(run_root: Path) -> dict[str, Any]:
    hierarchy = _tsv(run_root / "results" / "hierarchy.tsv")
    events = _tsv(run_root / "results" / "events.tsv")
    terminal_families, terminal_members = terminal_table_paths(run_root)
    families = _tsv(terminal_families)
    members = _tsv(terminal_members)
    species_by_family: dict[str, set[int]] = defaultdict(set)
    for row in members:
        species_by_family[row["family_id"]].add(int(row["species_id"]))
    terminal_species = {
        (int(row["component_id"]), int(row["cluster_id"])): species_by_family[row["family_id"]]
        for row in families
    }
    nodes = {(int(row["component_id"]), int(row["cluster_id"])): row for row in hierarchy}
    children: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    for key, row in nodes.items():
        if row["parent_id"]:
            children[(key[0], int(row["parent_id"]))].append(key)
    cache: dict[tuple[int, int], set[int]] = {}

    def species(key: tuple[int, int]) -> set[int]:
        if key in cache:
            return cache[key]
        if key in terminal_species:
            cache[key] = set(terminal_species[key])
        else:
            cache[key] = set().union(*(species(child) for child in children[key]))
        return cache[key]

    consistent = 0
    evaluated = 0
    for row in events:
        key = (int(row["component_id"]), int(row["cluster_id"]))
        child_keys = children.get(key, [])
        if len(child_keys) < 2:
            continue
        child_species = [species(child) for child in child_keys]
        overlaps = [
            bool(left & right)
            for left_index, left in enumerate(child_species)
            for right in child_species[left_index + 1 :]
        ]
        event = row["network_event"]
        expected = (
            (event == "SPECIATION_LIKE" and len(child_keys) == 2 and not any(overlaps))
            or (event == "POLYTOMY" and len(child_keys) > 2 and not any(overlaps))
            or (event in {"MIXED", "DUPLICATION_LIKE"} and any(overlaps))
        )
        consistent += int(expected)
        evaluated += 1
    result: dict[str, Any] = {
        "species_overlap_consistency": _ratio(consistent, evaluated),
        "evaluated_network_nodes": evaluated,
    }
    phylo_path = run_root / "evolution" / "phylogenetic" / "phylogenetic-events.tsv"
    if phylo_path.is_file():
        phylo = _tsv(phylo_path)
        comparable = [row for row in phylo if row["conflict_status"] != "UNRESOLVED"]
        resolved = [row for row in phylo if row["phylo_event"] != "UNRESOLVED"]
        result.update(
            {
                "gene_tree_concordance": _ratio(
                    sum(row["conflict_status"] == "CONCORDANT" for row in comparable),
                    len(comparable),
                ),
                "reconciliation_consistency": _ratio(len(resolved), len(phylo)),
                "mean_reconciliation_confidence": (
                    mean(float(row["event_confidence"]) for row in resolved) if resolved else None
                ),
                "refined_families": len(phylo),
            }
        )
    else:
        result.update(
            {
                "gene_tree_concordance": None,
                "reconciliation_consistency": None,
                "mean_reconciliation_confidence": None,
                "refined_families": 0,
            }
        )
    return result


def orthology_metrics(run_root: Path, truth: list[dict[str, str]]) -> dict[str, Any] | None:
    path = run_root / "results" / "ortholog_pairs.tsv.zst"
    if not path.is_file():
        return None
    by_family: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in truth:
        by_family[row["family"]].append(row)
    species_count = len({row["species"] for row in truth})
    scorable = {
        family_id: rows
        for family_id, rows in by_family.items()
        if len(rows) == species_count and len({row["species"] for row in rows}) == species_count
    }
    true_pairs: set[tuple[str, str]] = {
        _canonical_pair(left["protein_id"], right["protein_id"])
        for rows in scorable.values()
        for left_index, left in enumerate(rows)
        for right in rows[left_index + 1 :]
        if left["species"] != right["species"]
    }
    evaluated_proteins = {row["protein_id"] for rows in scorable.values() for row in rows}
    original_by_integer = {
        int(row["protein_id"]): row["original_id"]
        for row in _tsv(terminal_table_paths(run_root)[1])
    }
    predicted: set[tuple[str, str]] = set()
    ignored = 0
    with pa.input_stream(str(path), compression="zstd") as stream:
        for row in csv.DictReader(_decoded_lines(stream), delimiter="\t"):
            left = original_by_integer[int(row["protein_a_id"])]
            right = original_by_integer[int(row["protein_b_id"])]
            if left not in evaluated_proteins or right not in evaluated_proteins:
                ignored += 1
                continue
            predicted.add(_canonical_pair(left, right))
    intersection = len(predicted & true_pairs)
    precision = _ratio(intersection, len(predicted))
    recall = _ratio(intersection, len(true_pairs))
    return {
        "truth_policy": "complete-single-copy-ground-truth-families-v1",
        "scorable_families": len(scorable),
        "true_pairs": len(true_pairs),
        "predicted_pairs": len(predicted),
        "ignored_unscorable_predictions": ignored,
        "intersection": intersection,
        "precision": precision,
        "recall": recall,
        "f1": _f1(precision, recall),
    }


def evaluate_run(
    run_root: Path,
    ground_truth_path: Path,
    *,
    method: str,
    dataset: str,
    large_family_size: int = 20,
) -> dict[str, Any]:
    truth = _tsv(ground_truth_path)
    families = _tsv(run_root / "results" / "families.tsv")
    members = _tsv(run_root / "results" / "members.tsv")
    inputs = {
        "families.tsv": sha256_file(run_root / "results" / "families.tsv"),
        "members.tsv": sha256_file(run_root / "results" / "members.tsv"),
        "hierarchy.tsv": sha256_file(run_root / "results" / "hierarchy.tsv"),
        "events.tsv": sha256_file(run_root / "results" / "events.tsv"),
        "ground_truth.tsv": sha256_file(ground_truth_path),
    }
    inputs.update(
        {
            path.relative_to(run_root).as_posix(): sha256_file(path)
            for path in terminal_table_paths(run_root)
        }
    )
    result = {
        "schema_version": "scientific-benchmark-v1",
        "method": method,
        "dataset": dataset,
        "input_checksums": inputs,
        "family": family_metrics(families, members, truth, large_family_size=large_family_size),
        "evolution": evolution_metrics(run_root),
        "orthology": orthology_metrics(run_root, truth),
    }
    return result


def write_evaluation(path: Path, result: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_json(path, result)


def evaluate_legacy_membership(
    normalized_hierarchy_path: Path,
    ground_truth_path: Path,
    *,
    method: str,
    dataset: str,
    large_family_size: int = 20,
) -> dict[str, Any]:
    try:
        normalized = json.loads(normalized_hierarchy_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise InputError(
            f"Failed to read normalized legacy hierarchy {normalized_hierarchy_path}: {error}"
        ) from error
    membership = {str(key): str(value) for key, value in normalized["membership"].items()}
    truth = _tsv(ground_truth_path)
    species_by_protein = {row["protein_id"]: row["species"] for row in truth}
    species_names = sorted({row["species"] for row in truth})
    species_id = {name: index for index, name in enumerate(species_names)}
    by_family: dict[str, list[str]] = defaultdict(list)
    for protein_id, family_id in membership.items():
        by_family[family_id].append(protein_id)
    families: list[dict[str, str]] = []
    members: list[dict[str, str]] = []
    for rank, (_family_id, proteins) in enumerate(sorted(by_family.items())):
        exported_id = f"V1_{rank:09d}"
        families.append(
            {
                "family_id": exported_id,
                "n_genes": str(len(proteins)),
                "n_species": str(len({species_by_protein[value] for value in proteins})),
            }
        )
        for protein in sorted(proteins):
            species = species_by_protein[protein]
            members.append(
                {
                    "family_id": exported_id,
                    "protein_id": protein,
                    "species_id": str(species_id[species]),
                    "original_id": protein,
                }
            )
    return {
        "schema_version": "scientific-benchmark-v1",
        "method": method,
        "dataset": dataset,
        "input_checksums": {
            "normalized_hierarchy.json": sha256_file(normalized_hierarchy_path),
            "ground_truth.tsv": sha256_file(ground_truth_path),
        },
        "family": family_metrics(families, members, truth, large_family_size=large_family_size),
        "evolution": {
            "species_overlap_consistency": None,
            "evaluated_network_nodes": 0,
            "gene_tree_concordance": None,
            "reconciliation_consistency": None,
            "mean_reconciliation_confidence": None,
            "refined_families": 0,
        },
        "orthology": None,
    }
