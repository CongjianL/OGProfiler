"""Deterministic family selection for optional phylogenetic refinement."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from ogprofiler.exceptions import PhylogenyError


@dataclass(frozen=True, slots=True)
class RefinementFamily:
    family_id: str
    component_id: int
    cluster_id: int
    n_genes: int
    n_species: int
    network_event: str
    selection_reasons: tuple[str, ...]


def _tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise PhylogenyError(f"Missing refinement input: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def select_refinement_families(
    run_root: Path,
    *,
    explicit_family_ids: tuple[str, ...],
    selection_events: set[str],
    large_family_size: int,
    max_families: int,
) -> tuple[RefinementFamily, ...]:
    family_rows = _tsv(run_root / "results" / "families.tsv")
    hierarchy_rows = _tsv(run_root / "results" / "hierarchy.tsv")
    event_rows = _tsv(run_root / "results" / "events.tsv")
    known = {row["family_id"] for row in family_rows}
    unknown = sorted(set(explicit_family_ids) - known)
    if unknown:
        raise PhylogenyError(f"Unknown family selected for refinement: {unknown[0]}")

    parents = {
        (int(row["component_id"]), int(row["cluster_id"])): (
            None if not row["parent_id"] else int(row["parent_id"])
        )
        for row in hierarchy_rows
    }
    events = {
        (int(row["component_id"]), int(row["cluster_id"])): row["network_event"]
        for row in event_rows
    }

    def duplication_rich(component_id: int, cluster_id: int) -> bool:
        current: int | None = cluster_id
        visited: set[int] = set()
        while current is not None:
            if current in visited:
                raise PhylogenyError(f"Cycle in exported hierarchy component {component_id}")
            visited.add(current)
            if events.get((component_id, current)) == "DUPLICATION_LIKE":
                return True
            current = parents.get((component_id, current))
        return False

    explicit = set(explicit_family_ids)
    candidates: list[tuple[tuple[int, int, int, int, bytes], RefinementFamily, bool]] = []
    for row in family_rows:
        family_id = row["family_id"]
        component_id = int(row["component_id"])
        cluster_id = int(row["cluster_id"])
        n_genes = int(row["n_genes"])
        network_event = row["network_event"]
        reasons: list[str] = []
        is_explicit = family_id in explicit
        if explicit and not is_explicit:
            continue
        if is_explicit:
            reasons.append("USER_SELECTED")
        if network_event in selection_events:
            reasons.append(network_event)
        is_large = n_genes >= large_family_size
        if is_large:
            reasons.append("LARGE_FAMILY")
        is_duplication_rich = duplication_rich(component_id, cluster_id)
        if is_duplication_rich:
            reasons.append("DUPLICATION_RICH")
        if not reasons:
            continue
        family = RefinementFamily(
            family_id,
            component_id,
            cluster_id,
            n_genes,
            int(row["n_species"]),
            network_event,
            tuple(reasons),
        )
        priority = (
            0 if is_explicit else 1,
            {"AMBIGUOUS": 0, "MIXED": 1}.get(network_event, 2),
            0 if is_large else 1,
            0 if is_duplication_rich else 1,
            family_id.encode("utf-8"),
        )
        candidates.append((priority, family, is_explicit))

    candidates.sort(key=lambda value: value[0])
    selected: list[RefinementFamily] = []
    automatic = 0
    for _, family, is_explicit in candidates:
        if not is_explicit:
            if automatic >= max_families:
                continue
            automatic += 1
        selected.append(family)
    return tuple(selected)
