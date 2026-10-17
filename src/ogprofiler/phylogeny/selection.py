"""Deterministic family selection for optional phylogenetic refinement."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

from ogprofiler.core.manifest import sha256_file
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


def terminal_result_paths(run_root: Path) -> tuple[Path, Path, Path, Path]:
    """Refinement remains terminal-family evidence, independent of OG selection."""
    root = run_root / "results"
    if (root / "terminal_families.tsv").is_file():
        manifest_path = root / "export-manifest.json"
        try:
            manifest = json.loads(manifest_path.read_text())
            if manifest["parameters"]["strategy"] != "v1_compatible":
                raise PhylogenyError("Terminal diagnostic export is not complete")
            for name in (
                "terminal_families.tsv",
                "terminal_members.tsv",
                "hierarchy.tsv",
                "events.tsv",
            ):
                if sha256_file(root / name) != manifest["output_checksums"][name]:
                    raise PhylogenyError("Terminal diagnostics are corrupt; rerun export")
            for name, checksum in manifest["input_checksums"].items():
                if name.startswith(
                    ("hierarchy/", "evolution/components/", "components/", "input/proteins.parquet")
                ):
                    if sha256_file(run_root / name) != checksum:
                        raise PhylogenyError("Terminal diagnostics are stale; rerun export")
        except (OSError, ValueError, KeyError, TypeError) as error:
            raise PhylogenyError(f"Invalid terminal diagnostic export: {error}") from error
        return (
            root / "terminal_families.tsv",
            root / "terminal_members.tsv",
            root / "hierarchy.tsv",
            root / "events.tsv",
        )
    diagnostic = root / "terminal-families"
    if (diagnostic / "export-manifest.json").is_file():
        return (
            diagnostic / "families.tsv",
            diagnostic / "members.tsv",
            diagnostic / "hierarchy.tsv",
            diagnostic / "events.tsv",
        )
    # Explicitly retain pre-P4 terminal exchange fixtures/legacy workspaces only.
    families = root / "families.tsv"
    rows = _tsv(families)
    if rows and not {"cluster_id", "terminal_reason", "network_event"}.issubset(rows[0]):
        raise PhylogenyError("Refinement requires terminal-family diagnostics; rerun export")
    return families, root / "members.tsv", root / "hierarchy.tsv", root / "events.tsv"


def select_refinement_families(
    run_root: Path,
    *,
    explicit_family_ids: tuple[str, ...],
    selection_events: set[str],
    large_family_size: int,
    max_families: int,
) -> tuple[RefinementFamily, ...]:
    families_path, _, hierarchy_path, events_path = terminal_result_paths(run_root)
    family_rows = _tsv(families_path)
    hierarchy_rows = _tsv(hierarchy_path)
    event_rows = _tsv(events_path)
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
