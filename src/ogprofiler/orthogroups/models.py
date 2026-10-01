"""Typed V1 event view; no changes to network-event annotations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from ogprofiler.exceptions import InputError

V1Event = Literal["I", "II", "III-1", "III-2", "III-3"]


@dataclass(frozen=True, slots=True)
class V1EventAnnotation:
    component_id: int
    cluster_id: int
    reference_order: int
    undirected_degree: int
    eligible_child_count: int
    v1_event: V1Event | None

    @property
    def selection_event(self) -> str:
        """Explicit unrefined hnn_analysis normalization, not raw event state."""
        return "None" if self.v1_event is None else self.v1_event


@dataclass(frozen=True, slots=True)
class Orthogroup:
    component_id: int
    local_group_id: int
    source_cluster_id: int | None
    selection_type: Literal["EVENT_I", "RESIDUAL_NONE", "SSN_ISOLATE"]
    v1_event: V1Event | None
    processing_level: int
    protein_ids: tuple[int, ...]
    n_species: int
    membership_hash: str

    @property
    def n_genes(self) -> int:
        return len(self.protein_ids)


@dataclass(frozen=True, slots=True)
class SelectionTrace:
    cluster_id: int | None
    processing_level: int
    selection_event: str
    status: Literal["SELECTED", "SKIPPED_CONSUMED", "DESCENDANT_CONSUMED"]
    consumed_by: int | None


@dataclass(frozen=True, slots=True)
class UnassignedProtein:
    protein_id: int
    terminal_cluster_id: int
    reason: str


@dataclass(frozen=True, slots=True)
class OrthogroupResult:
    component_id: int
    groups: tuple[Orthogroup, ...]
    trace: tuple[SelectionTrace, ...]
    unassigned: tuple[UnassignedProtein, ...]
    remaining_cluster_ids: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class OrthogroupConfig:
    strategy: str = "v1_compatible"
    species_overlap_count: int = 0
    refinement: bool = False

    def __post_init__(self) -> None:
        if self.strategy != "v1_compatible":
            raise InputError("orthogroups.strategy must be v1_compatible")
        value = self.species_overlap_count
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise InputError("orthogroups.species_overlap_count must be a non-negative integer")
        if not isinstance(self.refinement, bool):
            raise InputError("orthogroups.refinement must be boolean")
        if self.refinement:
            raise InputError("orthogroups.refinement=true is not implemented")
