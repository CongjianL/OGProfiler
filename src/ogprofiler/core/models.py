"""Serializable, implementation-independent core data models."""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields
from typing import Any, TypeVar

T = TypeVar("T", bound="SerializableModel")


class SerializableModel:
    """Mixin for deterministic JSON-compatible dataclass serialization."""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)  # type: ignore[call-overload,no-any-return]

    @classmethod
    def from_dict(cls: type[T], values: dict[str, Any]) -> T:
        allowed = {field.name for field in fields(cls)}  # type: ignore[arg-type]
        unknown = set(values) - allowed
        if unknown:
            names = ", ".join(sorted(unknown))
            raise ValueError(f"Unknown {cls.__name__} fields: {names}")
        return cls(**values)


@dataclass(frozen=True, slots=True)
class Protein(SerializableModel):
    protein_id: int
    species_id: int
    original_id: str
    length: int


@dataclass(frozen=True, slots=True)
class Species(SerializableModel):
    species_id: int
    species_name: str
    source_file: str


@dataclass(frozen=True, slots=True)
class DatasetManifest(SerializableModel):
    dataset_sha256: str
    input_checksums: dict[str, str]
    n_species: int
    n_proteins: int
    id_algorithm: str = "sorted-path+sorted-original-id-v1"


@dataclass(frozen=True, slots=True)
class RunManifest(SerializableModel):
    ogprofiler_version: str
    algorithm_version: str
    command: list[str]
    random_seed: int
    resolved_config_sha256: str
    dataset: DatasetManifest

    def to_dict(self) -> dict[str, Any]:
        values = asdict(self)
        return values

    @classmethod
    def from_dict(cls, values: dict[str, Any]) -> RunManifest:
        copied = dict(values)
        dataset = copied.get("dataset")
        if isinstance(dataset, dict):
            copied["dataset"] = DatasetManifest.from_dict(dataset)
        allowed = {field.name for field in fields(cls)}
        unknown = set(copied) - allowed
        if unknown:
            names = ", ".join(sorted(unknown))
            raise ValueError(f"Unknown {cls.__name__} fields: {names}")
        return cls(**copied)


@dataclass(frozen=True, slots=True)
class ComponentSummary(SerializableModel):
    component_id: int
    n_vertices: int
    n_edges: int
    n_species: int


@dataclass(frozen=True, slots=True)
class HierarchyNode(SerializableModel):
    cluster_id: int
    parent_id: int | None
    component_id: int
    depth: int
    n_genes: int
    n_species: int
    resolution: float | None = None
    quality: float | None = None
    child_count: int = 0
    split_status: str = "PENDING"
    network_event: str | None = None
    phylo_event: str | None = None
    terminal_reason: str | None = None
    search_status: str | None = None
    termination_kind: str | None = None
    failure_codes: tuple[str, ...] = ()
    selection_phase: str | None = None
    selection_kind: str | None = None
    refinement_truncated: bool = False


@dataclass(frozen=True, slots=True)
class SplitResult(SerializableModel):
    membership: tuple[int, ...]
    resolution: float
    quality: float
    child_count: int
    min_child_size: int
    largest_child_fraction: float
    accepted: bool
    reason: str
