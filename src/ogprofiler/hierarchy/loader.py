"""Memory-conscious component-local graph loading for hierarchy workers."""

from __future__ import annotations

import os
import uuid
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path

import igraph as ig
import numpy as np
import numpy.typing as npt
import pyarrow.parquet as pq

from ogprofiler.exceptions import ComponentError
from ogprofiler.graph.components import Component
from ogprofiler.graph.partition import load_component_edge_table

Int64Array = npt.NDArray[np.int64]
Int32Array = npt.NDArray[np.int32]


@dataclass(frozen=True, slots=True)
class ComponentLoadConfig:
    """Choose bounded-overhead component representations."""

    memory_map_arrow: bool = True
    numpy_memmap: bool = True


class LocalArrayLookup(Mapping[int, int]):
    """Read-only global-protein lookup backed by sorted dense local arrays."""

    __slots__ = ("_global_ids", "_values")

    def __init__(self, global_ids: Int64Array, values: Int64Array | Int32Array) -> None:
        if global_ids.ndim != 1 or values.ndim != 1 or len(global_ids) != len(values):
            raise ValueError("Local lookup arrays must be aligned one-dimensional vectors")
        if len(global_ids) > 1 and bool(np.any(global_ids[1:] <= global_ids[:-1])):
            raise ValueError("Global protein IDs must be unique and sorted")
        self._global_ids = global_ids
        self._values = values

    def __getitem__(self, protein_id: int) -> int:
        local_id = int(np.searchsorted(self._global_ids, protein_id))
        if local_id >= len(self._global_ids) or int(self._global_ids[local_id]) != protein_id:
            raise KeyError(protein_id)
        return int(self._values[local_id])

    def __iter__(self) -> Iterator[int]:
        return (int(value) for value in self._global_ids)

    def __len__(self) -> int:
        return len(self._global_ids)


class SpeciesBitmapLookup(Mapping[int, int]):
    """Compatibility view that computes bitmaps without storing Python integers."""

    __slots__ = ("_species",)

    def __init__(self, species: LocalArrayLookup) -> None:
        self._species = species

    def __getitem__(self, protein_id: int) -> int:
        return 1 << self._species[protein_id]

    def __iter__(self) -> Iterator[int]:
        return iter(self._species)

    def __len__(self) -> int:
        return len(self._species)


@dataclass(frozen=True, slots=True)
class LoadedComponentGraph:
    component: Component
    graph: ig.Graph
    local_to_global: tuple[int, ...]
    global_ids: Int64Array
    species_ids: Int32Array
    global_to_local: Mapping[int, int]
    species_by_protein: Mapping[int, int]
    species_bitmap_by_protein: Mapping[int, int]


def _atomic_save(path: Path, values: npt.NDArray[np.generic]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        with temporary.open("wb") as stream:
            np.save(stream, values, allow_pickle=False)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


class ComponentGraphLoader:
    def __init__(
        self, run_root: Path, config: ComponentLoadConfig | None = None
    ) -> None:
        self.run_root = run_root
        self.config = config or ComponentLoadConfig()

    def _local_arrays(
        self, component_id: int, local_to_global: tuple[int, ...], species: dict[int, int]
    ) -> tuple[Int64Array, Int32Array]:
        global_values = np.asarray(local_to_global, dtype=np.int64)
        species_values = np.asarray(
            [species[protein_id] for protein_id in local_to_global], dtype=np.int32
        )
        if not self.config.numpy_memmap:
            return global_values, species_values

        cache = self.run_root / "hierarchy" / "local-index"
        global_path = cache / f"component-{component_id:08d}-global.npy"
        species_path = cache / f"component-{component_id:08d}-species.npy"
        if not global_path.is_file() or not species_path.is_file():
            _atomic_save(global_path, global_values)
            _atomic_save(species_path, species_values)
        cached_global = np.load(global_path, mmap_mode="r", allow_pickle=False)
        cached_species = np.load(species_path, mmap_mode="r", allow_pickle=False)
        if not np.array_equal(cached_global, global_values) or not np.array_equal(
            cached_species, species_values
        ):
            _atomic_save(global_path, global_values)
            _atomic_save(species_path, species_values)
            cached_global = np.load(global_path, mmap_mode="r", allow_pickle=False)
            cached_species = np.load(species_path, mmap_mode="r", allow_pickle=False)
        return cached_global, cached_species

    def load(self, component_id: int) -> LoadedComponentGraph:
        edge_table = load_component_edge_table(
            self.run_root / "components",
            component_id,
            memory_map=self.config.memory_map_arrow,
        )
        graph, local_to_global = edge_table.to_igraph()
        protein_rows = pq.read_table(
            self.run_root / "input" / "proteins.parquet",
            columns=["protein_id", "species_id"],
            filters=[("protein_id", "in", list(local_to_global))],
            memory_map=self.config.memory_map_arrow,
        ).to_pylist()
        species = {int(row["protein_id"]): int(row["species_id"]) for row in protein_rows}
        missing = set(local_to_global).difference(species)
        if missing:
            raise ComponentError(
                f"Missing protein metadata for component {component_id}: {min(missing)}"
            )
        global_ids, species_ids = self._local_arrays(component_id, local_to_global, species)
        local_ids = np.arange(len(global_ids), dtype=np.int64)
        species_lookup = LocalArrayLookup(global_ids, species_ids)
        component = Component(
            component_id=component_id,
            vertices=edge_table.vertices,
            edges=edge_table.edges,
            original_ids=edge_table.original_ids,
        )
        return LoadedComponentGraph(
            component=component,
            graph=graph,
            local_to_global=local_to_global,
            global_ids=global_ids,
            species_ids=species_ids,
            global_to_local=LocalArrayLookup(global_ids, local_ids),
            species_by_protein=species_lookup,
            species_bitmap_by_protein=SpeciesBitmapLookup(species_lookup),
        )
