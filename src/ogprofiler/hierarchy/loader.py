"""Component-local graph loader for production hierarchy workers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import igraph as ig
import pyarrow.parquet as pq

from ogprofiler.graph.components import Component
from ogprofiler.graph.partition import load_component_edge_table


@dataclass(frozen=True, slots=True)
class LoadedComponentGraph:
    component: Component
    graph: ig.Graph
    local_to_global: tuple[int, ...]
    global_to_local: dict[int, int]
    species_by_protein: dict[int, int]
    species_bitmap_by_protein: dict[int, int]


class ComponentGraphLoader:
    def __init__(self, run_root: Path) -> None:
        self.run_root = run_root

    def load(self, component_id: int) -> LoadedComponentGraph:
        edge_table = load_component_edge_table(self.run_root / "components", component_id)
        graph, local_to_global = edge_table.to_igraph()
        wanted = set(local_to_global)
        protein_rows = pq.read_table(
            self.run_root / "input" / "proteins.parquet",
            columns=["protein_id", "species_id"],
            filters=[("protein_id", "in", sorted(wanted))],
        ).to_pylist()
        species = {int(row["protein_id"]): int(row["species_id"]) for row in protein_rows}
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
            global_to_local={protein_id: index for index, protein_id in enumerate(local_to_global)},
            species_by_protein=species,
            species_bitmap_by_protein={
                protein_id: 1 << value for protein_id, value in species.items()
            },
        )
