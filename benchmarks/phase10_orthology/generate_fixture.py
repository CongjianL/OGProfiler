#!/usr/bin/env python3
"""Generate a deterministic four-child orthology streaming fixture."""

from __future__ import annotations

import argparse
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from ogprofiler.evolution.stage import EVENT_SCHEMA


def build_fixture(root: Path, genes_per_species: int = 100) -> None:
    input_root = root / "input"
    hierarchy_root = root / "hierarchy/components/component=00000000"
    evolution_root = root / "evolution/components/component=00000000"
    input_root.mkdir(parents=True, exist_ok=True)
    hierarchy_root.mkdir(parents=True, exist_ok=True)
    evolution_root.mkdir(parents=True, exist_ok=True)

    proteins = []
    memberships = []
    for species_id in range(4):
        for offset in range(genes_per_species):
            protein_id = species_id * genes_per_species + offset
            proteins.append(
                {
                    "protein_id": protein_id,
                    "species_id": species_id,
                    "original_id": f"species_{species_id}_gene_{offset:04d}",
                    "length": 100,
                }
            )
            memberships.append(
                {"protein_id": protein_id, "terminal_cluster_id": species_id + 1}
            )
    pq.write_table(pa.Table.from_pylist(proteins), input_root / "proteins.parquet")

    nodes = [
        {
            "cluster_id": 0,
            "parent_id": None,
            "component_id": 0,
            "depth": 0,
            "n_genes": 4 * genes_per_species,
            "n_species": 4,
            "resolution": 1.0,
            "quality": 1.0,
            "child_count": 4,
            "split_status": "SPLIT",
            "terminal_reason": None,
        }
    ]
    for species_id in range(4):
        nodes.append(
            {
                "cluster_id": species_id + 1,
                "parent_id": 0,
                "component_id": 0,
                "depth": 1,
                "n_genes": genes_per_species,
                "n_species": 1,
                "resolution": None,
                "quality": None,
                "child_count": 0,
                "split_status": "TERMINAL",
                "terminal_reason": "ONE_SPECIES",
            }
        )
    pq.write_table(pa.Table.from_pylist(nodes), hierarchy_root / "nodes.parquet")
    pq.write_table(pa.Table.from_pylist(memberships), hierarchy_root / "members.parquet")

    events = [
        {
            "component_id": 0,
            "cluster_id": 0,
            "child_count": 4,
            "network_event": "POLYTOMY",
            "overlap_count": 0,
            "overlap_score": 0.0,
            "pairwise_overlap_summary": "[]",
            "confidence": 1.0,
            "legacy_event": "III-3",
        }
    ]
    for species_id in range(4):
        events.append(
            {
                "component_id": 0,
                "cluster_id": species_id + 1,
                "child_count": 0,
                "network_event": "SPECIES_SPECIFIC",
                "overlap_count": 0,
                "overlap_score": 0.0,
                "pairwise_overlap_summary": "[]",
                "confidence": 1.0,
                "legacy_event": "",
            }
        )
    pq.write_table(
        pa.Table.from_pylist(events, schema=EVENT_SCHEMA),
        evolution_root / "events.parquet",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--genes-per-species", type=int, default=100)
    args = parser.parse_args()
    build_fixture(args.out, args.genes_per_species)


if __name__ == "__main__":
    main()
