"""Convert a reference SSN into the Phase 4 retained-edge schema for Phase 5 tests."""

from __future__ import annotations

import argparse
from pathlib import Path

import pyarrow.parquet as pq

from ogprofiler.graph.legacy import import_legacy_ssn
from ogprofiler.similarity.io import write_retained_edges
from ogprofiler.similarity.models import RetainedEdge


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--ssn", required=True, type=Path)
    args = parser.parse_args()

    protein_rows = pq.read_table(
        args.run / "input" / "proteins.parquet", columns=["protein_id", "species_id"]
    ).to_pylist()
    species = {int(row["protein_id"]): int(row["species_id"]) for row in protein_rows}
    edge_table = import_legacy_ssn(args.ssn)
    if set(edge_table.vertices) != set(species):
        raise SystemExit("Reference SSN vertex IDs do not match the prepared protein index")
    edges = [
        RetainedEdge(
            u=edge.source,
            v=edge.target,
            u_species=species[edge.source],
            v_species=species[edge.target],
            score_uv=edge.weight,
            score_vu=edge.weight,
            weight=edge.weight,
            coverage=100.0,
            edge_type="REFERENCE",
        )
        for edge in edge_table.edges
    ]
    write_retained_edges(args.run / "edges" / "retained_edges.parquet", edges)
    print(f"Prepared {len(edges)} reference edges for {len(edge_table.vertices)} proteins")


if __name__ == "__main__":
    main()
