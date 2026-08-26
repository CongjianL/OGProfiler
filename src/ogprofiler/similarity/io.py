"""Fixed Arrow schemas and batched readers/writers for Phase 4."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from ogprofiler.exceptions import EdgeConstructionError
from ogprofiler.search.hits import HIT_SCHEMA
from ogprofiler.similarity.models import DirectionalHit, NormalizedHit, RetainedEdge

NORMALIZED_HIT_SCHEMA = HIT_SCHEMA.append(pa.field("normalized_score", pa.float64()))

RETAINED_EDGE_SCHEMA = pa.schema(
    [
        ("u", pa.int64()),
        ("v", pa.int64()),
        ("u_species", pa.int32()),
        ("v_species", pa.int32()),
        ("score_uv", pa.float64()),
        ("score_vu", pa.float64()),
        ("weight", pa.float64()),
        ("coverage", pa.float32()),
        ("edge_type", pa.string()),
    ]
)


def iter_directional_hits(path: Path, batch_size: int = 65_536) -> Iterator[DirectionalHit]:
    parquet = pq.ParquetFile(path)
    if parquet.schema_arrow != HIT_SCHEMA:
        raise EdgeConstructionError(
            f"Directional hit schema mismatch in {path}: {parquet.schema_arrow}"
        )
    for batch in parquet.iter_batches(batch_size=batch_size):
        for row in batch.to_pylist():
            yield DirectionalHit(
                query_id=int(row["query_id"]),
                target_id=int(row["target_id"]),
                query_species=int(row["query_species"]),
                target_species=int(row["target_species"]),
                bitscore=float(row["bitscore"]),
                identity=float(row["identity"]),
                query_coverage=float(row["query_coverage"]),
                target_coverage=float(row["target_coverage"]),
                evalue=float(row["evalue"]),
            )


def write_normalized_hits(path: Path, hits: list[NormalizedHit]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pylist(
        [
            {
                "query_id": item.hit.query_id,
                "target_id": item.hit.target_id,
                "query_species": item.hit.query_species,
                "target_species": item.hit.target_species,
                "bitscore": item.hit.bitscore,
                "identity": item.hit.identity,
                "query_coverage": item.hit.query_coverage,
                "target_coverage": item.hit.target_coverage,
                "evalue": item.hit.evalue,
                "normalized_score": item.normalized_score,
            }
            for item in hits
        ],
        schema=NORMALIZED_HIT_SCHEMA,
    )
    pq.write_table(table, path, compression="zstd")


def write_retained_edges(path: Path, edges: list[RetainedEdge]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pylist(
        [
            {
                "u": edge.u,
                "v": edge.v,
                "u_species": edge.u_species,
                "v_species": edge.v_species,
                "score_uv": edge.score_uv,
                "score_vu": edge.score_vu,
                "weight": edge.weight,
                "coverage": edge.coverage,
                "edge_type": edge.edge_type,
            }
            for edge in edges
        ],
        schema=RETAINED_EDGE_SCHEMA,
    )
    pq.write_table(table, path, compression="zstd")
