"""Stream backend hit output into the standard directional Parquet schema."""

from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from ogprofiler.exceptions import SearchError

HIT_COLUMNS = (
    "query_id",
    "target_id",
    "query_species",
    "target_species",
    "bitscore",
    "identity",
    "query_coverage",
    "target_coverage",
    "evalue",
)

HIT_SCHEMA = pa.schema(
    [
        ("query_id", pa.int64()),
        ("target_id", pa.int64()),
        ("query_species", pa.int32()),
        ("target_species", pa.int32()),
        ("bitscore", pa.float64()),
        ("identity", pa.float32()),
        ("query_coverage", pa.float32()),
        ("target_coverage", pa.float32()),
        ("evalue", pa.float64()),
    ]
)

HitRow = tuple[int, int, int, int, float, float, float, float, float]


def _protein_id(value: str) -> int:
    prefix = "OGP2P"
    if not value.startswith(prefix) or not value[len(prefix) :].isdigit():
        raise SearchError(f"Unexpected prepared protein identifier: {value!r}")
    return int(value[len(prefix) :])


def _hit_table(rows: list[HitRow]) -> pa.Table:
    columns = list(zip(*rows, strict=True)) if rows else [()] * len(HIT_COLUMNS)
    return pa.table(
        {
            name: pa.array(values, type=field.type)
            for name, values, field in zip(HIT_COLUMNS, columns, HIT_SCHEMA, strict=True)
        },
        schema=HIT_SCHEMA,
    )


def parse_diamond_hits(
    raw_path: Path,
    proteins_path: Path,
    output_path: Path,
    *,
    batch_size: int = 100_000,
) -> int:
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    proteins = pq.read_table(proteins_path, columns=["protein_id", "species_id"]).to_pylist()
    species_by_protein = {int(item["protein_id"]): int(item["species_id"]) for item in proteins}

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(output_path.suffix + ".tmp")
    rows: list[HitRow] = []
    hit_count = 0
    writer: pq.ParquetWriter | None = None
    try:
        writer = pq.ParquetWriter(temporary_path, HIT_SCHEMA, compression="zstd")
        with raw_path.open(encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.rstrip("\n")
                if not line.strip():
                    continue
                fields = line.split("\t")
                if len(fields) != 8:
                    raise SearchError(
                        f"DIAMOND row {line_number} has {len(fields)} fields; expected 8"
                    )
                (
                    query_raw,
                    target_raw,
                    identity,
                    aligned,
                    query_length,
                    target_length,
                    evalue,
                    bitscore,
                ) = fields
                query_id = _protein_id(query_raw)
                target_id = _protein_id(target_raw)
                if query_id not in species_by_protein or target_id not in species_by_protein:
                    raise SearchError(f"DIAMOND row {line_number} references an unknown protein")
                try:
                    aligned_length = float(aligned)
                    query_length_value = float(query_length)
                    target_length_value = float(target_length)
                    if query_length_value <= 0 or target_length_value <= 0:
                        raise ValueError("sequence length must be positive")
                    row: HitRow = (
                        query_id,
                        target_id,
                        species_by_protein[query_id],
                        species_by_protein[target_id],
                        float(bitscore),
                        float(identity),
                        100.0 * aligned_length / query_length_value,
                        100.0 * aligned_length / target_length_value,
                        float(evalue),
                    )
                except ValueError as error:
                    raise SearchError(
                        f"Invalid numeric value in DIAMOND row {line_number}: {error}"
                    ) from error
                rows.append(row)
                hit_count += 1
                if len(rows) >= batch_size:
                    writer.write_table(_hit_table(rows))
                    rows.clear()
        if rows or hit_count == 0:
            writer.write_table(_hit_table(rows))
        writer.close()
        writer = None
        temporary_path.replace(output_path)
    except OSError as error:
        raise SearchError(f"Failed to stream DIAMOND hits from {raw_path}: {error}") from error
    finally:
        if writer is not None:
            writer.close()
        temporary_path.unlink(missing_ok=True)
    return hit_count
