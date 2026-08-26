"""Streaming connected-component indexing and partitioned edge access."""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path

import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

from ogprofiler.exceptions import ComponentError
from ogprofiler.graph.edges import EdgeTable
from ogprofiler.graph.union_find import UnionFind
from ogprofiler.similarity.io import RETAINED_EDGE_SCHEMA

COMPONENT_INDEX_SCHEMA = pa.schema(
    [("protein_id", pa.int64()), ("component_id", pa.int64())]
)
COMPONENT_STATISTICS_SCHEMA = pa.schema(
    [
        ("component_id", pa.int64()),
        ("n_vertices", pa.int64()),
        ("n_edges", pa.int64()),
        ("n_species", pa.int64()),
    ]
)
SINGLETON_SCHEMA = pa.schema(
    [
        ("protein_id", pa.int64()),
        ("component_id", pa.int64()),
        ("terminal_reason", pa.string()),
    ]
)
PARQUET_ROW_GROUP_SIZE = 262_144


@dataclass(frozen=True, slots=True)
class ComponentBuildResult:
    component_count: int
    singleton_count: int
    protein_count: int
    edge_count: int


def _protein_metadata(path: Path) -> tuple[list[int], dict[int, int]]:
    parquet = pq.ParquetFile(path)
    ids: list[int] = []
    species_by_id: dict[int, int] = {}
    for batch in parquet.iter_batches(columns=["protein_id", "species_id"]):
        for row in batch.to_pylist():
            protein_id = int(row["protein_id"])
            if protein_id in species_by_id:
                raise ComponentError(f"Duplicate protein ID: {protein_id}")
            ids.append(protein_id)
            species_by_id[protein_id] = int(row["species_id"])
    if not ids:
        raise ComponentError("Protein index is empty")
    return ids, species_by_id


def _validate_edge_schema(path: Path) -> pq.ParquetFile:
    parquet = pq.ParquetFile(path)
    if parquet.schema_arrow != RETAINED_EDGE_SCHEMA:
        raise ComponentError(f"Retained edge schema mismatch in {path}")
    return parquet


def build_component_artifacts(
    proteins_path: Path,
    edges_path: Path,
    output_directory: Path,
    *,
    batch_size: int = 65_536,
    max_open_files: int = 64,
) -> ComponentBuildResult:
    """Scan edges twice: first for DSU, then for bounded-memory partition output."""
    protein_ids, species_by_id = _protein_metadata(proteins_path)
    edge_file = _validate_edge_schema(edges_path)
    union_find = UnionFind(protein_ids)
    edge_count = 0
    for batch in edge_file.iter_batches(columns=["u", "v"], batch_size=batch_size):
        for left, right in zip(
            batch.column(0).to_pylist(), batch.column(1).to_pylist(), strict=True
        ):
            union_find.union(int(left), int(right))
            edge_count += 1

    vertices_by_root: dict[int, list[int]] = {}
    for protein_id in protein_ids:
        vertices_by_root.setdefault(union_find.find(protein_id), []).append(protein_id)
    ordered = sorted(
        (sorted(vertices) for vertices in vertices_by_root.values()),
        key=lambda vertices: (-len(vertices), vertices[0]),
    )
    component_by_protein: dict[int, int] = {}
    for component_id, vertices in enumerate(ordered):
        for protein_id in vertices:
            component_by_protein[protein_id] = component_id

    edge_counts = [0] * len(ordered)
    for batch in edge_file.iter_batches(columns=["u"], batch_size=batch_size):
        for left in batch.column(0).to_pylist():
            edge_counts[component_by_protein[int(left)]] += 1

    output_directory.mkdir(parents=True, exist_ok=True)
    build_root = output_directory / ".building"
    if build_root.exists():
        shutil.rmtree(build_root)
    edge_root = build_root / "edges"
    edge_root.mkdir(parents=True)

    index_rows = [
        {"protein_id": protein_id, "component_id": component_by_protein[protein_id]}
        for protein_id in sorted(protein_ids)
    ]
    stats_rows = []
    singleton_rows = []
    for component_id, vertices in enumerate(ordered):
        species = {species_by_id[protein_id] for protein_id in vertices}
        stats_rows.append(
            {
                "component_id": component_id,
                "n_vertices": len(vertices),
                "n_edges": edge_counts[component_id],
                "n_species": len(species),
            }
        )
        if len(vertices) == 1:
            singleton_rows.append(
                {
                    "protein_id": vertices[0],
                    "component_id": component_id,
                    "terminal_reason": "SINGLETON",
                }
            )
    pq.write_table(
        pa.Table.from_pylist(index_rows, schema=COMPONENT_INDEX_SCHEMA),
        build_root / "index.parquet",
        compression="zstd",
    )
    pq.write_table(
        pa.Table.from_pylist(stats_rows, schema=COMPONENT_STATISTICS_SCHEMA),
        build_root / "statistics.parquet",
        compression="zstd",
    )
    pq.write_table(
        pa.Table.from_pylist(singleton_rows, schema=SINGLETON_SCHEMA),
        build_root / "singleton_terminal_families.parquet",
        compression="zstd",
    )

    partition_schema = pa.schema([("component", pa.string())])
    partitioning = ds.partitioning(partition_schema, flavor="hive")
    edge_file = _validate_edge_schema(edges_path)
    for batch_number, batch in enumerate(edge_file.iter_batches(batch_size=batch_size)):
        components = [
            f"{component_by_protein[int(protein_id)]:08d}"
            for protein_id in batch.column(batch.schema.get_field_index("u")).to_pylist()
        ]
        table = pa.Table.from_batches([batch]).append_column(
            "component", pa.array(components, type=pa.string())
        )
        ds.write_dataset(
            table,
            edge_root,
            format=ds.ParquetFileFormat(),
            partitioning=partitioning,
            basename_template=f"part-{batch_number:08d}-{{i}}.parquet",
            existing_data_behavior="overwrite_or_ignore",
            max_open_files=max_open_files,
            file_options=ds.ParquetFileFormat().make_write_options(compression="zstd"),
            min_rows_per_group=PARQUET_ROW_GROUP_SIZE,
            max_rows_per_group=PARQUET_ROW_GROUP_SIZE,
        )

    for name in (
        "index.parquet",
        "statistics.parquet",
        "singleton_terminal_families.parquet",
        "edges",
    ):
        destination = output_directory / name
        if destination.is_dir():
            shutil.rmtree(destination)
        elif destination.exists():
            destination.unlink()
        (build_root / name).replace(destination)
    build_root.rmdir()
    return ComponentBuildResult(
        component_count=len(ordered),
        singleton_count=len(singleton_rows),
        protein_count=len(protein_ids),
        edge_count=edge_count,
    )


def read_component_edges(
    output_directory: Path, component_id: int, *, memory_map: bool = False
) -> pa.Table:
    """Read one edge partition without scanning partitions for other components."""
    path = output_directory / "edges" / f"component={component_id:08d}"
    if not path.is_dir():
        return pa.Table.from_pylist([], schema=RETAINED_EDGE_SCHEMA)
    if memory_map:
        fragments = [
            pq.read_table(fragment, schema=RETAINED_EDGE_SCHEMA, memory_map=True)
            for fragment in sorted(path.glob("*.parquet"))
        ]
        table = (
            pa.concat_tables(fragments)
            if fragments
            else pa.Table.from_pylist([], schema=RETAINED_EDGE_SCHEMA)
        )
    else:
        table = ds.dataset(path, format="parquet", schema=RETAINED_EDGE_SCHEMA).to_table()
    return table.sort_by([("u", "ascending"), ("v", "ascending")])


def load_component_edge_table(
    output_directory: Path, component_id: int, *, memory_map: bool = False
) -> EdgeTable:
    """Load one component and remap only when a later hierarchy worker requests it."""
    index = pq.read_table(
        output_directory / "index.parquet",
        filters=[("component_id", "=", component_id)],
        memory_map=memory_map,
    )
    vertices = [int(value) for value in index["protein_id"].to_pylist()]
    if not vertices:
        raise ComponentError(f"Unknown component ID: {component_id}")
    rows = read_component_edges(
        output_directory, component_id, memory_map=memory_map
    ).to_pylist()
    return EdgeTable.canonicalize(
        vertices,
        [(int(row["u"]), int(row["v"]), float(row["weight"])) for row in rows],
    )


def non_singleton_component_ids(output_directory: Path) -> tuple[int, ...]:
    """Return only components eligible for hierarchy scheduling, largest first."""
    rows = pq.read_table(
        output_directory / "statistics.parquet",
        columns=["component_id", "n_vertices"],
    ).to_pylist()
    return tuple(int(row["component_id"]) for row in rows if int(row["n_vertices"]) > 1)
