# Phase 5 connected-component engine

The production component stage indexes every prepared protein, streams the
canonical retained edge table through a Union-Find data structure, and writes
component-local edge partitions without constructing a global igraph object.

```bash
ogprofiler components --run run/
```

## Deterministic identity

Union-Find uses integer protein IDs, path compression, and union by size.
Component IDs do not depend on Union-Find root identity or edge order. After
the edge scan, components are sorted by descending vertex count and then by
their smallest protein ID. The largest component is therefore `0`, with stable
tie-breaking.

## Artifacts

```text
components/
├── index.parquet
├── statistics.parquet
├── singleton_terminal_families.parquet
├── component-manifest.json
└── edges/
    ├── component=00000000/part-*.parquet
    └── component=00000001/part-*.parquet
```

`index.parquet` maps `protein_id` to `component_id`. Statistics contain
`component_id`, `n_vertices`, `n_edges`, and `n_species`. Isolated proteins are
present in the index and are published once as terminal families with reason
`SINGLETON`; they have no edge partition and bypass Leiden in Phase 6.

`read_component_edges()` addresses one Hive-style partition directly.
`load_component_edge_table()` joins that partition with only the selected
component's protein IDs and returns the bounded component-local graph input
used by the hierarchy worker.

## Streaming and resume

The retained edge Parquet file is scanned in bounded record batches. The first
pass performs unions, a second counts edges after final component assignment,
and a third writes partitioned Parquet fragments with a bounded number of open
files. Edge rows are never represented as a global Python edge list.

The manifest records algorithm and parameter identities, input hashes, counts,
and every output fragment hash. An identical verified artifact set is reused;
changed input, configuration, or artifact content triggers reconstruction.
