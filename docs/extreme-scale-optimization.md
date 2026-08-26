# Phase 14 extreme-scale optimization

Phase 14 starts only after the scientific behavior frozen in Phases 12–13. It
changes storage and scheduling, not clustering semantics.

## Dense component-local loading

`ComponentGraphLoader` now memory-maps Arrow/Parquet inputs and materializes two
aligned dense vectors per component:

- `global_ids: int64[n]`
- `species_ids: int32[n]`

The vectors are atomically cached as NumPy `.npy` files and reopened read-only
with `mmap_mode="r"`. Global-ID lookup uses binary search over the sorted ID
vector instead of one Python dictionary entry per protein. The legacy mapping
interfaces remain read-only views so the hierarchy algorithm is unchanged.

## Deterministic subtree identity

Subtrees use `(component_id, child_ordinal_path)` identities. Children are
ordered by their minimum global protein ID, as in the frozen hierarchy engine.
Lexicographic path order is therefore deterministic depth-first order and is
converted to final integer cluster IDs during merge independent of worker
completion order. Once a split is accepted, every child is immediately
released to the process pool; further accepted splits recursively release their
children. The production stage activates this path when both
`hierarchy.subtree_workers > 1` and component size reaches
`hierarchy.subtree_release_size`.

## I/O profiling

`benchmarks/run_phase14_io.py` measures the Cartesian product of:

- compression: none, Snappy, Zstandard;
- row groups: 16,384, 65,536, 262,144 rows;
- access: full sequential memory-mapped scan and randomized row-group reads.

Raw timings must be interpreted on the target filesystem. The selected layout
must balance file size and both access modes rather than choosing from a single
local run.

Formal Slurm job `1404222` profiled two independent 2-million-row tables.
Zstandard with 262,144-row groups was simultaneously the smallest compressed
layout and the fastest compressed sequential/random layout in both tasks; it
is therefore the durable component-partition default.

## Slurm task protocol

`build_phase14_task_manifest.py` produces a checksummed, largest-first Parquet
manifest. `phase14_extreme_scale.sh` consumes it as two strided array shards,
each requesting 4 CPU, 16 GiB, and 30 minutes. Component execution retains the
existing checksummed hierarchy resume behavior, so retrying a failed shard
reuses verified component outputs. `aggregate_phase14_components.py` verifies
every component artifact checksum and publishes a merged result index; it is
the completeness barrier before downstream annotation/export.

## Distributed-runtime decision

Ray or Dask is deferred. The current scheduler already has deterministic task
identity, immutable Slurm boundaries, and verified resume. A distributed
runtime is justified only if scale measurements show that recursive subtree
coordination—not Leiden, graph construction, or filesystem throughput—is the
dominant bottleneck.

## Phase exit

Job `1404222` completed both C/E tasks with exit `0:0`, verified all seven
component result manifests, and reused every output during the second pass.
The formal evidence is archived in
`benchmarks/phase14_extreme_scale/20260826T130055Z_a377421e/`.
