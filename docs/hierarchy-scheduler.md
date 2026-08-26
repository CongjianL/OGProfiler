# Phase 7 component scheduler and resume

Run every non-singleton Phase 5 component through the production hierarchy
worker pool with:

```bash
ogprofiler hierarchy-all --run run/ --set runtime.workers=4
```

To revisit only failed or invalidated components:

```bash
ogprofiler hierarchy-all --run run/ --failed-only
```

## Scheduling and isolation

Components are read from `components/statistics.parquet`, already ordered by
descending vertex count with deterministic component-ID tie-breaking. The main
process creates a spawn-based `ProcessPoolExecutor`; immutable run/config state
is installed by the worker initializer, and each submitted work item contains
only one integer `component_id`.

Each worker loads its own component partition, creates its own igraph, validates
the inferred hierarchy, and publishes directly to its component-local output.
The main process receives only a bounded `WorkerSummary`. There is no Manager,
shared graph, or large-object queue.

## Checkpoint state machine

The main process is the sole owner of `run.db` writes. A task row contains:

```text
stage, task_id, status, started_at, completed_at, input_hash,
output_path, error, attempts, algorithm_version
```

Transitions are appended to `task_events` for audit. Resume behavior is:

- `DONE` plus verified manifest and hashes: skip;
- stale `RUNNING`: reset to `PENDING`, then accept already-published verified
  output or dispatch again;
- `FAILED`: retry through `runtime.component_retries`;
- changed input, algorithm, or hierarchy configuration: record `INVALID` and
  rebuild;
- damaged or partial output: record `INVALID` and rebuild.

Workers stage all hierarchy files in a private temporary directory. Individual
artifacts are atomically replaced and the checksummed manifest is published
last. A process interruption therefore leaves either a verified result or an
artifact set that the next scheduler invocation detects and reconstructs.

## Fault isolation

A component exception is converted to a compact failure summary. Other futures
continue, the final scheduler manifest lists failed component IDs, and the CLI
returns a failing stage status only after all runnable components finish.
`--failed-only` then targets those tasks without recomputing verified work.

Singleton components remain in the Phase 5 singleton terminal-family artifact
and never enter the Leiden worker pool.
