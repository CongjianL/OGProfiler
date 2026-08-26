# Phase 8 network evolution annotation regression

## Provenance

- Source hierarchy: Phase 7 Slurm run `1404182`
- Source run ID: `20260826T060816Z_e46af3ed522d_55098951_14986`
- Source Git commit: `e46af3ed522d4d8b21870a23409791e4d08e1953`
- Annotation algorithm: `network-overlap-v1`
- Network overlap threshold: `0.5`
- Remote environment: project `ogprofiler` Conda environment
- Execution class: short remote integration regression on completed C/E
  hierarchy artifacts; no additional Slurm workload was needed

The Phase 8 implementation was deployed from a captured dirty development
state. The source commit above therefore identifies the base revision rather
than claiming that it contains the uncommitted Phase 8 changes.

## Results

| Dataset | Components | Nodes | Internal event | Child count | Overlap count | Overlap score | Confidence | Legacy label | Resume reused |
|---|---:|---:|---|---:|---:|---:|---:|---|---:|
| C | 6 | 9 | `DUPLICATION_LIKE` | 3 | 4 | 1.0 | 1.0 | `III-3` | 6/6 |
| E | 1 | 7 | `DUPLICATION_LIKE` | 6 | 4 | 1.0 | 1.0 | `III-3` | 1/1 |

C contains one three-way internal hierarchy node and eight multi-species
terminal nodes. E contains one six-way internal hierarchy node and six
multi-species terminal nodes. The terminal nodes are intentionally
`AMBIGUOUS`; this preserves the distinction between a terminal network family
and a one-species `SPECIES_SPECIFIC` family.

The root child species sets in both fixtures overlap across all four species,
so their normalized overlap is 1.0 and their k-way event is
`DUPLICATION_LIKE`. The compatibility export assigns `III-3` without driving
the internal classification.

## Verification

- all seven component event-table hashes match their component manifests;
- manifest node counts equal Parquet row counts (C: 9, E: 7);
- every hierarchy node has a closed-schema `network_event`;
- the immediate second invocation checksum-verified and reused all components;
- local and remote suites both passed: 52 tests;
- Ruff, mypy (46 source files), Phase 0 fixture checksums, and
  `git diff --check` passed before archiving.

The archive retains each component's `events.parquet` and
`network-event-manifest.json`, plus the global event manifest for both
references.
