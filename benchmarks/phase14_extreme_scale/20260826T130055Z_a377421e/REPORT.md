# Phase 14 extreme-scale optimization benchmark

## Provenance

- Slurm array job: `1404222` (`0=C`, `1=E`)
- Run ID: `20260826T130055Z_c1da58cee4ee_a377421e_95874`
- Immutable source SHA-256: `a377421e044d837f7525c8f8cdbc122a2aa4f68892b21c7b6f4f9e38844342e4`
- Git commit: `c1da58cee4ee897eb3f30df0f4bf6c27bf171a8f` (`dirty=1`, captured snapshot)
- Resources per task: 4 CPUs, 16 GiB, 30-minute limit
- Both tasks: `COMPLETED`, exit `0:0`; elapsed 25 s (C) and 20 s (E)

## Recursive hierarchy and resume

| Dataset | Components | Nodes | Terminal families | Parallel wall s | Resume wall s | Parallel max RSS MiB |
|---|---:|---:|---:|---:|---:|---:|
| C | 6 | 9 | 8 | 6.58 | 1.13 | 72.6 |
| E | 1 | 7 | 6 | 3.22 | 0.63 | 75.9 |

Every accepted split was eligible for recursive process release (`release_size=1`).
The merged result indexes verified all seven component manifests and artifact
checksums. The immediate second pass reused all verified hierarchy outputs.

## Parquet I/O selection

Each task profiled nine layouts on a deterministic 2,000,000-row edge-shaped
table. Zstandard with 262,144-row groups was the smallest compressed layout
(25.9 MiB) and the fastest compressed layout for both sequential and randomized
reads in both runs. It is the Phase 14 default recommendation. Uncompressed
reads were faster, but required about 48.5 MiB and are not selected for durable
partition storage.

This run validates recursive scheduling, deterministic merge, mmap loading,
verified resume, and the I/O selection protocol on the frozen C/E regressions.
It is engineering evidence rather than a claim about a particular future
proteome-scale dataset.
