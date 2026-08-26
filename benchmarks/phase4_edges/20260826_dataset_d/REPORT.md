# Phase 4 Dataset D coverage sensitivity

Date: 2026-08-26  
Environment: remote `ogprofiler` Conda environment  
Search backend: DIAMOND 2.2.5 (`--sensitive --max-target-seqs 0 --max-hsps 1`)  
Dataset: `benchmarks/datasets/D_fusion_multidomain` (48 proteins, 4 species)

## Comparison

Both runs use legacy NBS, LRB retention, a best-hit tolerance of `1e-3`, and
`max` symmetrization. Only the directional query and target coverage thresholds
change.

| threshold | raw directional hits | normalized hits | coverage-eligible hits | retained directional hits | canonical edges |
|---:|---:|---:|---:|---:|---:|
| 0% | 244 | 196 | 196 | 132 | 66 |
| 50% | 244 | 196 | 132 | 132 | 66 |

The 50% cutoff removes 64 of 196 normalized directional candidates (32.7%).
Those candidates are not selected by the LRB threshold at 0%, so the final 66
canonical edges are byte-identical in this fixture. This is the expected useful
sensitivity result: coverage filtering has a measurable candidate-level effect,
while the retained graph is robust for the present Dataset D score structure.

## Reproducibility

Remote run directories:

- `/home/mselab/licj/projects/running/ogprofiler-runs/phase4-dataset-d-20260826-cov0`
- `/home/mselab/licj/projects/running/ogprofiler-runs/phase4-dataset-d-20260826-cov50`

Archived evidence is under `cov0/` and `cov50/`. Each directory contains the
search and edge manifests and the retained edge table; `cov0/hits.parquet`
preserves the common raw hit table used for local audit.

Key SHA-256 identities:

- common hit table: `4aa2e3c89459e643842b7e081d384317de4dff594a356981f4b8fd5b65c49535`
- common retained edge table: `a94a1b6e74414b4e5c91d8d0d2b6b87ce7a500050a508f3e2c3efd064e732245`
- 0% edge manifest: `48e8a40674a966c23c9dccbe930a1029442dcd833f15df779d4634b8e81dd1d1`
- 50% edge manifest: `432a41e1d0ff28a7ff8bbae94ae1974565084aefa1daf8bce3a7bfae7f45497c`

The remote integration suite passed `37` tests before the sensitivity run. This
fixture is an integration/regression check rather than a runtime scale benchmark.
