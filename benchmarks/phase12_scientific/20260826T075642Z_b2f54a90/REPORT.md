# Phase 12 scientific parameter benchmark

## Execution identity

- Slurm job: `1404184` (`0-15%2`)
- Run ID: `20260826T075642Z_e46af3ed522d_b2f54a90_94464`
- Source SHA-256: `b2f54a90928376ef7ab1a960083cd4253df080ec818f0a1cce088254f935b68f`
- Git commit: `e46af3ed522d4d8b21870a23409791e4d08e1953`, exact dirty state snapshotted
- Dataset: `A_small_sanity`, 48 proteins in four species
- Search input: frozen DIAMOND hits, SHA-256
  `2ce33504812e26d0da23614124b7b0cc7dd7229e16b0e4b5b09cd26d2e8ac066`
- Resources: 16 array tasks, at most two concurrent; 4 CPU, 16 GiB, 30 minutes per task

## Execution audit

All 16 array tasks reached `COMPLETED` with exit code `0:0`. There are 16
scientific metric documents, 16 run summaries, and 16 verified export
manifests. All 16 stdout and 16 stderr files were inspected; stderr contains
normal structured stage logging and no `Traceback`, `ERROR`, `FAILED`, or
`Exception` marker.

Observed task elapsed time was 6–11 seconds. Pipeline time recorded inside each
row ranged from 5.995 to 8.629 seconds (median 6.963 seconds). These timings are
reported for execution context only: the tiny two-at-a-time run is not a
performance-ranking experiment.

## Scientific result

Every baseline and non-baseline row produced the same perfect Dataset A family
result:

| Metric | All 16 rows |
|---|---:|
| Pairwise family F1 | 1.0 |
| Exact-family recovery | 1.0 |
| Complete single-copy-family recovery | 1.0 |
| Scientific score | 1.0 |

The tested changes to normalization, directional coverage, symmetrization,
Leiden objective, gamma search, split acceptance, and seed produced no Dataset
A regression and no scientific improvement over baseline. Dataset A contains
no evaluable internal network-event node in these terminal outputs, so
species-overlap consistency is `NA`. Pairwise ortholog generation was outside
this matrix, so this run supplies no orthology ranking signal.

## Default parameter decision

Retain the existing baseline as the Phase 12 default:

```yaml
similarity:
  normalization: legacy_nbs
edges:
  min_query_coverage: 50.0
  min_target_coverage: 50.0
  symmetrization: max
hierarchy:
  method: rber
  resolution_strategy: adaptive
  min_family_size: 2
  max_child_fraction: 0.95
  seed: 42
```

This is a conservative evidence-backed decision: the baseline exactly recovers
all Dataset A reference families, remains stable across the tested seed and
algorithm axes, preserves V1-compatible normalization, and no competing row
improves the declared scientific score. Runtime is not used as a tie-breaker.

The decision is scoped to the current deterministic reference suite. Reference
C/E already expose expansion and fragmentation behavior, while broader
evolution-event and orthology discrimination requires Phase 13 ground-truth
genealogies. A future default change should therefore require a declared
multi-dataset improvement rather than a Dataset A timing difference.
