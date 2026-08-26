# Phase 12 scientific benchmark and algorithm selection

Phase 12 separates scientific selection from engineering performance. A run is
ranked from predeclared family, evolutionary-consistency, and orthology metrics;
runtime remains visible but is not allowed to replace scientific evidence.

## Parameter matrix

`ogprofiler benchmark plan --out OUT` writes a deterministic, checksummed
one-factor-at-a-time matrix. The baseline is compared with every non-baseline
value of these closed axes:

- normalization: legacy NBS, raw bitscore, length-scaled bitscore;
- directional coverage: 0%, 50%, 70%;
- symmetrization: maximum, minimum, arithmetic mean, geometric mean;
- Leiden objective: RBConfiguration, RBER, CPM, modularity;
- gamma search: adaptive and fixed log grid;
- split acceptance: permissive, default, and conservative;
- random seeds: 7, 42, and 104729.

The current matrix contains 16 rows including the baseline. OFAT makes each
effect attributable; interactions belong in a separately declared follow-up
matrix rather than being inferred from these runs.

## Metrics

`ogprofiler benchmark metrics` emits `scientific-benchmark-v1` JSON with input
checksums and the following dimensions:

- family count, size distribution, species and membership coverage;
- pairwise clustering precision/recall/F1 and exact-family recovery;
- complete single-copy-family recovery and large-family fragmentation;
- species-overlap consistency for network events;
- gene-tree concordance and reconciliation consistency when Phase 11 evidence
  exists;
- cross-species orthology precision/recall/F1 over complete single-copy
  reference families.

Compressed orthology output is decoded incrementally. The evaluator retains
only the evaluated pair set, not the compressed file or all parsed records.

Frozen V1 normalized membership can be evaluated through
`--legacy-normalized`. `ogprofiler benchmark compare` puts V1, V2, and imported
method results into one closed comparison schema. OrthoFinder and SonicParanoid
have explicit planned-import registry entries; their adapters are outside the
Phase 12 execution core.

## Formal execution

`slurm/phase12_scientific_matrix.sh` runs the 16 rows as a Slurm array capped at
two concurrent tasks. Each task requests 4 CPUs, 16 GiB, and 30 minutes and
executes an immutable source snapshot. Shared Dataset A DIAMOND hits isolate the
scientific effects downstream of search. Every row records input checksums,
resolved overrides, stage timings, and headline metrics.

After completion:

```bash
python benchmarks/aggregate_phase12_matrix.py \
  --root RUN_ROOT/results \
  --out RUN_ROOT/results/summary
```

The aggregator publishes every completed row, the unchanged baseline, and the
highest scientific score. A default change is accepted only after the complete
matrix succeeds and the selected row improves the declared scientific score
without an unexplained reference regression.

## Preflight evidence

The versioned `benchmarks/phase12_scientific/preflight/` artifacts establish
schema and regression readiness before formal submission:

- all 16 matrix rows are materialized;
- a real Dataset A tracer exactly recovers its 18 reference families;
- reference C/E exports and frozen V1 membership use the same family schema;
- C recovers 5/6 exact families and all 5 complete single-copy families;
- E exposes the expected unresolved large-family fragmentation signal;
- species-overlap annotations are internally consistent on C and E.

These are preflight observations, not the final parameter recommendation. The
formal Slurm matrix supplies the data-driven default decision.

## Completed Phase 12 decision

Formal Slurm job `1404184` completed all 16 rows with successful exit status.
Every row achieved family pair F1, exact-family recovery, and single-copy-family
recovery of 1.0 on Dataset A. Since no tested change improved the scientific
score, the selected action is to retain the baseline (`legacy_nbs`, 50%
directional coverage, maximum symmetrization, RBER, adaptive gamma search,
minimum family size 2, maximum child fraction 0.95, seed 42). Runtime from this
small concurrent run is not used to break the scientific tie. The archived
report and machine-readable evidence are under
`benchmarks/phase12_scientific/20260826T075642Z_b2f54a90/`.
