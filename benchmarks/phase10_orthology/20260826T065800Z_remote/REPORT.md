# Phase 10 hierarchy-driven orthology streaming regression

## Provenance

- Run ID: `20260826T065800Z_phase10_orthology`
- Base Git commit: `e46af3ed522d4d8b21870a23409791e4d08e1953`
- Deployment: `dirty=1`, captured at `2026-08-26T06:58:04Z`
- Algorithm: `hierarchy-cross-child-v1`
- Remote environment: project `ogprofiler` Conda environment
- Fixture: four disjoint child families, four species, 100 genes per species
- Chunk size: 4,096 candidates

This was a short remote integration/streaming regression rather than a Slurm
workload. The generator is preserved as
`benchmarks/phase10_orthology/generate_fixture.py`; the archived Parquet inputs
capture the exact executed fixture.

## Results

| Metric | Value |
|---|---:|
| Prepared proteins | 400 |
| Orthology-supporting nodes | 1 `POLYTOMY` |
| Candidate pairs | 60,000 |
| Output chunks | 15 |
| Compressed size | 39,642 bytes |
| Wall time | 0.44 s |
| Maximum RSS | 63,824 KiB |
| Output SHA-256 | `11aa5e8623b1c08b768a46a36b9f1b03fea2c3ba4556e625909d5dbba2d5b164` |

The expected count is six species pairs multiplied by 10,000 gene pairs per
species pair. Decompression found exactly 60,000 data rows and 60,000 unique
canonical protein pairs. Every row is cross-species and satisfies
`protein_a_id < protein_b_id`.

The immediate second invocation checksum-verified and reused the compressed
output with the same digest. The default disabled invocation and damaged-output
reconstruction are covered by integration tests.

## Validation

- archived regression invocation: 58 tests passed;
- final local and remote Conda suites after the all-singleton boundary test:
  59 tests passed;
- Ruff and mypy (53 source files) passed;
- the 100,000-pair unit scale test remained below 2 MB peak Python allocation;
- Phase 0 fixture checksum verification passed.
