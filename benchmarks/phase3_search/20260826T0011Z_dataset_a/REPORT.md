# Phase 3 DIAMOND Dataset A smoke test

- Date: 2026-08-26
- Host/environment: `campus-server`, conda environment `ogprofiler`
- Dataset: `A_small_sanity`, 48 prepared proteins from four species
- Backend: DIAMOND 2.2.5
- Mode: global directional all-vs-all search
- Parameters: sensitive, e-value 0.001, two threads, unlimited target count
- Result: 168 directional hit rows in nine standard columns
- Parser: bounded streaming conversion to Parquet
- `hits.parquet` SHA-256:
  `2ce33504812e26d0da23614124b7b0cc7dd7229e16b0e4b5b09cd26d2e8ac066`
- Manifest/output checksum verification: passed
- Verified resume: identical invocation checked the current DIAMOND version and
  reused the search artifact without rebuilding the database or rerunning
  `blastp`
- Local and remote regression suites: 27 passed

This is a short environment/integration validation rather than a performance
benchmark. The captured Parquet and manifest are regression evidence for the
Phase 3 Diamond MVP path.

