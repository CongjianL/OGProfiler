# Phase 3 DIAMOND Dataset A smoke test

- Date: 2026-08-26
- Host/environment: `campus-server`, conda environment `ogprofiler`
- Dataset: `A_small_sanity`, 48 prepared proteins from four species
- Backend: DIAMOND 2.2.5
- Mode: global directional all-vs-all search
- Parameters: sensitive, e-value 0.001, two threads, unlimited target count
- Result: 168 directional hit rows
- Output schema: query/target protein and species IDs, bitscore, identity,
  query/target coverage, and e-value
- `hits.parquet` SHA-256:
  `30bc5bde6db22947a8ed452dc15856bfbe343b6656c16bd41b978b13b600f727`
- Manifest/output checksum verification: passed
- Second identical invocation: verified manifest/checksum resume; database build
  and search were skipped
- Remote regression suite: 27 passed

This is a short environment/integration validation rather than a performance
benchmark. The captured Parquet and manifest are regression evidence for the
Phase 3 Diamond MVP path.

