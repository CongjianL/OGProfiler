# Phase 3 three-backend Dataset A smoke test

- Date: 2026-08-26
- Host/environment: `campus-server`, conda environment `ogprofiler`
- Source: commit `4cd85257c1213006623bf6cda1862c0ae6691c6c` plus the captured dirty backend implementation
- Dataset: `A_small_sanity`, 48 prepared proteins from four species
- Backends: DIAMOND 2.2.5, MMseqs2 18.8cc5c, BLAST+ 2.17.0+
- Mode: global directional all-vs-all search
- Result: every backend emitted 168 normalized hit rows
- Schema: all three outputs used the same nine-column directional-hit schema
- Execution class: short remote environment/integration smoke test; no Slurm job
- Remote output root: `/tmp/ogprofiler-search-backends-20260826T1431`

The run validates real executable discovery, database construction, search,
streaming parsing, and normalized Parquet publication for all three backends.
It is regression evidence rather than a runtime or sensitivity comparison.
