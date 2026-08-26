# Changelog

## 2.0.0a1 - 2026-08-26

### Phase 14

- Added Arrow memory-mapped component reads, NumPy memmapped global/species
  vectors, and binary-search mapping views to remove per-protein dictionary
  storage from hierarchy workers.
- Added deterministic path-based subtree identities and checksummed,
  largest-first component task manifests. Accepted splits now recursively
  release child process tasks and merge completion-order-independent topology.
- Added Parquet compression/row-group sequential/random I/O profiling plus a
  two-shard Slurm runner with component-level verified resume.

### Phase 15

- Added `ogprofiler run`, `status`, and `inspect`, including dry-run planning,
  stage-bounded resume, JSON reporting, and component summaries.
- Added runtime/external-tool provenance capture and comprehensive installation,
  quick-start, output, scientific-caveat, reproducibility, and release docs.

### Phases 0-13

- Freeze the OGProfiler 1 reference implementation.
- Add the OGProfiler 2 package, configuration, workspace, core models, and
  deterministic FASTA prepare stage.
- Accept the component-centric hierarchy architecture ADR.
- Complete the Phase 2 C/E reference-SSN regression and scale benchmark.
- Persist resolution-search candidate traces and strengthen hierarchy
  invariants for component identity, depth, and descendant cardinality.
- Add the Phase 3 `SearchBackend` protocol and DIAMOND all-vs-all backend with
  standard directional Parquet hits, checksummed provenance, and verified
  resume.
- Add the Phase 4 numeric edge engine: fixed batched hit schema, legacy NBS,
  coverage filtering, best-hit/RBH/LRB retention, AR/ARB compatibility,
  canonical direction joins, four symmetrization policies, and verified edge
  manifests.
- Add the Phase 5 streaming Union-Find component engine with deterministic
  largest-first IDs, component statistics, singleton terminal-family handoff,
  Hive-style edge partitions, component-local loading, and verified resume.
- Promote the component-centric prototype to the Phase 6 production
  Hierarchical Leiden core with component-local loading, explicit split
  acceptance, multi-seed stability/ARI/NMI, complete terminal reasons,
  species bitmaps, DFS graph lifetime control, and checksummed component stores.
- Add the Phase 7 largest-first process scheduler, scheduler-owned SQLite task
  state and audit events, stale-run recovery, configurable retries,
  algorithm/config invalidation, atomic component publication, fault isolation,
  failed-only execution, and checksummed resume.
- Add Phase 8 species-bitmap network evolution annotations for binary and k-way
  hierarchy nodes, continuous overlap evidence, confidence, closed event
  taxonomy, legacy I/II/III export mapping, and component-level verified resume.
- Add Phase 9 deterministic terminal-family IDs and final families, members,
  hierarchy, and event TSV exports, including singleton handoff, verified
  manifests, selected-family FASTA, and opt-in component GraphML.
- Add the Phase 10 hierarchy-driven orthology engine with speciation-like and
  disjoint-polytomy traversal, explicit co-ortholog candidate semantics,
  same-species filtering, bounded-memory chunked Zstandard output, progress,
  checksummed resume, and opt-in pair generation.
- Add Phase 11 selected-family phylogenetic refinement with deterministic
  selection policy, MAFFT and FastTree adapters, Newick/species-tree support,
  midpoint/outgroup/species-aware rooting, LCA reconciliation, and independent
  network/phylogenetic event evidence with conflict status.
- Add the Phase 12 scientific benchmark framework with a deterministic 16-row
  parameter matrix, configurable normalization and gamma-search strategies,
  family/evolution/orthology metrics, frozen-V1 evaluation, closed method
  comparison, and reproducible Slurm matrix aggregation.
- Add the Phase 13 deterministic synthetic evolution framework with a
  full-factorial scenario matrix, explicit genealogy/event/orthology/domain
  truth, recovery metrics, applicability thresholds, and strided Slurm runner.
