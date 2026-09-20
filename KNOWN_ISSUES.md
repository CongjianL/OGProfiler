# OGProfiler 1 known issues

This inventory freezes the architectural and correctness risks of the legacy
implementation. Each item is a regression target or an explicit behavior
change to classify during OGProfiler 2 development.

| Area | Known issue | V2 disposition |
| --- | --- | --- |
| Architecture | The pipeline is coupled in one large script. | Replace with explicit, layered modules. |
| Hierarchy | Recursive hierarchy construction repeatedly calls `subgraph()`. | Materialize only the active component/cluster representation. |
| Resolution | Resolution search repeatedly calls `partition.subgraphs()`. | Search over memberships and bounded candidate metrics. |
| Parallelism | Pools and `Manager.Queue` objects are repeatedly created. | Use a persistent component-level scheduler. |
| IPC | Complete graph objects are pickled between workers. | Pass `component_id`; load component edges inside the worker. |
| Storage | Internal hierarchy nodes duplicate complete descendant gene lists. | Store membership once at terminal families. |
| Orthology | Pair combinations plus repeated shortest paths have high complexity. | Traverse hierarchy and stream pairs. |
| Serialization | CM dumps in `sorted_connected_rbh`/`sorted_connected_rh` have suspect indentation. | Cover with focused V1 regression fixtures. |
| Reciprocal hits | Directional hit loading/indexing can disagree. | Canonicalize directed hits before reciprocal joins. |
| Edge weights | Undirected edge weights can depend on input direction. | Define and test an order-invariant symmetrization policy. |
| Determinism | `os.listdir()` order influences generated IDs and execution order. | Sort normalized input paths before ID assignment. |
| Determinism | Random seeds are not centrally controlled. | Resolve and record one run seed. |
| External tools | FastTree has a hard-coded path. | Discover binaries through configuration and `PATH`. |
| Orchestration | Pipeline stages are dispatched through `eval()` and `methodcaller()`. | Use explicit function calls. |
| Subprocesses | External commands use `shell=True` and shell redirection. | Use argument-vector `subprocess.run(..., check=True)`. |
| Pairwise output | `WriteOGPairwise` assumes gene IDs contain `|` (`split('|')[1]`) and crashes with `IndexError` when they do not (e.g. Prodigal IDs `NZ_XXX.1_1`). | Not reproduced in V2; V1 pairwise output is not part of the S0 comparison. |
| Similarity | Coverage is calculated but does not consistently filter retained hits. | Disabled by default (`edges.apply_coverage_filter=false`), matching V1/OrthoFinder3. |

## V1 ↔ V2 有意差异（公示）

These differences are intentional V2 behavior changes, not regression targets.
They are recorded here so V1↔V2 result differences can be attributed.

| Area | Difference | Disposition |
| --- | --- | --- |
| Input | V2 sorts input paths and gene IDs deterministically; V1 used `os.listdir()` / file order. | Keep V2 sorted order; do not reproduce V1 non-determinism. |
| Input | V2 uppercases and validates sequences; V1 kept raw sequences. | Keep V2 validation. |
| Hierarchy | V2 uses stability-based hierarchical Leiden; V1 used target-community-count gamma search. | Keep V2 algorithm. |
| Evolution | V2 annotates events via species bitmaps/overlap score; V1 used degree+neighbor genome sets (`I/II/III-*`). | Keep V2 algorithm. |
| OG export | V2 exports terminal families; V1 used coalescence extraction (`OGFile_coalescence_SameGenome.txt`). | Keep V2 export; V1 output format not reproduced. |
| `-d ar` | V1 `sorted_connected_rh` had a direction/indexing bug (kept all hits); V2 fixes it. | Keep V2 fix. |
| blastp multi-HSP | V1 kept the last HSP (matrix overwrite); V2/OrthoFinder3 keep the max score. | Keep V2 max dedup. |

