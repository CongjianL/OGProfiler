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
| Similarity | Coverage is calculated but does not consistently filter retained hits. | Make coverage a validated edge-filter parameter. |

