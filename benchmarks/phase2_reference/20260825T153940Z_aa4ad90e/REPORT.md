# Phase 2 final reference-SSN regression and scale benchmark

## Provenance

- Run ID: `20260825T153940Z_e46af3ed522d_aa4ad90e_67933`
- Slurm job: `1404170`
- Source SHA-256: `aa4ad90ecfb1327e4f23277d1f3641ce2860dc3a7ca548eb816ff4c12bfae9da`
- Three repeats per dataset; frozen V1 and V2 consume the identical reference SSN.
- Resources per repeat: 8 CPUs, 16 GiB, one hour.

## Scale and performance

| Dataset | Nodes | Edges | Components | Largest | V1 wall s | V2 wall s | Wall V1/V2 | V1 core s | V2 core s | Core V1/V2 | V1 RSS MiB | V2 RSS MiB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| C_gene_family_expansion | 80 | 602 | 6 | 60 | 5.320 | 1.690 | 3.148 | 3.9817 | 0.0338 | 117.849 | 108.1 | 59.2 |
| E_large_connected_component | 240 | 4685 | 1 | 240 | 14.590 | 0.830 | 17.578 | 13.3119 | 0.0353 | 376.875 | 111.4 | 63.0 |

Ratios above 1 favor V2. Wall values include command startup and artifact writing; core values are implementation-instrumented hierarchy computation.

## Regression and topology

| Dataset | V1 deterministic | V2 deterministic | Common IDs | Pair Jaccard | Node Jaccard | Edge Jaccard | Leiden calls | Candidates | Subgraphs | V2 nodes | V2 terminals |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| C_gene_family_expansion | false | true | 80 | 1.000 | 0.900 | 0.167 | 81 | 81 | 3 | 9 | 8 |
| E_large_connected_component | true | true | 240 | 1.000 | 0.636 | 0.000 | 61 | 61 | 6 | 7 | 6 |

## Phase 2 exit gate

### C_gene_family_expansion

- [x] three completed repeats
- [x] V2 normalized topology deterministic
- [x] V2 median RSS below V1
- [x] largest reference component exercised

### E_large_connected_component

- [x] three completed repeats
- [x] V2 normalized topology deterministic
- [x] V2 median RSS below V1
- [x] largest reference component exercised

Successful V2 command completion also means component-level hierarchy invariants passed before each artifact was written. Membership and topology differences remain explicit compatibility measurements rather than an exact-parity claim.

C produced two frozen-V1 semantic topology fingerprints across three repeats, while V2 produced one. E was deterministic in both implementations. Candidate counts in this run equal Leiden calls because the prototype evaluates exactly one candidate per call; the post-benchmark implementation additionally persists `candidates.parquet`.
