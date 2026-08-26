# Phase 6 production Hierarchical Leiden reference benchmark

## Provenance

- Slurm array job: `1404177` (`0=C`, `1=E`)
- Run ID: `20260826T023943Z_e46af3ed522d_484dc7c7_81066`
- Immutable source SHA-256: `484dc7c78ac448798b33ddeb2f6554d94502dada05d0ca74d71b689ee73ea75a`
- Git commit: `e46af3ed522d4d8b21870a23409791e4d08e1953` (`dirty=1`, captured snapshot)
- Resources per task: 4 CPUs, 16 GiB, 30-minute limit
- Both array tasks: `COMPLETED`, exit `0:0`, elapsed 10 seconds
- Stability mode: `robust` (three deterministic seeds per candidate)

Both runs consume the frozen reference SSN through the fixed retained-edge
schema. The hierarchy stage reads only component 0 from its Phase 5 partition.

## Results

| Dataset/component | Proteins | Edges | Nodes | Terminal families | Candidates | Leiden calls | Subgraphs | Core s | Wall s | Peak RSS MiB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| C largest component | 60 | 572 | 4 | 3 | 31 | 93 | 3 | 0.0761 | 1.26 | 71.7 |
| E component 0 | 240 | 4,685 | 7 | 6 | 61 | 183 | 6 | 0.1399 | 1.05 | 75.8 |

The C fixture contains 80 prepared proteins in six Phase 5 components; this
Phase 6 core benchmark deliberately processes only its 60-protein largest
component. Phase 7 will schedule the other components. E is one 240-protein
component and therefore exercises the complete fixture.

All C/E terminal memberships match the planted three- and six-family
structures. Robust mode accounts exactly for three Leiden calls per resolution
candidate. Only one induced subgraph is constructed per child; gamma and seed
searches reuse that cluster graph.

Peak process RSS rises by only 4.1 MiB from the 60-protein C component to the
240-protein E component. On identical reference E input, the earlier frozen-V1
median was 14.59 s wall, 13.31 s hierarchy core, and 111.4 MiB RSS. Phase 6
records 1.05 s wall, 0.140 s core, and 75.8 MiB RSS while additionally running
three-seed stability checks.

## Regression gates

- hierarchy invariants validated before artifact publication;
- k-way roots recovered as 3-way C and 6-way E splits;
- terminal membership covers every protein in the selected component once;
- full candidate traces include acceptance, stability, ARI/NMI, tiny-fragment,
  quality, and edge-separation fields;
- second invocation reused verified artifacts without changing `nodes.parquet`;
- downloaded node/member/candidate/metrics hashes match each hierarchy manifest.

Key artifact identities:

- C members: `727da690a74438b7f8f6071f501f51298ddaff2244d5d2467cccb656b6964fca`
- E members: `343572b64fa1f4212bbca3a37e209831c27a14398d7bc42cda4887673faf5dc8`
- C manifest: `c6ec7229f9c41f6c78138e83c4e22a24c0e38828730ec50bff9f99bff8c94ccb`
- E manifest: `f994155f91fa137b32c9322092a598dc90dc4ab4f19283184633b544eddfa98b`
