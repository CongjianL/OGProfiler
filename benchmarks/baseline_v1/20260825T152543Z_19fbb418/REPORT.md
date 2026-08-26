# Frozen V1 A–E baseline and same-SSN V1/V2 comparison

## Run identity

- Run ID: `20260825T152543Z_e46af3ed522d_19fbb418_53770`
- Slurm job: `1404165` (array tasks 0–4: `COMPLETED`)
- Source SHA-256: `19fbb418d6de4c0a59032e843de154d84ee5147dcf173b18a8a66dff20e024ba`
- Dataset tree SHA-256: `e4abe89f62f28fa4351cffc453d66c416719ecc1d4b8a72165e89c2f013a13d5`
- Resources per task: 8 CPUs, 16 GiB, 2 hours; maximum two concurrent tasks.
- V1 full pipeline: DIAMOND search → SSN → hierarchy → OG output.
- V1/V2 hierarchy comparison: both implementations consumed the exact V1-produced `ssn.gml`.

## Frozen V1 full-pipeline baseline

| Dataset | Input proteins | SSN nodes | SSN edges | Components | Largest | Wall s | RSS MiB |
|---|---:|---:|---:|---:|---:|---:|---:|
| A_small_sanity | 48 | 48 | 60 | 18 | 4 | 14.96 | 112.6 |
| B_paralog | 48 | 48 | 59 | 19 | 6 | 16.36 | 112.7 |
| C_gene_family_expansion | 80 | 80 | 23 | 64 | 3 | 12.10 | 112.8 |
| D_fusion_multidomain | 48 | 48 | 66 | 15 | 4 | 11.52 | 112.6 |
| E_large_connected_component | 240 | 240 | 136 | 144 | 19 | 12.16 | 112.8 |

## Same-SSN hierarchy runtime and RSS

| Dataset | V1 wall s | V2 wall s | Wall V1/V2 | V1 core s | V2 core s | Core V1/V2 | V1 RSS MiB | V2 RSS MiB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A_small_sanity | 1.96 | 4.37 | 0.45 | 0.4550 | 0.0326 | 13.96 | 108.3 | 58.8 |
| B_paralog | 2.06 | 4.67 | 0.44 | 0.4725 | 0.0306 | 15.42 | 108.3 | 58.6 |
| C_gene_family_expansion | 1.47 | 5.40 | 0.27 | 0.1418 | 0.0101 | 14.00 | 108.2 | 59.2 |
| D_fusion_multidomain | 1.54 | 2.28 | 0.68 | 0.2492 | 0.0323 | 7.72 | 108.2 | 59.3 |
| E_large_connected_component | 1.52 | 9.58 | 0.16 | 0.1967 | 0.0666 | 2.95 | 108.3 | 59.6 |

`wall` is GNU `/usr/bin/time -v` around the complete hierarchy command, including environment startup and artifact serialization. `core` is implementation-instrumented hierarchy computation. Ratios above 1 favor V2.

## Membership and hierarchy topology

| Dataset | Common members | Pair precision | Pair recall | Pair Jaccard | Node Jaccard | Edge Jaccard |
|---|---:|---:|---:|---:|---:|---:|
| A_small_sanity | 40 | 0.000 | NA | 0.000 | 0.128 | 0.000 |
| B_paralog | 38 | 0.000 | NA | 0.000 | 0.117 | 0.000 |
| C_gene_family_expansion | 24 | 0.000 | NA | 0.000 | 0.083 | 0.000 |
| D_fusion_multidomain | 44 | 0.000 | NA | 0.000 | 0.136 | 0.000 |
| E_large_connected_component | 107 | 0.790 | 0.861 | 0.701 | 0.154 | 0.250 |

Membership metrics compare co-clustered protein pairs over IDs present in both normalized hierarchies. Topology metrics compare descendant-member-set identities and parent→child edges, so transient numeric cluster IDs do not affect the result.

## Result reading

- The V2 core hierarchy routine used less memory and less instrumented compute time in all five cases.
- End-to-end hierarchy wall time was higher for V2 in all five cases because the prototype command starts a fresh environment and serializes per-component Parquet artifacts.
- A–D show structural divergence: the V1 hierarchy produces nested singleton terminal leaves, whereas V2 terminates these small components at their component roots.
- E exercises recursive splitting and is the informative topology case: pair Jaccard 0.701, node Jaccard 0.154, and edge Jaccard 0.250.
- The V1 hierarchy membership export excludes SSN isolates; V2 records them as singleton components. Pairwise metrics therefore use only the common ID set.

The Slurm state establishes successful execution. The reported membership/topology values establish measured prototype behavior; they do not yet establish V2 parity with frozen V1.
