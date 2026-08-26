# Phase 5 reference E component regression and scale check

Date: 2026-08-26  
Environment: remote `ogprofiler` Conda environment  
Input: frozen Dataset E `reference_ssn.gml`  
Remote run: `/home/mselab/licj/projects/running/ogprofiler-runs/phase5-reference-e-20260826`

The reference SSN was converted to the fixed retained-edge schema so this run
measures the component engine independently of search and edge-retention
variation.

## Result

| metric | value |
|---|---:|
| proteins | 240 |
| retained edges | 4,685 |
| components | 1 |
| largest-component vertices | 240 |
| largest-component edges | 4,685 |
| species | 4 |
| singleton terminal families | 0 |
| partition fragments | 37 |
| component-stage wall time | 0.41 s |
| component-stage peak RSS | 69,804 KiB |

The partition loader recovered exactly 240 vertices and 4,685 edges from the
single requested component. The 37 fragments result from the deliberately
small 128-row benchmark batch size and verify multi-chunk local loading.

The full component artifact set contains 40 checksummed files totaling 123,217
bytes. Local verification checked every recorded output hash, the retained-edge
input hash, statistics, and the component-local round trip.

Key SHA-256 identities:

- retained edges: `2a17ff0c98421c6f124aed4947eebdd6196f2362ba5610c8c6d49cd364da66ae`
- component manifest: `bd699f74bef302e01cace85f9ef4b3ec16e1a33e1d31124c82143178fca2e75f`

The timing includes Python process startup and artifact hashing. Dataset E is a
regression-scale fixture rather than a production extrapolation.

The final local and remote suites both passed all 40 tests, including randomized
igraph equivalence, disconnected record batches, singleton handling,
component-local reads, verified resume, and corrupt-partition reconstruction.
