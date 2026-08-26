# Phase 6 production Hierarchical Leiden core

The production hierarchy command processes one Phase 5 component at a time:

```bash
ogprofiler hierarchy --run run/ --component-id 0
```

Component-level scheduling and parallel execution remain Phase 7 concerns.

## ComponentGraphLoader

`ComponentGraphLoader` addresses one Hive edge partition, reads only the
selected component's protein metadata, creates global-to-local and
local-to-global mappings, and constructs one root igraph inside the worker.
All descendants are induced from that root during an explicit DFS. Parent
subgraphs are not retained by queued siblings.

Species identity is represented by integer bitmaps. Cluster `n_species` is the
population count of the bitwise union; genome-name lists are not persisted.

## Split acceptance

Resolution search and acceptance are separate. Every tested candidate records:

- child count, minimum child size, and largest-child fraction;
- tiny-fragment fraction;
- Leiden quality;
- stability, adjusted Rand index, and normalized mutual information;
- inter-community and intra-community edge fractions;
- validity, selection, and rejection reason.

The lowest accepted resolution is selected. `min_split_quality` is optional;
when unset, quality remains recorded evidence rather than an uncalibrated hard
cutoff.

## Stochasticity modes

- `fast`: one deterministic seed;
- `robust`: three derived seeds and mean pairwise ARI/NMI;
- `publication`: configurable five-or-more seeds, with the most-supported
  observed partition selected as the consensus representative.

An accepted split must meet `stability_threshold`. The complete terminal-reason
set is `SINGLETON`, `NO_SPLIT`, `ONE_SPECIES`, `MIN_SIZE`, `UNSTABLE`,
`LOW_QUALITY`, `MAX_DEPTH`, `NO_EDGES`, and `GAMMA_LIMIT`.

## Component-local store

```text
hierarchy/components/component=00000000/
├── nodes.parquet
├── members.parquet
├── candidates.parquet
├── metrics.json
└── hierarchy-manifest.json
```

Nodes store bounded structural fields and statistics, including `child_count`.
Only terminal membership is persisted. The manifest binds the component edge
fragments, protein metadata, full resolved hierarchy configuration, and output
checksums; verified results are reused and damaged artifacts are reconstructed.
