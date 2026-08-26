# Phase 11 phylogenetic refinement smoke regression

## Provenance

- Run ID: `20260826T071500Z_phase11_phylogeny`
- Source family: reference E, `OG000000000`
- Family size: 40 genes across four species
- Base Git commit: `e46af3ed522d4d8b21870a23409791e4d08e1953`
- Deployment: captured dirty source state
- Alignment: MAFFT `v7.526 (2024/Apr/26)`
- Gene tree: FastTree `2.2.0 Double precision`
- Rooting: species-tree-aware minimum-duplication edge search
- Reconciliation: `lca-reconciliation-v1`

This was a short real-tool integration smoke on the completed Phase 9 reference
E outputs. Its runtime and memory measurements describe execution success and
scale only; they are not a biological accuracy benchmark.

## Result

| Metric | Value |
|---|---|
| Selected families | 1 |
| Family network event | `AMBIGUOUS` |
| Root phylogenetic event | `DUPLICATION` |
| Event confidence | 0.5 |
| Conflict status | `UNRESOLVED` |
| Supporting gene-tree node | `gene_node_000038` |
| Wall time | 2.56 s |
| Maximum RSS | 64,648 KiB |
| Immediate resume | 1/1 reused |

The supplied species tree contained one extra leaf. The archived pruned tree
contains exactly the four selected-family species. The original
`network_event=AMBIGUOUS` remains present beside `phylo_event=DUPLICATION`, so
the second evidence layer does not overwrite the network interpretation.

All family output hashes match `phylogeny-manifest.json`; the global event-table
hash matches `phylogenetic-manifest.json`. The archive includes input FASTA,
alignment, unrooted and rooted gene trees, pruned species tree, node-level
reconciliation, tool commands/versions, and deployment provenance.

## Validation

- local suite: 63 tests passed;
- final remote Conda suite: 63 tests passed;
- Ruff and mypy (59 source files) passed;
- midpoint, outgroup, and species-tree-aware rooting are covered;
- missing-species validation, pruning, selection, conflict status, corruption
  reconstruction, and verified resume are covered.
