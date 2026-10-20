# Phase 10 hierarchy-driven orthology engine

Pairwise ortholog output is an explicit, optional projection of the rooted
network hierarchy:

```bash
ogprofiler orthologs --run run/ --emit-pairwise-orthologs
```

The default configuration sets `output.emit_pairwise_orthologs=false`, so the
potentially quadratic artifact is absent unless requested.

## Candidate semantics

For every `SPECIATION_LIKE` node, the engine collects genes below each direct
child and emits the Cartesian product between distinct children. Phase 8
defines `POLYTOMY` as a k-way node whose child species-overlap matrix is all
zero, so those nodes use the same cross-child rule without manufacturing a
binary tree.

Pairs from the same species are discarded. Multiple retained genes from one
species may pair with genes from another child and are labeled
`CO_ORTHOLOG_CANDIDATE`. Same-child genes are not emitted at that node;
in-paralog status is not inferred by this stage. Each pair is supported by the
unique lowest hierarchy node whose direct child branches separate the genes.

The output columns are:

```text
protein_a_id
protein_b_id
species_a_id
species_b_id
component_id
supporting_cluster_id
relationship
```

Protein IDs are canonicalized so `protein_a_id < protein_b_id`.

## Streaming and resume

The generator processes one component at a time and holds descendant gene
lists, not generated pairs. A bounded byte buffer is flushed every
`output.ortholog_pair_chunk_size` candidates through Arrow's Zstandard stream
to:

```text
results/ortholog_pairs.tsv.zst
```

Peak Python memory is therefore bounded by the component hierarchy, terminal
membership, and configured chunk size rather than total pair count. Progress
is reported after each output chunk.

`results/ortholog-manifest.json` records upstream hierarchy, membership,
event, and protein checksums; algorithm parameters; output checksum; pair,
chunk, component, and supporting-event counts. Matching runs reuse the verified
compressed output, while corruption or changed inputs trigger atomic rebuild.


## P4 grouping dependency contract

The policy version is `hierarchy-cross-child-v2`. Pairwise candidates continue
to consume the terminal hierarchy and `network_event`; the manifest records
`grouping_dependency=none` and the event algorithm version. OG memberships,
OG coverage, and `v1_event` are not pair-generation inputs. In particular,
cross-species membership within one OG does not imply a pairwise ortholog.
The standard `run` pipeline appends this stage only with the explicit
`output.emit_pairwise_orthologs=true` opt-in.
