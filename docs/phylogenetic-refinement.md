# Phase 11 optional phylogenetic refinement

Phylogenetic refinement is a second evidence layer for selected terminal
network families. It does not replace the network hierarchy or rewrite
`network_event`.

```bash
ogprofiler annotate --run run/ --phylogenetic-refinement \
  --family OG000000123 \
  --species-tree species.nwk \
  --rooting species-tree-aware
```

## Selection policy

When one or more explicit `--family` values are supplied, only those families
are selected. Without explicit IDs, automatic selection is deterministic and
prioritizes `AMBIGUOUS`, then `MIXED`, large families, and families below a
`DUPLICATION_LIKE` ancestor. The configurable automatic cap prevents an
optional refinement request from expanding unexpectedly:

```text
phylogeny.selection_events
phylogeny.large_family_size
phylogeny.max_families
```

## Backend pipeline

The closed Phase 11 backend set is MAFFT plus FastTree. Executables come from
configuration, are discovered through `PATH` or an explicit path, and run
through one captured subprocess boundary. Commands, versions, stdout/stderr
provenance, input identity, and output checksums are retained per family.

```text
family FASTA -> MAFFT -> FastTree -> rooting -> LCA reconciliation
```

Supported rooting modes are:

- `midpoint`: branch-length midpoint of the gene-tree diameter;
- `outgroup`: root on the edge adjacent to the requested gene-tree leaf label;
- `species-tree-aware`: examine gene-tree edges and minimize the inferred
  duplication count against the pruned species tree.

Family gene-tree leaves use stable labels such as `P000000000123`. The same
label is used for `--outgroup`.

## Species tree and reconciliation

The built-in Newick parser supports branch lengths, internal labels, quoted
leaf names, validation, deterministic writing, and pruning. Species-tree
leaves map to prepared `species_name` values. Species absent from a selected
family are pruned; a selected species missing from the supplied tree is an
error.

LCA reconciliation labels gene-tree internal nodes as `SPECIATION` or
`DUPLICATION`, records confidence and the mapped species-tree node, and reports
the rooted gene-tree event as the family summary.

Outputs are stored under:

```text
evolution/phylogenetic/
├── phylogenetic-events.tsv
├── phylogenetic-manifest.json
└── family=OG000000123/
    ├── input.faa
    ├── alignment.faa
    ├── gene_tree.unrooted.nwk
    ├── gene_tree.rooted.nwk
    ├── species_tree.pruned.nwk
    ├── reconciliation.tsv
    └── phylogeny-manifest.json
```

The summary contains both `network_event` and `phylo_event`, plus
`event_confidence`, `supporting_node`, and `conflict_status`. Ambiguous network
evidence remains visible as `UNRESOLVED`; comparable disagreement is retained
as `CONFLICT`.
