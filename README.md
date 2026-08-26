# OGProfiler 2

OGProfiler 2 infers hierarchical protein families from sequence-similarity
networks, with optional phylogenetic refinement. The V2 implementation is
developed under `src/ogprofiler`; the frozen V1 reference is under `legacy/`.

## Development preview

```bash
python -m pip install -e '.[dev]'
ogprofiler prepare --proteomes testdata/ --out run/
ogprofiler search --run run/ --backend diamond
ogprofiler edges --run run/ --method lrb
ogprofiler components --run run/
ogprofiler hierarchy --run run/ --component-id 0
ogprofiler hierarchy-all --run run/ --set runtime.workers=4
ogprofiler annotate-network --run run/
ogprofiler export --run run/
ogprofiler orthologs --run run/ --emit-pairwise-orthologs
ogprofiler annotate --run run/ --phylogenetic-refinement --family OG000000123
ogprofiler benchmark plan --out benchmark-plan/
ogprofiler benchmark metrics --run run/ --ground-truth truth.tsv \
  --dataset reference-a --method ogprofiler2 --out scientific-metrics.json
```

The prepare stage assigns deterministic integer species/protein IDs, writes
normalized FASTA and Parquet metadata, and captures the resolved configuration
and input checksums in the run workspace.

The search stage builds one global DIAMOND database, runs directional
all-vs-all search, writes `search/hits.parquet`, and publishes a checksummed
search manifest with verified resume semantics.

The edge stage reproduces legacy NBS normalization, applies explicit coverage
thresholds, performs RBH/LRB retention, and writes order-invariant canonical
edges without species-pair sparse matrices.

The export stage assigns deterministic dataset-scoped OG IDs from canonical
terminal memberships and writes `families.tsv`, `members.tsv`, `hierarchy.tsv`,
and `events.tsv`. Per-family FASTA and component GraphML are explicit opt-in
exports.

Pairwise ortholog candidates are disabled by default because output can be
quadratic. When enabled, the hierarchy-driven orthology stage streams
cross-child, cross-species candidates directly to `ortholog_pairs.tsv.zst`.

Selected terminal families can receive a second, explicitly separate evidence
layer through MAFFT, FastTree, configurable rooting, and LCA reconciliation.
The resulting `phylo_event` is stored beside, rather than over, the original
`network_event`.

Scientific algorithm selection uses a checksummed parameter matrix and common
family, evolutionary-consistency, and orthology metrics. Frozen V1 membership
and other methods can be compared through the same versioned result schema.
