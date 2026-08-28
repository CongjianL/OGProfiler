# B1 comparator output mapping (frozen before formal scoring)

The target is the Open Orthobench root-level orthogroup: genes descended from one
ancestral gene at the root of the sampled species. Pairwise orthologue tables are not
terminal partitions and are excluded from primary scoring.

| Method | Frozen primary output | Semantic reason |
|---|---|---|
| OrthoFinder 3.1.5 | `Orthogroups/Orthogroups.tsv` | One row is a conventional OrthoFinder orthogroup across all input species. Pairwise orthologues and node-specific HOG tables answer different questions. |
| FastOMA 0.5.1 | `RootHOGs.tsv` | RootHOGs are the deepest/root-level gene families. `OrthologousGroups.tsv` contains strict marker-like groups with at most one representative per species; nested HOGs are taxonomic-level groups. |
| SonicParanoid2 2.0.9 | run-level multi-species ortholog-group table | SonicParanoid merges pairwise clusters into an MCL-derived multi-species OG partition. Pairwise species tables are excluded. The exact observed filename will be recorded after smoke validation without changing this semantic choice. |
| Proteinortho 6.3.6 | `<project>.proteinortho.tsv` | Standard protein-group output. PoFF/synteny output is excluded because the benchmark supplies protein FASTA only. |

`FASTOMA_PRIMARY_OUTPUT=RootHOGs.tsv`

## FastOMA species tree

Open Orthobench contains per-RefOG gene trees but no single bundled species tree whose
leaves correspond exactly to the 12 input proteome filenames. Stage 2A therefore freezes
`02_configs/fastoma/open_orthobench_ncbi_taxonomy_tree.nwk`, a topology derived from NCBI
taxonomy ranks and independent of RefOG labels or competitor scores. Branch lengths are
omitted, as supported by FastOMA. No candidate species trees will be scored and selected.
