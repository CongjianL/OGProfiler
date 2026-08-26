# Phase 13 synthetic evolution benchmark

Phase 13 provides deterministic protein-evolution fixtures with explicit family,
genealogy, event, orthology, and domain-architecture truth. Its purpose is to
measure where the network hierarchy remains informative, rather than merely
confirming that a pipeline invocation exits successfully.

## Scenario design

The default matrix is a 243-row full factorial:

| Axis | Levels |
|---|---|
| Sequence divergence | 0.05, 0.20, 0.40 |
| Terminal duplication probability | 0.00, 0.25, 0.50 |
| Gene-loss probability | 0.00, 0.20, 0.40 |
| Ancestral lineage expansion | 1, 2, 4 |
| Fusion probability | 0.00, 0.10, 0.25 |

Replicates are an explicit additional dimension. Horizontal transfer remains a
future axis and is not silently approximated by fusion.

```bash
ogprofiler benchmark synthetic-plan --replicates 1 --out phase13-plan/
ogprofiler benchmark synthetic-generate --index 0 --out scenario-000/
```

The generator uses a scenario-derived random seed, a fixed four-species tree,
and amino-acid substitutions from a closed alphabet. Each ancestral family may
expand into independently inherited lineages. Loss is recorded as an explicit
genealogy node; terminal duplication creates co-ortholog copies; fusion appends
a donor half-domain while retaining an explicit primary and secondary family.

## Ground truth contract

Each scenario contains:

```text
proteomes/*.faa
ground_truth.tsv
genealogy.tsv
true_events.tsv
true_orthologs.tsv
domain_architecture.tsv
species_tree.nwk
manifest.json
```

`manifest.json` records the complete scenario, cardinalities, and artifact
checksums. `ground_truth.tsv` uses the Phase 12 family-metric schema. Genealogy
nodes have stable parent links and the closed event set `ANCESTRAL_FAMILY`,
`LINEAGE`, `SPECIATION`, `DUPLICATION`, `LOSS`, and `GENE`.

A fused sequence receives one primary family for partition metrics and a
secondary domain family in `domain_architecture.tsv`. This deliberately exposes
bridge-induced family merging without introducing ambiguous multi-membership
into the terminal-family partition.

## Recovery metrics

```bash
ogprofiler benchmark synthetic-metrics \
  --run RUN --dataset-root SCENARIO --out synthetic-metrics.json
```

The `synthetic-recovery-v1` result contains:

- terminal-family pairwise precision, recall, F1, and exact recovery;
- hierarchy exact-clade precision, recall, F1, and mean best-truth Jaccard;
- event exact-clade coverage, confusion counts, matched accuracy, and
  end-to-end accuracy;
- streamed orthology precision, recall, and F1 against `true_orthologs.tsv`.

Event accuracy requires both recovery of the true descendant set and the
correct event label. Unrecovered truth events therefore reduce end-to-end
accuracy instead of disappearing from the denominator.

## Applicability map

The default applicability rule requires all of:

```text
family F1       >= 0.90
hierarchy F1    >= 0.80
event accuracy  >= 0.80
orthology F1    >= 0.90
```

Thresholds and separate family, hierarchy, event, orthology, and overall
pass/fail decisions are serialized in `leiden-applicability-v1`; they are not
inferred after viewing results. This separation prevents strong terminal-family
recovery from being hidden by weak genealogy reconstruction, or vice versa.

```bash
python benchmarks/aggregate_phase13_applicability.py \
  --root RUN_ROOT/results --out RUN_ROOT/results/summary
```

The formal Slurm script uses 16 array tasks with at most two concurrent. Each
task executes a deterministic stride through the 243 scenarios, avoiding 243
separate scheduler allocations while preserving one result directory per
scenario.

## Tracer interpretation

The first zero-duplication, zero-loss, zero-fusion scenario is a pipeline tracer,
not an applicability conclusion. It verifies generation, DIAMOND search, edge
construction, hierarchy, event export, streaming orthology, and metric parsing.
The full factorial is required before drawing a boundary over divergence,
duplication, loss, expansion, or fusion.

## Formal result

Slurm job `1404203` completed all 243 scenarios. Family F1 passed the declared
0.90 threshold in 63 scenarios, exact-clade hierarchy F1 passed 0.80 in four,
and no scenario passed the event or orthology threshold. Expansion was the
dominant family-recovery failure: mean family F1 was 0.944 for one lineage and
approximately 0.69 for two or four lineages. The supported claim is therefore
terminal-family recovery in a restricted non-expanded, low-to-moderate
divergence region. Network topology, event annotation, and cross-child
orthology remain distinct heuristic outputs rather than recovered gene
genealogy. The complete evidence is archived under
`benchmarks/phase13_synthetic/20260826T120325Z_529425f5/`.
