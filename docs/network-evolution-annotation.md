# Phase 8 network evolution annotation

Network events are species-overlap heuristics on the network hierarchy. They
are not gene-tree reconciliation results and do not populate `phylo_event`.

```bash
ogprofiler annotate-network --run run/ \
  --set evolution.network_overlap_threshold=0.5
```

## Binary nodes

For child species bitmaps `S1` and `S2`, the continuous score is:

```text
|S1 & S2| / min(|S1|, |S2|)
```

Zero overlap is `SPECIATION_LIKE`. A positive score at or above the configured
threshold is `DUPLICATION_LIKE`; a positive score below it is `MIXED`. The
absolute overlap count, normalized score, and bounded confidence are retained.

## Polytomies

Every child pair is recorded as deterministic compact JSON containing child
cluster IDs, overlap count, and normalized overlap. An all-zero matrix is
`POLYTOMY`; uniformly positive high overlap is `DUPLICATION_LIKE`; heterogeneous
or partial patterns are `MIXED`. This preserves k-way Leiden structure rather
than manufacturing binary nodes.

Terminal one-species nodes are `SPECIES_SPECIFIC`. Other terminal nodes are
`AMBIGUOUS`. Consequently every stored hierarchy node has one event from the
closed schema:

```text
SPECIATION_LIKE, DUPLICATION_LIKE, MIXED,
POLYTOMY, AMBIGUOUS, SPECIES_SPECIFIC
```

## Legacy export mapping

Compatibility labels are derived independently of the internal event name:

- disjoint binary children → `I`;
- child intersection equals the parent species set → `II`;
- one child equals the parent species set → `III-1`;
- other partial binary overlap → `III-2`;
- k-way internal node → `III-3`.

The annotation algorithm itself does not branch on these labels.

## Artifacts and resume

Each hierarchy component produces:

```text
evolution/components/component=00000000/
├── events.parquet
└── network-event-manifest.json
```

The fixed event schema includes component/cluster identity, child count,
`network_event`, overlap count, overlap score, pairwise summary, confidence, and
legacy export label. Input and output hashes provide component-level verified
resume and reconstruct damaged event tables.
