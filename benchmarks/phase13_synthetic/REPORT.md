# Phase 13 synthetic evolution tracer report

## Implemented benchmark contract

- Deterministic 243-scenario full factorial over divergence, duplication, loss,
  ancestral expansion, and fusion.
- Explicit family membership, rooted genealogy, duplication/speciation/loss
  events, cross-species ortholog pairs, domain architecture, and species tree.
- Terminal-family, exact-clade hierarchy, end-to-end event, and streamed
  orthology recovery metrics.
- Separate family/hierarchy/event/orthology applicability decisions plus a
  strict overall decision.
- Sixteen-task strided Slurm runner, capped at two concurrent allocations.

## Remote end-to-end tracers

Both tracers ran in the configured remote `ogprofiler` Conda environment using
real DIAMOND search and the complete OGProfiler 2 pipeline.

| Scenario | Proteins | Family F1 | Hierarchy F1 | Event accuracy | Orthology F1 | Pipeline seconds |
|---|---:|---:|---:|---:|---:|---:|
| 0.05 divergence; no duplication/loss/expansion/fusion | 48 | 1.000 | 0.500 | 0.000 | 0.000 | 4.054 |
| 0.40 divergence; 0.50 duplication; 0.40 loss; 4× expansion; 0.25 fusion | 189 | 0.705 | 0.498 | 0.006 | 0.515 | 5.106 |

The baseline tracer perfectly recovers terminal families but does not reproduce
the within-family species genealogy. Its terminal components therefore produce
no cross-child ortholog candidates. This is a scientific result of the declared
metric, not an execution defect: a network family partition and a gene tree are
different structures.

The stress tracer completes successfully and shows the expected degradation in
family and orthology recovery under simultaneous high divergence, duplication,
loss, expansion, and fusion. It confirms that the matrix contains discriminating
conditions, while two endpoints alone do not define an applicability boundary.

## Next execution boundary

The full 243-scenario matrix completed as Slurm job `1404203`. The formal report,
243 compact metric documents, axis summary, scheduler accounting, and provenance
are archived under `20260826T120325Z_529425f5/`.
