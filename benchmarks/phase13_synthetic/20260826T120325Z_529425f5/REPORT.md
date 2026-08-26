# Phase 13 formal synthetic evolution benchmark

## Execution identity

- Slurm job: `1404203`, 16 strided array tasks, at most two concurrent
- Run ID: `20260826T120325Z_c1da58cee4ee_529425f5_54287`
- Source SHA-256: `529425f53e103f7c505a1f91bbe3e986ffa567180944a7f8b41635e7bb6e06e9`
- Git base: `c1da58cee4ee897eb3f30df0f4bf6c27bf171a8f`, exact dirty source snapshotted
- Resources: 4 CPUs, 16 GiB, 30 minutes per task
- Design: 243-scenario full factorial, one deterministic replicate

All 16 tasks completed with exit code `0:0` in 2:31–3:42. The archive contains
243 run summaries and 243 synthetic recovery documents. All 16 stdout and 16
stderr logs were audited, with zero traceback/error/failure/exception marker.
Individual scenario pipelines took 7.41–23.84 seconds (median 9.96 seconds).

## Applicability thresholds and counts

| Dimension | Threshold | Passing scenarios |
|---|---:|---:|
| Terminal-family F1 | ≥0.90 | 63/243 |
| Exact-clade hierarchy F1 | ≥0.80 | 4/243 |
| End-to-end event accuracy | ≥0.80 | 0/243 |
| Orthology F1 | ≥0.90 | 0/243 |
| All four dimensions | all thresholds | 0/243 |

This result separates capabilities rather than treating “applicable” as one
undifferentiated claim. OGProfiler 2 is robust as a terminal-family partitioner
in a restricted region, while its network hierarchy and heuristic event labels
do not reconstruct the complete simulated gene genealogy under this contract.

## Main scientific effects

Ancestral lineage expansion is the dominant terminal-family failure mode:

| Expansion | Mean family F1 | Family-applicable |
|---:|---:|---:|
| 1× | 0.944 | 61/81 |
| 2× | 0.698 | 1/81 |
| 4× | 0.690 | 1/81 |

Within the 1×-expansion region, mean family F1 is 0.977 at divergence 0.05,
0.963 at 0.20, and 0.892 at 0.40. Thus the evidence-backed terminal-family
applicability region is primarily non-expanded families at low-to-moderate
divergence; high divergence causes fragmentation, while expanded ancestral
lineages are usually merged or partitioned inconsistently.

Marginal fusion effects on family F1 are smaller (0.800, 0.768, 0.764 for
fusion probabilities 0, 0.10, 0.25) than the expansion effect. Duplication and
loss also have smaller marginal effects on family F1 than expansion. These are
descriptive full-factorial means from one deterministic replicate, not
confidence intervals.

Orthology F1 rises with divergence (mean 0.025, 0.114, 0.293) because higher
divergence induces more hierarchy splits and therefore more cross-child
candidates. It never reaches the declared 0.90 threshold. This is not evidence
that higher divergence improves true genealogy reconstruction: hierarchy F1
remains low and end-to-end event accuracy is near zero.

## Algorithm decision

The Phase 12 baseline remains the default for terminal-family inference; Phase
13 supplies no evidence for a parameter change. The scientific claim must be
narrowed as follows:

1. Terminal-family recovery is supported mainly for non-expanded ancestral
   families at low-to-moderate divergence.
2. The current Hierarchical Leiden topology is a network subdivision, not a
   gene-tree estimator.
3. `network_event` remains heuristic evidence and must not be reported as
   recovered duplication/speciation genealogy.
4. Hierarchy-derived pairwise candidates need phylogenetic refinement or a
   revised traversal/selection model before orthology-level claims.

The next scientific iteration should use multiple replicates and test a
hierarchy objective aligned to duplication-defined lineage clades. The current
243-scenario archive is the regression baseline for that work.
