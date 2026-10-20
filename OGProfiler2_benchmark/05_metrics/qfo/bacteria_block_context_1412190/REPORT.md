# Fixed legal-node source-context decomposition — job 1412190

## Provenance and validation

- Slurm COMPLETED, exit 0:0, elapsed 00:02:40; 4 CPU /16G /2h.
- RUN_ID `20261008T144003Z_dbb1ae780742_19746445_22299`.
- Source `dbb1ae7807426de25127ae08e976d0774aa1fd66`, dirty=0.
- Source SHA256 `197464459b1db67e3e48307e0a3b8cb27681721483f96f3935e9bc892f20e2a0`.
- Remote output `/home/mselab/licj/projects/running/ogprofiler-runs/20261008T144003Z_dbb1ae780742_19746445_22299/species-block-context`.
- Frozen 1205 complete admitted sources, 24462 legal internal nodes; same induced graphs,
  source denominators, complete-child choices and repaired exact DP as job 1412157.
- Exact per-block and complete immediate-child identities asserted during execution.
  Coverage and frozen input hashes verified. Local post-harvest check confirms all
  `strata` are identical to 1412157; context aggregates sum within float-output tolerance.
- Reference labels used only post hoc; production soft42 unchanged. No new partition,
  threshold, epsilon, qualification rule, gamma/resolution search or source-external edges.

## Definition

For each source block and legal immediate-child cut, split endpoint strengths into
I (edges with both endpoints inside the node) and X (edges crossing the node boundary,
but with both endpoints in the source). Keeping the original source block denominator,
`E = EII + EIX + EXX`; EIX includes both orientations. Actual is the original
between-child weight. Thus `gain = (EII-actual)/Wsource + (EIX+EXX)/Wsource`.
This is an explanatory decomposition, not a local null refit or alternate DP objective.
Edges wholly outside the node affect the fixed denominator but contribute no node endpoints.

## Actual active_split nodes

| Source / local reference cohort | Legal nodes | Internal nonpositive → full positive | Sum internal gain | Sum external gain | Unique LCA removed TP / FP |
|---|---:|---:|---:|---:|---:|
| clean / pure |279|44 (15.77%)|6.817752|2.080236|7561 /0|
| polluted / pure |137|36 (26.28%)|0.965724|0.866087|6062 /0|
| polluted / mixed |2029|324 (15.97%)|64.808068|8.941157|34226 /620135|

Gain sums normalize each node by its own source total. Nested-node sums are diagnostic
instances, not an independent sample, pooled graph score or unique graph-edge fraction.
The TP/FP column belongs to all nodes in each row, not just the sign-changing subset.
Sign changes refer to direct complete-child cuts; child-optimal DP paths may differ.

## Cancelling versus decisive blocks

For clean/pure active_split nodes:
- 1684 connected-equal blocks; 325 have external endpoint strength. Summed internal
  gain −0.910091 and external gain +0.910091 cancel exactly per block.
- 349 positive decisive blocks; 88 internally nonpositive become positive due to context.
  Internal gain +8.464546, external gain +1.081115.
- 143 negative decisive blocks; internal gain −0.736702, external gain +0.089029.

For polluted/pure active_split nodes:
- 380 connected-equal blocks, 81 with external endpoints: −0.143807 +0.143807.
- 277 positive decisive blocks, 91 context-induced positives: +1.455606 +0.672220.
- 248 negative decisive blocks: −0.346075 +0.050060.

For polluted/mixed active_split nodes:
- 19850 connected-equal blocks, 3255 with external endpoints: −2.893091 +2.893091.
- 24065 positive decisive blocks, 2222 context-induced positives: +72.379447 +4.897624.
- 6515 negative decisive blocks: −4.678288 +1.150442.

These context block counts include same-species blocks; the previous one-sided
cross-species degeneracy counts have a different definition and denominator.
External terms are nonnegative algebraically and occur in both pure and mixed strata;
blanket removal is not supported by this observational audit. Removing them alone would
also change currently cancelling blocks into negative evidence for cutting.

## Key nodes

### 21396 — real positive boundary, context sign reversal

16 connected-equal blocks contribute zero, all with EIX=EXX=0.
The sole decisive block (1,11):
- actual 0.5360556758019136;
- EII 0.46736499897321476;
- EIX 0.12122695965386061; EXX 0;
- internal-only deficit −0.06869067682869885;
- node-internal block weight 1.8104164021539222;
- node-boundary weight 0.26608481225368336;
- wholly node-external source-internal weight 0.

Normalized direct complete-cut gain:
`−0.00016032917480086512 + 0.0002829527863496145 = +0.0001226236115487494`.
The exact DP actively splits [21397,21444,21445], removing 47 reference TP /0 FP.
This is not numerical tie noise. Previous job 1412154 traced the singleton 31901
boundary edge of weight 0.1390442158172523; its weight is a subset of this block's
0.26608481225368336 boundary total. Only cross-child endpoint products enter EIX;
not every boundary edge contributes identically. The present aggregate does not
replace that earlier endpoint-level evidence.

### 10052 — genuine cancelling tie

One block (15,16), actual=EII=0.9655215666683362, EIX=EXX=0;
complete direct gain exactly zero, active_keep. Counterfactual direct split would
remove 9 TP, but this is not an actual loss.

## Next decision

Prioritize endpoint allocation and complete DP-path influence in the sign-changing
subset, with pure and mixed comparison cases matched on source/species/copy/tree
structure. Attribute unique LCA TP/FP losses to that subset explicitly. Do not infer
its losses from the full-stratum totals above. Separate within-child boundary
products (which cancel) from between-child products driving EIX, and retain complete
legal child partitions when tracing downstream DP choices. These results support
further diagnosis, not a degeneracy-ratio filter or automatic production score change.
No additional job was submitted when harvesting this run.
