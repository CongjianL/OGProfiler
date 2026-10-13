# Frozen endpoint allocation / complete DP paths — job 1412208

## Provenance and checks

Slurm COMPLETED /0:0; elapsed 00:03:06, 4 CPU /16G /2h.
RUN_ID `20261008T170938Z_665fc6874aac_722d001b_53277`.
Source `665fc6874aac4efd3e71a5fbec8cac17cde24d22`, dirty=0.
Source SHA256 `722d001b8daae5fd3c8ca2eaea1a4750b9abfbd73b0754e230ca2d59617e6471`.
Remote output `/home/mselab/licj/projects/running/ogprofiler-runs/20261008T170938Z_665fc6874aac_722d001b_53277/species-endpoint-paths`.

1205 fixed admitted sources; original graph, tree, legal nodes and exact objective.
All predictions replayed protein-by-protein against 1412190. Coverage, input hashes,
child-pair exact EII/EIX/EXX identities, direct+descendant DP margin identities and
actual subtree/LCA traversal loss checks passed during execution. Post-harvest,
strata/context_strata/degeneracy_strata/degeneracy_source_summary exactly match
1412190; subset totals and inactive losses/frontier counts independently checked.
Soft42 production unchanged. No extra run, score rule, qualification or extraction edit.

## Subset definition and actual unique LCA losses

Label-free subset: immediate complete-child internal-only gain <=0, full gain >0;
source denominators remain fixed. This differs from arbitrary per-block sign changes.
Reference labels enter only after features and exact cuts have been fixed.

| Source / local reference | Active split | Actual direct LCA TP / FP | Immediate frontier / deeper DP frontier | Inactive, covered by selected ancestor |
|---|---:|---:|---:|---:|
| clean / pure |44|744 /0|29 /15|81|
| polluted / pure |36|607 /0|32 /4|247|
| polluted / mixed |324|3189 /17100|166 /158|361|

404 active nodes, 689 inactive nodes (1093 total). Inactive direct actual losses
are all zero, irrespective of locally preferred split. No active_keep target exists:
a strictly positive direct gain already guarantees the local exact DP chooses split.

The active subset removes 4540 TP and 17100 FP, respectively 9.4882% of total
47849 TP losses and 2.7575% of total 620135 FP removals under this diagnostic cut.
Pure-node subset removes 1351 TP /0 FP. These are unique direct LCA attribution,
not a causal estimate of recovery after changing scores or forcing parent keeps.
43309 TP losses occur outside this particular sign-reversal subset. Thus the
observed mechanism explains some errors, but does not establish a dominant error
cause or support blanket external-term suppression.

227 active targets stop at their immediate children, 177 continue to deeper optimal
frontiers. In mixed targets, 158/324 use deeper frontiers; pure active targets 19/80.
Per-target subtree totals are nested and should not be summed across targets.
The summary deliberately aggregates only unique direct LCA losses, not subtree totals.

## 21396: exact endpoint allocation and direct-cut error

Fixed root path: 21394 split →21395 split →21396 split. No selected ancestor blocks
21396. Local and actual frontier both [21397,21444,21445]; each child locally keeps.
Composition: core 21397 has 24 proteins across 23 species (species18 has two copies);
21444 and21445 each contain one species11 protein. All 25 reference-assigned proteins
belong to one reference family. No deeper optimization gain contributes here.

The only decisive species block (1,11) has these endpoint strengths:

| Child | Internal species1 | Internal species11 | External species1 | External species11 |
|---|---:|---:|---:|---:|
|21397 core|1.8104164021539222|1.2743607263520085|0|0.12704059643643106|
|21444 singleton|0|0.5360556758019136|0|0|
|21445 singleton|0|0|0|0.1390442158172523|

External means node-external but source-internal, not source-external edges.
For core—21444, actual .5360556758019136 exceeds EII .46736499897321476:
deficit −.06869067682869885, EIX=EXX=0.
For core—21445, actual=EII=0, EIX=.12122695965386061, EXX=0.
No other unordered child pair contributes to the block.

Core's external .12704059643643106 gives within-child IX .11076149531342624,
which is excluded from between-child cut expectation. The singleton21445 external
.1390442158172523 creates cross-child EIX .12122695965386061. Accordingly, total
node-boundary weight .26608481225368336 is not one undifferentiated driver.
The other16 connected blocks cancel exactly and have no external expectation.

Source-normalized direct gain:
`−.00016032917480086512 + .0002829527863496145 = +.0001226236115487494`.
Descendant optimization gain exactly0; optimized split-minus-keep equals this direct
gain. Actual subtree loss and unique direct LCA loss both47 TP /0 FP. The error
arises at the legal 24+1+1 boundary, not from a deeper frontier, selected-ancestor
consumption or floating-point tie. This does not create a legal partial child union.

## Interpretation / next evidence

- Distinguish local score preference from actual consumption: all689 inactive targets
  illustrate why a local sign alone is not evidence of an extracted error.
- External endpoint position relative to children matters, not its total weight alone.
- The sign-reversal subset accounts for about9.5% of direct TP losses; most losses
  have a different mechanism under the current complete objective.
- Prioritize structurally matched pure/mixed nodes outside this subset where internal
  evidence already favors cutting, while retaining full legal merge/cut and downstream
  DP accounting. Test a clearly stated complete-objective hypothesis rather than an
  external-strength or degeneracy threshold. Current observations do not establish
  which new objective would improve a full partition.

Full endpoint-targets.json and cases.json remain remote; this archive preserves compact
summary and focal endpoint/path evidence. Harvesting this run did not submit another job.
