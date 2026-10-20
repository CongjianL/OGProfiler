# Same-attraction / cross-null complete DP intervention — job1553883

## Provenance and validation

Slurm COMPLETED /0:0, elapsed00:03:24; 4CPU /16G /2h requested.
RUN_ID `20261020T061804Z_302ff9b86602_56766154_15420`.
Source `302ff9b86602029738769a4af0b72ba19274bddf`, dirty=0.
Source SHA256 `56766154f2f8c0ec45af295536c4182c2ff6720599d3e4643921fff2cfd620bd`.
Remote output `/home/mselab/licj/projects/running/ogprofiler-runs/20261020T061804Z_302ff9b86602_56766154_15420/same-attraction-intervention`.

Algorithm `same_attraction_cross_null_v1`: S*(v)=(A_same+A_cross−E_cross)/W.
Original graph, tree, all legal nodes, source denominators and exact parent-keep ties.
1205 frozen source cuts replay1553882 protein-by-protein; all candidate cuts optimized
before reference labels are loaded. Independent raw-edge/root objective, feasible-cut
ordering, complete coverage, coarsening, unique-LCA pair transitions and global-hybrid
joint-partition identities pass during execution. Original input/production members
hashes recorded and checked. Post-harvest all four existing strata summaries equal
1553882; metric formulas, fixed universes and cohort totals independently verified.
No production/default or main checkout changes. No second candidate or parameter search.

## Explicit all-reference hybrid universe

Reference assigned62466 proteins,385205 reference pairs; total production proteins82507.
Two comparison arms replace the same1205 source clades (30399 total proteins) with their
respective original or candidate complete cuts. Outside proteins keep frozen production
families, but source and outside namespaces are identically separated in both arms.
This common replacement can itself split original production groups spanning sources
or source/outside. Therefore this is a full hybrid diagnostic comparison, not an
end-to-end production run. The frozen production comparator below is separate.

| Arm | TP | FP | FN | Precision % | Recall % | Pair F1 % |
|---|---:|---:|---:|---:|---:|---:|
|Original exact DP hybrid|239907|144118|145298|62.471714|62.280344|62.375882|
|Candidate DP hybrid|255992|179573|129213|58.772399|66.456043|62.378498|
|Frozen soft42 production, separate comparator|209075|64005|176130|76.561813|54.276294|63.521119|

Candidate versus original hybrid: precision −3.699315 percentage points, recall
+4.175699pp, F1 +0.002616pp. This tiny numerical F1 increase is essentially a tradeoff,
not a demonstrated robust improvement. No uncertainty/held-out replication is claimed.
Candidate-versus-production differences include the shared source replacement and
namespace operation; they are not attributable solely to the null ablation.

## Target reference universe

29372 assigned target proteins,170196 reference pairs (including reference pairs
spanning different sources, which remain unrecovered in both arms).

| Arm | TP | FP | Precision % | Recall % | Pair F1 % |
|---|---:|---:|---:|---:|---:|
|Original exact DP|117958|118128|49.963996|69.307152|58.067057|
|Candidate|134043|153583|46.603228|78.758020|58.556819|

Target F1 +0.489762pp has a different denominator and is not the full hybrid F1 effect.
Both reference systems here use frozen OF3 assignments, not biological ground truth
or official QFO evaluation scores.

## Unique pair transitions and ancestor propagation

- Recovered TP16085, new merge FP35455; approximately2.204 new FP per recovered TP.
- New split TP loss0, removed FP0: the exact coarsening invariant is observed.
- Changed sources234/1205: clean52/440 recover2783 TP /0 added FP; polluted182/765
  recover13302 TP with35455 added FP.
- Local decisions change918 nodes; actual active states change1304 nodes;700 state
  changes have unchanged local decisions and arise from ancestor coverage changes.
  These counts overlap: they are distinct diagnostics, not disjoint quantities to sum.

Unique LCA TP/FP transitions equal full partition changes; no nested subtree loss sum
is used. These are observed intervention effects in the specified frozen hybrid, not
causal claims about biological families or deployment behavior on other datasets.

## Focal boundaries

21396 remains active_split with [21397,21444,21445]. Direct and optimized margin both
remain .0001226236115487494; removed same-null boundary gain0. Its47 TP /0FP loss
therefore remains. Full source21394 still has5 groups, TP278/FP0/lost352.
This boundary was cross-null driven, so removing only same null leaves it unchanged.

10052 remains active_keep, exact margin0, removed same-null gain0. Source10040 retains
6 groups, TP102/FP0/lost333. The earlier exact tie repair remains intact.

## Decision

The experiment establishes the predicted direction: suppressing same null reduces
some overcutting, but adds substantially more merge FP and yields nearly unchanged
all-reference hybrid F1. It does not support promoting this candidate to production.
Soft42 defaults remain unchanged; no follow-up computation is submitted on harvest.

If pursuing the mechanism, freeze these two cuts and attribute new FP to the responsible
merge nodes and ancestor-consumption paths, alongside the TP recovered at those same
nodes. Separate direct same-null removal from descendants it consumes. This is a
posthoc diagnostic of a fixed intervention, not a proposal to sweep coefficients or
fit a purity/gain filter. Historical strict matching common support0 remains a limit
on controlled pure/mixed score comparisons.

Full cases and complete predictions stay remote. This archive preserves compact
summary and metric/cohort comparison evidence; production results are not overwritten.
