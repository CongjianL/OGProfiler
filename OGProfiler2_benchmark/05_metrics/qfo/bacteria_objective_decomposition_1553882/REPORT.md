# Complete-objective channel decomposition — job1553882

## Provenance and checks

Slurm COMPLETED /0:0, elapsed00:04:33; 4CPU /16G /2h requested.
RUN_ID `20261020T004942Z_617e3f7508b8_64c81293_91868`.
Source `617e3f7508b8e2bea4c104c58e00ac63f140b0c7`, dirty=0.
Source SHA256 `64c812930c09721bd58df4e67bdc40c89bf138485aec0e2370d660719c202786`.
Remote output `/home/mselab/licj/projects/running/ogprofiler-runs/20261020T004942Z_617e3f7508b8_64c81293_91868/species-objective-decomposition`.

1205 fixed sources,24462 legal internal nodes. Exact predictions replayed protein-by-protein
against1485398. Stored raw LCA edge flows, exact keep scores, same/cross-species null,
direct+descendant=forced-split, optimal-frontier and actual-LCA/root telescope identities
verified during execution. Coverage and immutable hashes pass. Post-harvest four existing
strata summaries equal1485398 exactly. Serialized gain decompositions, sign counts and
source-channel versus active direct-LCA sums independently checked within float tolerance.
2445 actual split nodes recover47849TP /620135FP removals; inactive/keep direct actual
losses all0. Reference labels are posthoc; no matching or objective modification occurred.

## Interpretation of units

Below values are sums of node contributions individually divided by their source graph
weight W. Raw pooled expected−actual can have a different sign due to heterogeneous W;
compare normalized expected/actual with normalized gains, not raw pooled differences.
Only active direct-LCA boundary terms telescope without double counting. Descendant,
forced and optimal summaries below overlap across nested targets; they describe where
local decisions obtain gain, not independent graph evidence or unique reference losses.
Same/cross-species graph channels are not reference TP/FP channels.

## Active internally positive boundaries

| Source/local reference | Nodes | Direct same gain | Direct cross gain | Descendant same | Descendant cross | Optimal same | Optimal cross |
|---|---:|---:|---:|---:|---:|---:|---:|
|clean/pure|114|4.445844|4.189785|.343484|.726183|4.789328|4.915968|
|polluted/pure|59|1.124936|.572868|.032976|.123443|1.157912|.696312|
|polluted/mixed|1119|10.004028|62.905336|8.397714|40.870777|18.401742|103.776113|

In the173 pure nodes together, direct same gain5.570780 versus cross4.762654:
same contributes53.91% of summed direct gain, and51.45% of summed local optimal gain.
In1119 mixed nodes, same contributes13.72% direct and15.06% local optimal gain.
These are unmatched observational sums, not per-node gain fractions or controlled effects.
Strict common support was0 in1485398, and is not reconstructed or relaxed here.

The pure same-species direct gain comes from normalized null expectation9.688112 minus
actual4.117332; cross gain comes from13.011172 minus8.248518. Mixed corresponding
same expectation13.896805 minus actual3.892777; cross107.458013 minus44.552676.
Thus both channels favor cutting in pure-node aggregate. For pure targets, cross direct
gain is positive at106/173 nodes, negative9, zero58; removing same null does not address
all pure boundaries. These numbers do not assign their8823 TP losses to either channel.

## Direct nonpositive / descendant-driven active boundaries

| Source/local reference | Nodes | Direct same | Direct cross | Descendant same | Descendant cross |
|---|---:|---:|---:|---:|---:|
|clean/pure|121|−.145072|−.020072|1.609412|2.661345|
|polluted/pure|42|−.020656|−.048599|.671049|.459994|
|polluted/mixed|586|−.438330|−.318891|4.353235|17.303808|

Descendant gains overcome direct costs in both pure and mixed cohorts. A direct-boundary
threshold would miss propagation of optimal cuts through ancestors. A new objective
requires full reoptimization, including keep versus forced-split, not just relabeling
these local margins. Original unique direct-LCA loss coverage remains71.17% internally
positive,9.49% context reversal,19.34% descendant-driven; this diagnostic changes none.

## One proposed experiment, not a validated correction

Candidate hypothesis `same_attraction_cross_null_v1`: keep all observed edges as attraction,
but retain configuration-null expectation only for cross-species blocks:

`S*(v) = [A_same(v) + A_cross(v) − E_cross(v)] / W`.
`Q*(P) = sum_{v in legal complete cut P} S*(v)`.

Original `S(v)=[A_same+A_cross−E_same−E_cross]/W`; equivalently
`S*(v)=S(v)+E_same(v)/W`. Source denominators and all graph/tree/eligible nodes unchanged.
No coefficient, resolution search, epsilon, label-based weighting or boundary veto.
Scores remain node-additive, so the exact same complete-tree DP optimizes Q*.
For any fixed complete boundary, candidate gain is
`[E_cross_boundary − A_cross_boundary − A_same_boundary]/W`, i.e. original gain minus
`E_same_boundary/W`. Same-species edges remain attractive; removing all same-species
edges or dropping their observed weight is a different experiment and is not proposed.

Rationale: pure-node aggregate shows a relatively larger positive same-null contribution,
while mixed-node aggregate is cross-dominated. This only motivates a falsifiable ablation;
it neither proves biological misspecification nor predicts improvement. Mixed boundaries
also receive same-null gain, and real paralog separations may be suppressed. Pure cross
gain remains. The zero-common-support limitation forbids a controlled purity claim.

Root score is generally no longer0: `S*(root)=W_same/W`. Report partition gain relative
to that single fixed root baseline; do not subtract S*(root) from every selected-node
score (that introduces a hidden group-count penalty). W=0 gives all0/tie keep parent.
Cross-only graphs reproduce the original exact objective and cut; same-only graphs
prefer complete merge (ties keep parent). For any fixed cut its gain versus merge
cannot increase; after global reoptimization boundary choices need not simply coarsen.

Before running: implement and test the additive formula, independent cut objective,
nonzero root baseline, exact ties and those limiting cases; then compare one frozen
candidate against original complete DP. Report TP/FP/P/R/F1 over an explicitly fixed
reference universe, pair transitions and unique-LCA losses, local decision and ancestor
path changes. Nested subtree totals must not be added. No performance or recovery
estimate is established in this report; no candidate was implemented or submitted here.

Full cases.json stays remote. Compact summary and active channel strata are archived.
Production soft42 and current main checkout remain unchanged.
