# Internally positive complete-boundary comparison — job 1485398

## Provenance and validation

- Slurm COMPLETED /0:0, elapsed00:05:25; 4 CPU /16G /2h requested.
- RUN_ID `20261013T112225Z_323a8225e625_94e83f7f_88848`.
- Source `323a8225e625907c34181e080afebaf4cc790a5e`, dirty=0.
- Source SHA256 `94e83f7f5ac385da4f9825a04e8e20f7fd2ffc9f0d3a2385d76c76ec449e37b9`.
- Remote output `/home/mselab/licj/projects/running/ogprofiler-runs/20261013T112225Z_323a8225e625_94e83f7f_88848/species-positive-boundary`.

1205 fixed admitted sources, original induced graphs, legal nodes, source denominators
and exact DP. All cuts replayed protein-by-protein against1412190. Input hashes,
coverage, raw LCA flow/node score identity, complete-frontier boundary identity and
actual subtree/reference pair accounting verified during execution. Post-harvest,
strata/context_strata/degeneracy_strata/degeneracy_source_summary exactly match1412208.
All2445 actual split nodes partition exactly into three categories below and recover
47849 TP /620135 FP removals, including each source/local-reference stratum.
Mean direct+descendant gains match complete gains within serialized-float tolerance.

Labels enter post hoc; pure/mixed means relative to the frozen OF3 assignment, not
biological ground truth or an official QFO benchmark score. This diagnostic is not a
new full production partition. Production defaults and main checkout unchanged.
No additional job, scoring threshold, gamma, epsilon or eligibility rule was introduced.

## Unique direct LCA coverage of all actual split losses

| Category | Nodes | TP removed | TP share | FP removed | FP share |
|---|---:|---:|---:|---:|---:|
|Internal direct gain strictly positive|1292|34053|71.17%|580916|93.68%|
|Internal nonpositive, external context makes direct gain positive|404|4540|9.49%|17100|2.76%|
|Direct gain nonpositive, descendant optimum induces split|749|9256|19.34%|22119|3.57%|

The remaining90.5% outside the previous sign-reversal subset was not all internally
positive. The descendant-induced category is substantial. All totals above use
unique direct LCA losses, not sums of overlapping complete-subtree losses; they are
attribution under the existing cut, not causal recovery predictions for interventions.

Within pure reference nodes, internally positive cuts account for8823 of13623 TP
losses (64.77%):4776 in clean sources and4047 in polluted sources. Blanket suppression
of internally positive cuts would also target the mixed-node boundaries responsible
for580916 FP removals. These data locate a tradeoff; they do not establish the outcome
of a modified full objective and reoptimized DP.

## Internally positive active_split nodes

| Source / local reference | Nodes | Unique LCA TP / FP | Deeper frontier nodes | Mean direct gain | Mean complete gain | Mean descendant gain |
|---|---:|---:|---:|---:|---:|---:|
|clean / pure|114|4776 /0|34|.075751|.085134|.009383|
|polluted / pure|59|4047 /0|9|.028776|.031428|.002651|
|polluted / mixed|1119|25230 /580916|714|.065156|.109185|.044029|

43/173 pure targets (24.86%) and714/1119 mixed targets (63.81%) use a deeper optimal
frontier rather than their immediate children. Mixed-node mean complete gain includes
a larger descendant contribution, but these are unmatched, node-weighted descriptive
strata, not a controlled purity effect or an independent sample of biological families.

There are15 additional internally positive inactive nodes (4 pure,11 mixed, all in
polluted sources), with actual direct losses0. Their local preferences do not imply
actual extraction. Complete/subtree losses remain non-additive across nested targets.

## Strict common support: zero

Predeclared exact signature: protein count, species count, maximum copy count, sorted
immediate-child sizes, depth inside source, source protein count; additionally require
same original source cohort and actual DP state. Reference labels are not features.

Matched cells0, matched pure nodes0, matched mixed nodes0. The planned controlled
pure/mixed comparison therefore has no support in this dataset under this signature.
This is a limitation of the comparison, not evidence that purity is intrinsically
inseparable from graph structure. No matching criteria were relaxed after observing
this result. Mean scores above do not justify a classifier, gain threshold, or causal
claim. The full-boundary identity is validated; its ability to distinguish pure/mixed
families at fixed structure remains unestablished here.

## Interpretation and next decision

1. The main tradeoff lies in internally positive boundaries: they account for most
   current TP losses and most desired FP removals. A small-sign or external-only
   repair does not address this full-objective tradeoff.
2. Descendant optimization separately induces19.34% of direct TP loss. A proposed
   objective must be assessed through a complete legal DP cut, not immediate-child
   boundary scores alone.
3. Exact structural matching found no support. Before defining a new objective,
   explicitly delimit its comparison population and transfer assumptions; a new
   label-free complete merge/cut hypothesis needs a frozen full-partition comparison,
   not a threshold fitted to the present unmatched means. This report does not identify
   a validated replacement objective.

The full cases.json and positive-boundary-targets.json remain remote. Only compact
summary and comparison tables are archived. No further computation was submitted
while harvesting this run; production defaults remain unchanged.
