# ADR 0001: Component-centric hierarchy architecture

- Status: Accepted for the Phase 2 prototype
- Date: 2026-08-25
- Owners: OGProfiler maintainers
- Supersedes: implicit hierarchy behavior in frozen OGProfiler 1

## Proposed amendment (2026-10-02)

ADR 0002 proposes separating singleton admission, recursion size stops and
unresolved search outcomes, with deterministic bounded fallback. Its status is
Proposed; this ADR's implemented strict child-size behavior remains in effect
until the replacement design is confirmed and implemented. See
`0002-hierarchy-admission-bounded-fallback.md` and the P6 acceptance audit.

## Context

Frozen OGProfiler 1 constructs a global igraph SSN, repeatedly creates induced
subgraphs during recursive hierarchy construction, creates multiprocessing
pools at hierarchy levels, and stores complete descendant gene lists on
internal hierarchy nodes. These choices couple memory use to global and
intermediate graph size and make hierarchy behavior difficult to reproduce and
resume.

Phase 2 tested a component-centric prototype against the frozen implementation
on identical SSNs. The prototype demonstrated lower hierarchy-core runtime and
RSS, but it also exposed a semantic difference: V1 creates nested singleton
terminal leaves for small connected components, while V2 can terminate a
component at its root. Isolated SSN vertices are omitted by the V1 hierarchy
membership export but retained by V2.

The hierarchy is a network community hierarchy. It is not a gene tree and is
not required to be binary.

## Decision

### Execution and object lifetime

1. The scheduling and persistence unit is one connected component.
2. A worker receives a `component_id`, loads only that component's vertices and
   edge partition, processes it root-to-terminal with an explicit DFS stack,
   writes component-local artifacts, and releases its graph.
3. A global igraph SSN is an import/prototyping representation, not a required
   production object.
4. The hierarchy engine does not create process pools. Component-level
   concurrency belongs to the scheduler.

### Hierarchy node model

Every node contains structural identity and bounded statistics:

```text
cluster_id
parent_id
component_id
depth
n_genes
n_species
resolution
quality
split_status
terminal_reason
```

Internal nodes do not persist complete descendant protein lists. Protein
membership is written exactly once as:

```text
protein_id -> terminal_cluster_id
```

Cluster IDs are deterministic only inside one component artifact. Cross-run
comparison uses the hash of the sorted descendant member set rather than a
transient numeric cluster ID.

### Root, terminal, and singleton semantics

1. Every connected component has exactly one root.
2. A root may also be terminal; a one-node hierarchy is valid.
3. Every imported SSN vertex, including an isolate, belongs to exactly one
   terminal family.
4. An isolate terminates with `SINGLETON`.
5. A small connected component may terminate at its root with `MIN_SIZE`; V2
   does not manufacture singleton leaves solely to reproduce V1 topology.
6. V1 omission of isolates and forced singleton leaves are recorded as
   classified compatibility differences, not adopted as V2 invariants.

### Splitting and resolution search

1. Leiden returns membership and quality; it does not return or retain child
   subgraphs.
2. Stable k-way splits are preserved. Binary topology is delegated to optional
   phylogenetic refinement.
3. Resolution search uses bounded exponential search followed by a local
   logarithmic grid.
4. The selected candidate is the lowest tested gamma satisfying all split
   acceptance rules:
   - more than one child;
   - each child has at least `min_child_size` members;
   - the largest child fraction is at most `max_child_fraction`.
5. Candidate traces record gamma, child count, quality, minimum child size,
   largest-child fraction, validity, method, seed, and algorithm version in the
   production hierarchy store.
6. Supported prototype methods are RBER, RBConfiguration, CPM, and Modularity.
   Method, seed, weight attribute, and resolved search configuration are run
   provenance.

### Terminal reasons

The Phase 2 closed set is:

```text
SINGLETON
MIN_SIZE
MAX_DEPTH
ONE_SPECIES
NO_EDGES
GAMMA_LIMIT
```

Production Phase 6 may add reasons only through a schema/version change. A
terminal node has one reason; a split node has none.

### Required invariants

For each component:

1. exactly one root exists;
2. every non-root node has exactly one parent in the same component;
3. the hierarchy is acyclic;
4. a split has at least two children;
5. children are pairwise disjoint;
6. the union of child descendants equals the parent descendants;
7. every component protein has exactly one terminal membership;
8. membership targets only terminal nodes;
9. stored `n_genes` equals the descendant membership cardinality.

Artifact publication occurs only after these invariants pass.

### Reproducibility and comparison contract

The run manifest records input hash, algorithm version, resolved configuration
hash, seed, Leiden method, weight attribute, and dependency versions. Equal
input, configuration, seed, algorithm, and dependency versions must produce
identical normalized membership and topology artifacts.

V1 compatibility is evaluated on an identical SSN using:

- common-ID terminal co-clustering precision, recall, and Jaccard;
- descendant-member-set node Jaccard;
- member-set parent-to-child edge Jaccard;
- hierarchy depth and terminal-family counts;
- wall time, hierarchy-core time, peak RSS, Leiden calls, and subgraph
  constructions.

Exact V1 topology is not a design requirement. Any topology change must be
classified, covered by a regression fixture, and reported.

## Consequences

### Positive

- Peak hierarchy memory is bounded primarily by the active component.
- Internal hierarchy storage is linear in nodes plus terminal membership.
- Component failures can later be isolated and resumed.
- k-way network structure remains distinct from phylogenetic bifurcation.
- Isolates have complete and explicit membership semantics.

### Costs and follow-up work

- V2 output intentionally differs from V1 for small components and isolates.
- Component-local Parquet serialization currently dominates command wall time
  on tiny fixtures and requires Phase 6/7 batching and scheduler work.
- Candidate traces are defined by this ADR but still need production storage.
- Species-aware terminal rules require species metadata integration.
- Phase 6 must replace prototype in-memory component extraction with
  partitioned edge loading.

## Phase 2 exit gate

Phase 2 closes when reference C/E runs establish all of the following:

- invariants pass for every component;
- repeated V2 runs have identical normalized membership and topology;
- the intended 240-protein E component is actually exercised;
- hierarchy RSS is lower than frozen V1 on C and E;
- Leiden calls and subgraph constructions are recorded and bounded;
- every observed V1/V2 structural difference is represented in the final
  regression report.
