# Phase 4 numeric edge engine

The edge stage converts the standard directional search table into normalized,
filtered, canonical undirected edges without BH/RBH/CM sparse matrices:

```bash
ogprofiler edges --run run/ --method lrb
```

It writes:

```text
edges/normalized_hits.parquet
edges/retained_edges.parquet
edges/edge-manifest.json
```

## Legacy NBS (OF3.1.5-aligned)

The configuration name is retained, but the default model now follows OF3.1.5
full-run scoring rather than frozen V1 bins:

1. Remove self hits and nonpositive scores; max-bitscore deduplicate each direction.
2. Traverse hits in sparse query/target ID order, then stable-sort by length product.
3. Below 100 hits, use all samples. Otherwise use disjoint full bins of width
   20/200/1000, discard the incomplete tail, and keep scores at/above NumPy's 95th percentile.
4. Fit `log10(bitscore) = a * log10(length_product) + b` with SciPy curve_fit.
5. Normalize with the corresponding diagonal length scaling.

Equal length products also use curve_fit, not max-score normalization. Solver
failures surface rather than silently changing the model. With fewer than two
selected fit points, `similarity.nbs_fallback=v1_zero` (default) drops the group;
explicit `v2_max` uses `a=0,b=log10(max_bitscore)` instead.

## Coverage

Coverage filtering occurs after NBS so changing coverage thresholds does not
change the fitted normalization model. It is **disabled by default** to match
V1 and OrthoFinder3, which never filter by coverage:

```yaml
edges:
  apply_coverage_filter: false
  min_query_coverage: 0
  min_target_coverage: 0
  min_bidirectional_coverage: 0
```

When `apply_coverage_filter` is enabled, a direction is retained only if it
satisfies all three configured limits. `min_bidirectional_coverage` applies to
`min(query_coverage, target_coverage)` for each directional hit.

## Best hits, RBH, and LRB

Best hits are selected per query and target species using the absolute
`best_hit_tolerance` (default `1e-3`) and V1's strict comparison
`score > maximum - tolerance`.

Same-species paralogs are best hits when their score reaches the query's best
external-species score within the same tolerance. Reciprocal best hits require
both directed best-hit keys and are restricted to different species.

LRB reproduces OF's cutoff, including repeated-query NumPy assignment: within
each target-species RBH matrix the last target in CSR order supplies the score,
then minima are taken across target species. This differs from taking the true
minimum of all near-tied RBHs. Hits at/above that cutoff are retained, including
same-species paralogs. The no-RBH fallback is `best_external_score + 1e-6`;
a query with only same-species hits uses `1e-6`. OF's sentinel bounds are preserved.

Compatibility modes are explicit:

- `rbh` and `arb`: reciprocal best hits;
- `ar`: any hit for which the reverse direction exists;
- `lrb`: most-distant-RBH threshold.

`ar` uses a correct reverse-direction join, fixing the V1 direction/indexing
bug rather than reproducing it.

## Canonicalization and symmetrization

Each retained direction is joined under `u=min(query,target)` and
`v=max(query,target)`. For LRB, assemble complete directional
`Wuv=(Cuv+Cvu)*Buv`, where C indicates cutoff passage and B is the complete
normalized matrix. A direction that failed its own cutoff can still contribute
through reverse C; truly absent hits have score zero. Both C directions passed
means a multiplier of 2. `score_uv` and `score_vu` store W, not filtered B.
Non-LRB experiment modes retain their selected-B assembly behavior.

Supported projections are `forward`, `max`, `min`, `mean`, and `geometric_mean`.
The default is `mean`: `(Wuv+Wvu)/2` for LRB. Leiden consumes an undirected
graph, whereas OF MCL consumes directional W; this is a deliberate projection,
not exact MCL equivalence. Mean uses full precision before projection, not OF
writer's three-decimal text quantization. Explicit `forward` uses score_uv
with reverse fallback; other projection choices remain available.

The retained edge schema is:

```text
u, v                         int64
u_species, v_species         int32
score_uv, score_vu, weight   float64
coverage                     float32
edge_type                    string
```

The manifest binds algorithm version, parameters, input hashes, output hashes,
and row counts. Identical verified artifacts are reused; changed inputs,
parameters, or output hashes trigger reconstruction.


The repaired version is `of315-nbs-lrb-directional-mean-v4`; pre-fix edge
caches are invalidated even with unchanged inputs and explicit parameters.
