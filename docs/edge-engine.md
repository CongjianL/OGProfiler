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

## Legacy NBS

Normalization is performed separately for every directed species pair. Self
hits are removed. For each group, V2 reproduces frozen V1 behavior:

1. sort hits by query-length × target-length;
2. for fewer than 100 hits, use the complete group;
3. otherwise use V1 bin widths (20, 200, or 1000), inclusive `scale + 1`
   slices, and the linearly interpolated 95th bitscore percentile;
4. fit `log10(bitscore) = a * log10(length_product) + b`;
5. calculate `normalized_score = bitscore / (10**b * length_product**a)`.

The logarithmic model is linear, so V2 uses the closed-form least-squares
solution rather than a nonlinear optimizer. Degenerate groups (one usable
point, degenerate length products, or a non-finite fit) are controlled by
`similarity.nbs_fallback`:

- `v1_zero` (default): drop the group's hits, matching V1 and OrthoFinder3's
  all-zero matrix for too-few-hit groups;
- `v2_max`: `a=0`, `b=log10(max_bitscore)`, mapping the group maximum to 1.

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

LRB assigns each query its lowest outgoing RBH score and retains every hit at
or above that threshold, including qualifying same-species paralogs. The
legacy no-RBH fallback is `best_external_score + 1e-6`; consequently a query
with external hits but no RBH retains none of those external hits. A query
with only same-species hits uses `1e-6`.

Compatibility modes are explicit:

- `rbh` and `arb`: reciprocal best hits;
- `ar`: any hit for which the reverse direction exists;
- `lrb`: most-distant-RBH threshold.

`ar` uses a correct reverse-direction join, fixing the V1 direction/indexing
bug rather than reproducing it.

## Canonicalization and symmetrization

Each retained direction is joined under `u=min(query,target)` and
`v=max(query,target)`. Missing directions have score zero. `score_uv` and
`score_vu` therefore have stable meanings independent of input order.

Supported weights are `forward`, `max`, `min`, `mean`, and `geometric_mean`.
The default is `forward`, matching OrthoFinder3's directional edge weight
(`connect2 × B[i→j]`): use `score_uv` (u→v), falling back to `score_vu` when
only the reverse direction is retained. `max`/`min`/`mean`/`geometric_mean`
remain available and are subject to benchmark selection.

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

