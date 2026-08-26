# Hierarchy prototype benchmark witnesses

The prototype consumes each dataset's deterministic `reference_ssn.gml` via:

```bash
ogprofiler prototype-hierarchy \
  --ssn benchmarks/datasets/E_large_connected_component/reference_ssn.gml \
  --out RUN_DIRECTORY \
  --set hierarchy.gamma_min=0.001 \
  --set hierarchy.gamma_max=2.0 \
  --set hierarchy.min_family_size=2 \
  --set hierarchy.max_child_fraction=0.9
```

Local development witness on 2026-08-25 for Dataset E:

| Metric | Value |
| --- | ---: |
| vertices | 240 |
| edges | 4,685 |
| connected components | 1 |
| hierarchy nodes | 7 |
| terminal families | 6 |
| Leiden calls | 67 |
| subgraph constructions | 6 |
| runtime | 0.101 s |
| observed peak RSS | 61,931,520 bytes |

The six terminal families match the six planted subfamilies. Runtime and RSS
are machine-specific witnesses; hierarchy counts and memberships are the
regression surface. A formal V1/V2 comparison is recorded only after the V1
baseline SSN is generated with the frozen environment and command manifest.

Remote `campus-server` witness in the declared `ogprofiler` environment on the
same date preserved the same structural metrics (7 nodes, 6 terminal families,
67 Leiden calls, 6 subgraph constructions) and observed 0.0478 s runtime with
64,888,832 bytes peak RSS.
