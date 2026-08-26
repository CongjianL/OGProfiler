# OGProfiler 1 baseline outputs

For each dataset A–E, capture the following in `<dataset>/metrics.json`:

- input manifest SHA-256;
- exact V1 command and dependency/tool versions;
- SSN node, edge, and connected-component counts;
- hierarchy node count and depth;
- terminal family membership and count;
- ortholog pair count;
- wall time, peak RSS, and disk usage.

The frozen implementation is `../../legacy/OGProfiler_v1.py`. Baseline runs
must use immutable input files from `../datasets/` and preserve raw result
checksums even when the raw files are not committed.

Formal runs are stored in timestamped directories. Each run contains a
machine-readable `summary.json`, tabular `summary.tsv`, human-readable
`REPORT.md`, and per-dataset normalized hierarchy comparisons. Regenerate a
run report with `../aggregate_baseline_report.py`; the run identity and source
and dataset hashes are required arguments so provenance is never inferred.
