#!/usr/bin/env bash
# Official revised Orthobench interface verified from BENCHMARKS/benchmark.py.
set -euo pipefail
[[ $# -eq 3 ]] || { echo "Usage: $0 OGPROFILER_RUN ORTHOBENCH_ROOT OUTPUT_DIR" >&2; exit 2; }
RUN=$1; ROOT=$2; OUT=$3; HERE=$(cd "$(dirname "$0")/.." && pwd -P)
mkdir -p "$OUT" "$HERE/04_standardized/groups/OGProfiler2/OpenOrthobench"
GROUPS="$HERE/04_standardized/groups/OGProfiler2/OpenOrthobench/groups.tsv"
python "$HERE/scripts/converters/convert_ogprofiler.py" --run "$RUN" --groups-out "$GROUPS" --fasta "$ROOT/BENCHMARKS/Input"
# benchmark.py accepts a one-group-per-line file. Remove the tabular header and preserve group boundaries.
awk -F '\t' 'NR>1 {a[$1]=a[$1] " " $2} END {for (x in a) print x ":" a[x]}' "$GROUPS" | LC_ALL=C sort > "$OUT/orthogroups_results_file.txt"
"$HERE/workflows/run_timed.sh" OpenOrthobench OpenOrthobench 1 "$OUT/raw_official" -- python "$ROOT/BENCHMARKS/benchmark.py" "$OUT/orthogroups_results_file.txt"
python "$HERE/scripts/metrics/orthobench_metrics.py" --refogs "$ROOT/BENCHMARKS/RefOGs" --groups "$GROUPS" --out "$HERE/05_metrics/orthobench/OGProfiler2_refog_metrics.tsv" --summary "$HERE/05_metrics/orthobench/OGProfiler2_extended_summary.tsv"
