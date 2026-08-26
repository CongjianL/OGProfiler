#!/usr/bin/env bash
# Real OGProfiler 2 CLI wrapper: configuration overrides carry threads and seed.
set -euo pipefail
[[ $# -eq 4 ]] || { echo "Usage: $0 INPUT OUTPUT THREADS SEED" >&2; exit 2; }
INPUT=$1; OUTPUT=$2; THREADS=$3; SEED=$4
ROOT=$(cd "$(dirname "$0")/.." && pwd -P)
SOURCE_ROOT=$(cd "$ROOT/.." && pwd -P)
RUN_DIR="$OUTPUT/raw_ogprofiler"
mkdir -p "$OUTPUT"
export BENCHMARK_THREADS="$THREADS" BENCHMARK_SEED="$SEED"
"$ROOT/workflows/run_timed.sh" OGProfiler2 OpenOrthobench 1 "$RUN_DIR" -- \
  env "PYTHONPATH=$SOURCE_ROOT/src${PYTHONPATH:+:$PYTHONPATH}" python -m ogprofiler run --proteomes "$INPUT" --out "$OUTPUT/run" \
  --set "search.threads=$THREADS" --set "runtime.workers=$THREADS" --set "hierarchy.seed=$SEED"
