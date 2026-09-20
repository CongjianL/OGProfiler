#!/usr/bin/env bash
# S0d: three-way artifact comparison (V1 ↔ V2 ↔ OrthoFinder3).
#
# Requires s0a/s0b/s0c to have completed in the shared results root.
# Part of the 5-submission S0 workflow (s0p -> s0a/s0b/s0c -> s0d).

#SBATCH --job-name=ogp-s0d-compare
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=04:00:00

set -euo pipefail

: "${DEV_SOURCE_DIR:?Missing immutable source snapshot}"
: "${DEV_RUN_DIR:?Missing immutable run directory}"
: "${DEV_CONDA:?Missing micromamba executable from project remote config}"
: "${DEV_CONDA_ENV:?Missing micromamba environment from project remote config}"

ENV="$DEV_CONDA_ENV"
DATA_NAME="real_embleya"
SHARED_ROOT="$(dirname "$DEV_RUN_DIR")/s0_v1_v2_of_regression"
TASK_ROOT="$SHARED_ROOT/results/$DATA_NAME"
V1_ROOT="$TASK_ROOT/v1_full"
V2_ROOT="$TASK_ROOT/v2_edges"
OF_WORKDIR="$(cat "$TASK_ROOT/of_workdir.txt" 2>/dev/null || true)"
COMPARE_ROOT="$TASK_ROOT/comparison"
mkdir -p "$COMPARE_ROOT"

[[ -d "$V1_ROOT/WorkingDirectory" ]] || { echo "Missing V1 result: $V1_ROOT/WorkingDirectory" >&2; exit 3; }
[[ -d "$V2_ROOT" ]] || { echo "Missing V2 result: $V2_ROOT" >&2; exit 3; }
[[ -n "$OF_WORKDIR" && -d "$OF_WORKDIR" ]] || {
  echo "Missing OrthoFinder WorkingDirectory (run s0c first)" >&2
  exit 3
}

echo "==> [$DATA_NAME] Compare V1 / V2 / OrthoFinder3"
"$DEV_CONDA" run -n "$ENV" python \
  "$DEV_SOURCE_DIR/benchmarks/compare_v1_v2_orthofinder.py" \
    --v1-working-dir "$V1_ROOT/WorkingDirectory" \
    --v2-run "$V2_ROOT" \
    --of-working-dir "$OF_WORKDIR" \
    --out "$COMPARE_ROOT/report.json" \
  > "$COMPARE_ROOT/report.stdout" 2> "$COMPARE_ROOT/report.stderr"

echo "==> S0 comparison complete: $COMPARE_ROOT/report.json"
