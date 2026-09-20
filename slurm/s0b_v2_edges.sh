#!/usr/bin/env bash
# S0b: V2 prepare -> search -> edges (the aligned upstream stages only).
#
# Part of the 5-submission S0 workflow (s0p -> s0a/s0b/s0c -> s0d); writes into the shared results root.

#SBATCH --job-name=ogp-s0b-v2
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=56
#SBATCH --mem=128G
#SBATCH --time=24:00:00

set -euo pipefail

: "${DEV_SOURCE_DIR:?Missing immutable source snapshot}"
: "${DEV_RUN_DIR:?Missing immutable run directory}"
: "${DEV_CONDA:?Missing micromamba executable from project remote config}"
: "${DEV_CONDA_ENV:?Missing micromamba environment from project remote config}"

ENV="$DEV_CONDA_ENV"
DATA_NAME="real_embleya"
SHARED_ROOT="$(dirname "$DEV_RUN_DIR")/s0_v1_v2_of_regression"
INPUT_DIR="$SHARED_ROOT/proteomes"
TASK_ROOT="$SHARED_ROOT/results/$DATA_NAME"
V2_ROOT="$TASK_ROOT/v2_edges"
mkdir -p "$TASK_ROOT"

"$DEV_CONDA" list -n "$ENV" > "$TASK_ROOT/v2-conda-list.txt" 2>&1 || true

echo "==> [$DATA_NAME] V2 prepare/search/edges (env=$ENV)"
/usr/bin/time -v -o "$TASK_ROOT/v2.time" \
  "$DEV_CONDA" run -n "$ENV" bash -lc \
    "cd '$DEV_SOURCE_DIR' && PYTHONPATH=src python -m ogprofiler run \
      --proteomes '$INPUT_DIR' --out '$V2_ROOT' --until-stage edges \
      --set search.threads=$SLURM_CPUS_PER_TASK" \
  > "$TASK_ROOT/v2.stdout" 2> "$TASK_ROOT/v2.stderr"

echo "==> V2 edges complete: $V2_ROOT"
