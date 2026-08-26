#!/usr/bin/env bash
#SBATCH --job-name=ogp-phase13
#SBATCH --array=0-15%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00

set -euo pipefail

OGP_SOURCE="${OGP_SOURCE:-${DEV_SOURCE_DIR:?Missing immutable source snapshot}}"
OGP_RUN_ROOT="${OGP_RUN_ROOT:-${DEV_RUN_DIR:?Missing Phase 13 run directory}}"
OGP_CONDA="${OGP_CONDA:-${DEV_CONDA:?Missing Conda/Micromamba executable}}"
OGP_CONDA_ENV="${OGP_CONDA_ENV:-${DEV_CONDA_ENV:?Missing project environment}}"

cd "$OGP_SOURCE"
export PYTHONPATH="$OGP_SOURCE/src"
for INDEX in $(seq "$SLURM_ARRAY_TASK_ID" 16 242); do
  "$OGP_CONDA" run -n "$OGP_CONDA_ENV" python \
    benchmarks/run_phase13_scenario.py \
    --out "$OGP_RUN_ROOT/results" \
    --index "$INDEX" \
    --workers "$SLURM_CPUS_PER_TASK"
done
