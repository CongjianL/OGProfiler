#!/usr/bin/env bash
#SBATCH --job-name=ogp-phase12
#SBATCH --array=0-15%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/phase12-%A_%a.out
#SBATCH --error=logs/phase12-%A_%a.err

set -euo pipefail

OGP_SOURCE="${OGP_SOURCE:-${DEV_SOURCE_DIR:?Missing immutable source snapshot}}"
OGP_RUN_ROOT="${OGP_RUN_ROOT:-${DEV_RUN_DIR:?Missing Phase 12 run directory}}"
OGP_PROTEOMES="${OGP_PROTEOMES:-$OGP_SOURCE/benchmarks/datasets/A_small_sanity/proteomes}"
OGP_HITS="${OGP_HITS:-$OGP_SOURCE/benchmarks/phase3_search/20260826T0011Z_dataset_a/hits.parquet}"
OGP_GROUND_TRUTH="${OGP_GROUND_TRUTH:-$OGP_SOURCE/benchmarks/datasets/A_small_sanity/ground_truth.tsv}"
OGP_CONDA="${OGP_CONDA:-${DEV_CONDA:?Missing Conda/Micromamba executable}}"
OGP_CONDA_ENV="${OGP_CONDA_ENV:-${DEV_CONDA_ENV:?Missing project environment}}"

mkdir -p "$OGP_RUN_ROOT/logs"
cd "$OGP_SOURCE"
export PYTHONPATH="$OGP_SOURCE/src"
"$OGP_CONDA" run -n "$OGP_CONDA_ENV" python \
  benchmarks/run_phase12_matrix.py \
  --proteomes "$OGP_PROTEOMES" \
  --hits "$OGP_HITS" \
  --ground-truth "$OGP_GROUND_TRUTH" \
  --dataset A_small_sanity \
  --out "$OGP_RUN_ROOT/results" \
  --index "$SLURM_ARRAY_TASK_ID" \
  --workers "$SLURM_CPUS_PER_TASK"
