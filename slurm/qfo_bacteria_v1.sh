#!/usr/bin/env bash
#SBATCH --job-name=ogp-qfo-bacteria-v1
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
set -euo pipefail
exec env QFO_METHOD=v1 bash "${DEV_SOURCE_DIR:?}/benchmarks/qfo/run_methods.sh" full "$@"
