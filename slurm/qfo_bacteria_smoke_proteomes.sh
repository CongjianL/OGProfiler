#!/usr/bin/env bash
#SBATCH --job-name=ogp-qfo-bacteria-3proteomes
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail
exec bash "${DEV_SOURCE_DIR:?}/benchmarks/qfo/run_methods.sh" smoke-proteomes "$@"
