#!/usr/bin/env bash
#SBATCH --job-name=ogp-qfo-bacteria-preflight
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail
: "${DEV_SOURCE_DIR:?}" "${DEV_RUN_DIR:?}" "${DEV_CONDA:?}" "${DEV_CONDA_ENV:?}"
BACTERIA_DIR=${1:?Provide independent QFO bacteria directory}
OF_ENV=${2:?Provide existing OrthoFinder environment name}
V1_ENV=${3:?Provide existing frozen V1 environment name}
PREVIOUS_RUN=${4:-}
printf 'supersedes_job\t%s\nscope_change\tuser selected bacteria only; all excluded\n' "$PREVIOUS_RUN" > "$DEV_RUN_DIR/scope-change-provenance.tsv"
if [[ ${QFO_PREFLIGHT_READY:-0} != 1 ]]; then
    exec "$DEV_CONDA" run -n "$DEV_CONDA_ENV" env QFO_PREFLIGHT_READY=1 bash "$0" "$@"
fi
export PYTHONPATH="$DEV_SOURCE_DIR/src:$DEV_SOURCE_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export BENCHMARK_THREADS=4 BENCHMARK_SEED=42
TIMED="$DEV_SOURCE_DIR/OGProfiler2_benchmark/workflows/run_timed.sh"
"$DEV_CONDA" list -n "$DEV_CONDA_ENV" > "$DEV_RUN_DIR/environment.txt"
"$TIMED" QFO_INPUT_AUDIT QFO_bacteria 1 "$DEV_RUN_DIR/audit-timing" -- \
    python -m benchmarks.qfo.preflight --source "$BACTERIA_DIR" --collection bacteria --out "$DEV_RUN_DIR/preflight"
"$TIMED" OGProfiler2_prepare QFO_bacteria 1 "$DEV_RUN_DIR/prepare-timing" -- \
    python -m ogprofiler prepare --proteomes "$DEV_RUN_DIR/preflight/input/bacteria" \
    --out "$DEV_RUN_DIR/prepared-bacteria" --config "$DEV_RUN_DIR/preflight/soft42-default.yaml"
"$TIMED" OGProfiler2_search_smoke QFO_3x20_smoke 1 "$DEV_RUN_DIR/smoke-prepare-timing" -- \
    python -m ogprofiler prepare --proteomes "$DEV_RUN_DIR/preflight/smoke-input" \
    --out "$DEV_RUN_DIR/smoke-run" --set search.threads=4
"$TIMED" OGProfiler2_search_smoke QFO_3x20_smoke 1 "$DEV_RUN_DIR/smoke-search-timing" -- \
    python -m ogprofiler search --run "$DEV_RUN_DIR/smoke-run"
"$TIMED" OrthoFinder3_environment QFO_bacteria 1 "$DEV_RUN_DIR/of-environment-timing" -- \
    "$DEV_CONDA" run -n "$OF_ENV" orthofinder -h
"$TIMED" OGProfiler1First_environment QFO_bacteria 1 "$DEV_RUN_DIR/v1-environment-timing" -- \
    "$DEV_CONDA" run -n "$V1_ENV" python -c 'import igraph,leidenalg,numpy,progressbar,pyfasta,scipy; print("V1 environment imports OK")'
printf '{"preflight_completed":true,"full_scientific_runs_started":false,"official_qfo_scores_produced":false}\n' \
    > "$DEV_RUN_DIR/qfo-preflight-completion.json"
