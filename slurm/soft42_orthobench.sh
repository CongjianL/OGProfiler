#!/usr/bin/env bash
#SBATCH --job-name=ogp-soft42-orthobench
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=56
#SBATCH --mem=250G
#SBATCH --time=72:00:00
set -euo pipefail
: "${DEV_SOURCE_DIR:?Missing immutable source}"
: "${DEV_RUN_DIR:?Missing independent run directory}"
: "${DEV_CONDA:?Missing configured environment manager}"
: "${DEV_CONDA_ENV:?Missing configured environment}"
ORIGIN=${1:?Provide original P6 benchmark run}
if [[ ${OGP_H5_READY:-0} != 1 ]]; then
    exec "$DEV_CONDA" run -n "$DEV_CONDA_ENV" env OGP_H5_READY=1 bash "$0" "$@"
fi
export PYTHONPATH="$DEV_SOURCE_DIR/src:$DEV_SOURCE_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python -m benchmarks.og_extraction.soft42_orthobench prepare \
    --run "$DEV_RUN_DIR" --origin "$ORIGIN" \
    --preset "$DEV_SOURCE_DIR/presets/embleya-soft42.yaml" --workers 56
"$DEV_CONDA" list -n "$DEV_CONDA_ENV" > "$DEV_RUN_DIR/environment.txt"
echo '==> H5: all components freshly computed, frozen mean SSN, selected soft42'
set +e
/usr/bin/time -v -o "$DEV_RUN_DIR/h5-hierarchy.time" \
    python -m ogprofiler hierarchy-all --run "$DEV_RUN_DIR/new-hierarchy" \
    > "$DEV_RUN_DIR/h5-hierarchy.stdout" 2> "$DEV_RUN_DIR/h5-hierarchy.stderr"
H5_EXIT=$?
set -e
printf '%s\n' "$H5_EXIT" > "$DEV_RUN_DIR/h5-hierarchy-exit-code.txt"
python -m benchmarks.og_extraction.hierarchy_regression summary \
    --run "$DEV_RUN_DIR/new-hierarchy" --out "$DEV_RUN_DIR/h5-hierarchy-summary.json"
if [[ "$H5_EXIT" != 0 ]]; then
    printf '{"h5_hierarchy_exit":%s,"h5_status":"BLOCKED_BY_HIERARCHY","scoring_started":false}\n' "$H5_EXIT" > "$DEV_RUN_DIR/h4-h5-completion.json"
    exit "$H5_EXIT"
fi
python -m benchmarks.og_extraction.soft42_orthobench freeze --run "$DEV_RUN_DIR"
python -m ogprofiler annotate-network --run "$DEV_RUN_DIR/new-hierarchy" \
    > "$DEV_RUN_DIR/h5-events.stdout" 2> "$DEV_RUN_DIR/h5-events.stderr"
/usr/bin/time -v -o "$DEV_RUN_DIR/h5-parity.time" \
    python -m benchmarks.og_extraction.fixed_hierarchy --run "$DEV_RUN_DIR/new-hierarchy" \
    --out "$DEV_RUN_DIR/h5-parity" > "$DEV_RUN_DIR/h5-parity.stdout" 2> "$DEV_RUN_DIR/h5-parity.stderr"
python -m ogprofiler export --run "$DEV_RUN_DIR/new-hierarchy" --strategy terminal \
    > "$DEV_RUN_DIR/h5-terminal.stdout" 2> "$DEV_RUN_DIR/h5-terminal.stderr"
/usr/bin/time -v -o "$DEV_RUN_DIR/h5-scoring.time" \
    python -m benchmarks.og_extraction.hierarchy_regression score --run "$DEV_RUN_DIR" --origin "$ORIGIN" \
    > "$DEV_RUN_DIR/h5-scoring.stdout" 2> "$DEV_RUN_DIR/h5-scoring.stderr"
printf '{"fresh_full_hierarchy":true,"h5_status":"EVALUATED","accuracy_improvement_claimed":false}\n' > "$DEV_RUN_DIR/h4-h5-completion.json"
