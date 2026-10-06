#!/usr/bin/env bash
# Called by resource-specific Slurm entrypoints; no remote infrastructure literals.
set -euo pipefail
MODE=${1:?}; shift
ORIGIN=${1:?}; DIGEST=${2:?}; OF_ENV=${3:?}; V1_ENV=${4:?}
: "${DEV_SOURCE_DIR:?}" "${DEV_RUN_DIR:?}" "${DEV_CONDA:?}" "${DEV_CONDA_ENV:?}" "${SLURM_CPUS_PER_TASK:?}"
if [[ ${QFO_METHODS_READY:-0} != 1 ]]; then
    exec "$DEV_CONDA" run -n "$DEV_CONDA_ENV" env QFO_METHODS_READY=1 bash "$0" "$MODE" "$@"
fi
export PYTHONPATH="$DEV_SOURCE_DIR/src:$DEV_SOURCE_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export BENCHMARK_THREADS="$SLURM_CPUS_PER_TASK" BENCHMARK_SEED=42
B="$DEV_SOURCE_DIR/OGProfiler2_benchmark"
TIMED="$B/workflows/run_timed.sh"
cd "$DEV_RUN_DIR"
ROOT="$DEV_RUN_DIR/campaign"
DATASET="QFO_bacteria_${MODE}"
"$TIMED" QFO_stage "$DATASET" 1 "$DEV_RUN_DIR/stage-timing" -- \
    python -m benchmarks.qfo.campaign --origin "$ORIGIN" --out "$ROOT" \
    --expected-digest "$DIGEST" --mode "$MODE"
for env_name in "$DEV_CONDA_ENV" "$OF_ENV" "$V1_ENV"; do
    "$DEV_CONDA" list -n "$env_name" --explicit > "$DEV_RUN_DIR/environment-$env_name.txt"
done
[[ $(sha256sum "$B/00_env/ogprofiler_v1_first/OGProfiler.py" | awk '{print $1}') == 65ec43d269b410956fadc4f8215a3e0a0772eac8bcc1060ab011b090823cca41 ]]
V1_PREFIX=$("$DEV_CONDA" run -n "$V1_ENV" python -c 'import sys; print(sys.prefix)')
[[ -x "$V1_PREFIX/bin/python" && -x "$V1_PREFIX/bin/diamond" ]]
"$TIMED" OGProfiler2_soft42 "$DATASET" 1 "$ROOT/v2-timing" -- \
    python -m ogprofiler run --proteomes "$ROOT/input" --out "$ROOT/v2" \
    --config "$ROOT/soft42-default.yaml" --set "search.threads=$BENCHMARK_THREADS" \
    --set "runtime.workers=$BENCHMARK_THREADS"
"$TIMED" OrthoFinder3 "$DATASET" 1 "$ROOT/of3-timing" -- \
    "$DEV_CONDA" run -n "$OF_ENV" orthofinder -f "$ROOT/input" -o "$ROOT/of3" \
    -S diamond -t "$BENCHMARK_THREADS" -a "$BENCHMARK_THREADS" -og
OF_GROUPS=$(find "$ROOT/of3" -type f -path '*/Orthogroups/Orthogroups.tsv')
[[ -n "$OF_GROUPS" && $(printf '%s\n' "$OF_GROUPS" | wc -l) -eq 1 ]]
# Keep the assigned primary partition: converter-added unassigned singletons
# are for complete coverage validation only, not OF3 pairwise reference labels.
printf '%s\n' "$OF_GROUPS" > "$ROOT/of3-assigned-primary.txt"
"$TIMED" OrthoFinder3_conversion "$DATASET" 1 "$ROOT/of3-conversion-timing" -- \
    python "$B/scripts/converters/convert_competitor_groups.py" --tool orthofinder \
    --input "$OF_GROUPS" --fasta "$ROOT/input" --out "$ROOT/of3-groups.tsv"
export OGPROFILER_V1_DIAMOND_REAL="$V1_PREFIX/bin/diamond"
V1_PATH="$B/02_configs/ogprofiler_v1_first/bin:$V1_PREFIX/bin:$PATH"
BENCHMARK_SEED=unavailable "$TIMED" OGProfiler1First "$DATASET" 1 "$ROOT/v1-timing" -- \
    env "PATH=$V1_PATH" "$V1_PREFIX/bin/python" "$B/00_env/ogprofiler_v1_first/OGProfiler.py" \
    -i "$ROOT/v1-input" -o "$ROOT/v1.gml" -s diamond -e 1e-5 \
    -t "$BENCHMARK_THREADS" -d lrb -w NBS -m rber -a "$BENCHMARK_THREADS" -g 1.0 --so 0
PRIMARY="$ROOT/v1-input/WorkingDirectory/OGFile_coalescence_SameGenome.txt"
[[ -s "$PRIMARY" ]]
"$TIMED" OGProfiler1First_conversion "$DATASET" 1 "$ROOT/v1-conversion-timing" -- \
    python "$B/scripts/converters/convert_competitor_groups.py" --tool ogprofiler_v1 \
    --input "$PRIMARY" --id-map "$ROOT/v1-id-map.tsv" --fasta "$ROOT/input" --out "$ROOT/v1-groups.tsv"
printf '{"methods_completed":true,"mode":"%s","official_qfo_scores_produced":false,"comparison_evaluated":false}\n' "$MODE" \
    > "$DEV_RUN_DIR/qfo-methods-completion.json"
