#!/usr/bin/env bash
#SBATCH --job-name=b1-proteinortho6-post
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail
SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
OB=/home/mselab/licj/project_data/OGProfiler/Open_Orthobench-master
MM=/home/mselab/licj/micromamba/bin/micromamba
OUT=$ROOT/03_runs/B1_orthobench/Proteinortho6/rep1_retry1
GROUPS_TSV=$ROOT/04_standardized/groups/Proteinortho6/OpenOrthobench/groups_rep1.tsv
METRICS=$ROOT/05_metrics/orthobench/by_method/Proteinortho6
POST=$OUT/postprocess_retry1
[[ -f "$GROUPS_TSV" && -f "$OUT/raw_tool/metadata.tsv" && -f "$OUT/validation.log" ]]
grep -q $'^status\tOK$' "$OUT/raw_tool/metadata.tsv"
grep -q '^OK$' "$OUT/validation.log"
[[ ! -e "$POST" ]] || { echo "Refusing existing postprocess retry directory: $POST" >&2; exit 3; }
[[ ! -e "$METRICS" ]] || { echo "Refusing existing metrics directory: $METRICS" >&2; exit 3; }
mkdir -p "$POST" "$METRICS"
awk -F '\t' 'NR>1 {sub(/\r$/, "", $2); a[$1]=a[$1] " " $2} END {for (x in a) print x ":" a[x]}' "$GROUPS_TSV" | LC_ALL=C sort > "$POST/prediction.txt"
cat > "$POST/MANIFEST.tsv" <<EOF
postprocess_of_slurm_job_id\t1404268
failed_official_slurm_job_id\t1404268
failure\tCRLF_in_standardized_TSV_became_universal-newline_group_splits
slurm_job_id\t${SLURM_JOB_ID:?}
git_commit\t${DEV_GIT_COMMIT:?}
source_sha256\t${DEV_SOURCE_HASH:?}
EOF
"$SRC/workflows/run_timed.sh" OpenOrthobench Proteinortho6 1 "$POST/raw_official" -- "$MM" run -n ogprofiler python "$OB/BENCHMARKS/benchmark.py" "$POST/prediction.txt"
"$MM" run -n ogprofiler python "$SRC/scripts/metrics/parse_open_orthobench_official.py" "$POST/raw_official/stdout.log" --out "$METRICS/official.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/metrics/orthobench_metrics.py" --refogs "$OB/BENCHMARKS/RefOGs" --groups "$GROUPS_TSV" --out "$METRICS/refog_metrics.tsv" --summary "$METRICS/extended_summary_raw.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/diagnostics/orthobench_pairwise_attribution.py" --groups "$GROUPS_TSV" --refogs "$OB/BENCHMARKS/RefOGs" --outdir "$METRICS"
[[ $(($(wc -l < "$METRICS/refog_metrics.tsv")-1)) -eq 70 ]]
