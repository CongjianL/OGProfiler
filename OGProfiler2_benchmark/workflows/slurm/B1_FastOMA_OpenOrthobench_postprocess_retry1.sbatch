#!/usr/bin/env bash
#SBATCH --job-name=b1-fastoma-post
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
set -euo pipefail

SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
OB=/home/mselab/licj/project_data/OGProfiler/Open_Orthobench-master
MM=/home/mselab/licj/micromamba/bin/micromamba
OUT=$ROOT/03_runs/B1_orthobench/FastOMA/rep1_retry1
PRIMARY=$OUT/result/RootHOGs.tsv
GROUPS_TSV=$ROOT/04_standardized/groups/FastOMA/OpenOrthobench/groups_rep1.tsv
METRICS=$ROOT/05_metrics/orthobench/by_method/FastOMA
POST=$OUT/postprocess_retry1

[[ -f "$PRIMARY" && -d "$OUT/input/proteome" ]] || { echo "Completed FastOMA formal artifacts are absent" >&2; exit 4; }
[[ -f "$OUT/raw_tool/metadata.tsv" ]] || { echo "FastOMA timing metadata is absent" >&2; exit 4; }
grep -q $'^status\tOK$' "$OUT/raw_tool/metadata.tsv"
[[ "$(head -n 1 "$PRIMARY")" == $'RootHOG\tProtein\tOMAmerRootHOG' ]] || {
 echo "Unexpected FastOMA RootHOGs.tsv schema" >&2; exit 5;
}
[[ $(($(wc -l < "$PRIMARY") - 1)) -eq 197736 ]] || { echo "Unexpected preserved RootHOG assignment count" >&2; exit 5; }
[[ ! -e "$POST" ]] || { echo "Refusing existing postprocess directory: $POST" >&2; exit 3; }
[[ ! -e "$GROUPS_TSV" ]] || { echo "Refusing existing standardized output: $GROUPS_TSV" >&2; exit 3; }
[[ ! -e "$METRICS" ]] || { echo "Refusing existing metric output: $METRICS" >&2; exit 3; }
mkdir -p "$POST"

cat > "$POST/MANIFEST.tsv" <<EOF
postprocess_of_slurm_job_id\t1404543
inference_status\tcompleted_successfully
inference_wall_clock\t1:36:14
primary_schema\tRootHOG_Protein_OMAmerRootHOG
primary_assignment_rows\t197736
primary_root_hogs\t19472
failure\tgeneric_converter_treated_OMAmerRootHOG_annotations_as_proteins
slurm_job_id\t${SLURM_JOB_ID:?}
git_commit\t${DEV_GIT_COMMIT:?}
source_sha256\t${DEV_SOURCE_HASH:?}
EOF
sha256sum "$PRIMARY" "$OUT/result/OrthologousGroups.tsv" "$OUT/result/FastOMA_HOGs.orthoxml" > "$POST/PRIMARY_SHA256SUMS"
sha256sum "$OUT"/input/proteome/*.fa "$OUT/input/species_tree.nwk" > "$POST/INPUT_SHA256SUMS"

"$MM" run -n ogprofiler python "$SRC/scripts/converters/convert_competitor_groups.py" \
 --tool fastoma --input "$PRIMARY" --fasta "$OUT/input/proteome" --out "$GROUPS_TSV"
"$MM" run -n ogprofiler python "$SRC/scripts/converters/validate_groups.py" \
 "$GROUPS_TSV" --fasta "$OUT/input/proteome" > "$POST/validation.log"
awk -F '\t' 'NR>1 {sub(/\r$/, "", $2); a[$1]=a[$1] " " $2} END {for (x in a) print x ":" a[x]}' \
 "$GROUPS_TSV" | LC_ALL=C sort > "$POST/prediction.txt"

"$SRC/workflows/run_timed.sh" OpenOrthobench FastOMA 1 "$POST/raw_official" -- \
 "$MM" run -n ogprofiler python "$OB/BENCHMARKS/benchmark.py" "$POST/prediction.txt"
mkdir -p "$METRICS"
"$MM" run -n ogprofiler python "$SRC/scripts/metrics/parse_open_orthobench_official.py" \
 "$POST/raw_official/stdout.log" --out "$METRICS/official.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/metrics/orthobench_metrics.py" \
 --refogs "$OB/BENCHMARKS/RefOGs" --groups "$GROUPS_TSV" \
 --out "$METRICS/refog_metrics.tsv" --summary "$METRICS/extended_summary_raw.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/diagnostics/orthobench_pairwise_attribution.py" \
 --groups "$GROUPS_TSV" --refogs "$OB/BENCHMARKS/RefOGs" --outdir "$METRICS"
[[ $(($(wc -l < "$METRICS/refog_metrics.tsv") - 1)) -eq 70 ]]
printf 'primary_output\t%s\ngroups\t%s\nstatus\tPASS\n' "$PRIMARY" "$GROUPS_TSV" > "$POST/FORMAL_STATUS.tsv"
