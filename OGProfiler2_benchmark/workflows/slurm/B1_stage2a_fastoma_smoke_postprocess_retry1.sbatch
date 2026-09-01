#!/usr/bin/env bash
#SBATCH --job-name=b1-fastoma-smoke-post
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -euo pipefail

SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
MM=/home/mselab/licj/micromamba/bin/micromamba
OUT=$ROOT/03_runs/B1_orthobench/competitor_smoke_stage2a/fastoma-0.5.1_slurm3
PRIMARY=$OUT/result/RootHOGs.tsv
POST=$OUT/postprocess_retry1

[[ -f "$PRIMARY" && -d "$OUT/input/proteome" ]] || { echo "Completed FastOMA smoke artifacts are absent" >&2; exit 4; }
[[ -f "$OUT/result/OrthologousGroups.tsv" && -f "$OUT/result/FastOMA_HOGs.orthoxml" ]] || {
 echo "FastOMA companion primary artifacts are absent" >&2; exit 4;
}
[[ -f "$OUT/timing/metadata.tsv" ]] || { echo "FastOMA timing metadata is absent" >&2; exit 4; }
grep -q $'^status\tFAILED$' "$OUT/timing/metadata.tsv"
grep -q "KeyError: 'is_roothog'" "$OUT/timing/stdout.log"
[[ "$(head -n 1 "$PRIMARY")" == $'RootHOG\tProtein\tOMAmerRootHOG' ]] || {
 echo "Unexpected FastOMA RootHOGs.tsv schema" >&2; exit 5;
}
[[ $(wc -l < "$PRIMARY") -eq 1 ]] || { echo "Expected the audited header-only smoke primary" >&2; exit 5; }
[[ $(find "$OUT/input/proteome" -maxdepth 1 -type f -name '*.fa' -print0 | xargs -0 grep -h '^>' | wc -l) -eq 30 ]] || {
 echo "Expected exactly 30 proteins in the preserved smoke input" >&2; exit 5;
}
[[ ! -e "$POST" ]] || { echo "Refusing existing postprocess retry directory: $POST" >&2; exit 3; }
mkdir -p "$POST"

cat > "$POST/MANIFEST.tsv" <<EOF
postprocess_of_slurm_job_id\t1404535
primary_inference_status\tcompleted_successfully
slurm_terminal_status\tfailed_in_optional_report
failure\tFastOMA_report_KeyError_is_roothog_on_empty_assignment_table
primary_rows_excluding_header\t0
unassigned_policy\tdeterministic_singleton_per_input_protein
slurm_job_id\t${SLURM_JOB_ID:?}
git_commit\t${DEV_GIT_COMMIT:?}
source_sha256\t${DEV_SOURCE_HASH:?}
EOF
sha256sum "$PRIMARY" "$OUT/result/OrthologousGroups.tsv" "$OUT/result/FastOMA_HOGs.orthoxml" > "$POST/PRIMARY_SHA256SUMS"
sha256sum "$OUT"/input/proteome/*.fa "$OUT/input/species_tree.nwk" > "$POST/INPUT_SHA256SUMS"

"$MM" run -n ogprofiler python "$SRC/scripts/converters/convert_competitor_groups.py" \
 --tool fastoma --input "$PRIMARY" --fasta "$OUT/input/proteome" --out "$POST/groups.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/converters/validate_groups.py" \
 "$POST/groups.tsv" --fasta "$OUT/input/proteome" > "$POST/validation.log"
[[ $(($(wc -l < "$POST/groups.tsv") - 1)) -eq 30 ]] || { echo "Converted smoke partition does not cover 30 proteins" >&2; exit 6; }
[[ $(cut -f1 "$POST/groups.tsv" | tail -n +2 | sort -u | wc -l) -eq 30 ]] || {
 echo "Header-only FastOMA primary must produce 30 singleton groups" >&2; exit 6;
}
printf 'primary_output\t%s\ngroups\t%s\nreport_status\tFAILED_OPTIONAL\npartition_status\tPASS\nstatus\tPASS\n' \
 "$PRIMARY" "$POST/groups.tsv" > "$POST/SMOKE_STATUS.tsv"
