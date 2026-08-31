#!/usr/bin/env bash
#SBATCH --job-name=b1-sonic-smoke-post
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -euo pipefail

SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
MM=/home/mselab/licj/micromamba/bin/micromamba
OUT=$ROOT/03_runs/B1_orthobench/competitor_smoke_stage2a/sonicparanoid-2.0.9_slurm4
PRIMARY=$OUT/result/runs/smoke/ortholog_groups/ortholog_groups.tsv
POST=$OUT/postprocess_retry1

[[ -f "$PRIMARY" && -d "$OUT/input" ]] || { echo "Completed SonicParanoid smoke artifacts are absent" >&2; exit 4; }
[[ -f "$OUT/timing/metadata.tsv" ]] || { echo "SonicParanoid timing metadata is absent" >&2; exit 4; }
grep -q $'^status\tOK$' "$OUT/timing/metadata.tsv"
[[ ! -e "$POST" ]] || { echo "Refusing existing postprocess retry directory: $POST" >&2; exit 3; }
mkdir -p "$POST"

cat > "$POST/MANIFEST.tsv" <<EOF
postprocess_of_slurm_job_id\t1404434
inference_status\tcompleted_successfully
failure\tcompact_six_column_SonicParanoid_schema_not_supported_by_converter
group_identity\tdeterministic_one_based_primary_output_row_ordinal
slurm_job_id\t${SLURM_JOB_ID:?}
git_commit\t${DEV_GIT_COMMIT:?}
source_sha256\t${DEV_SOURCE_HASH:?}
EOF
sha256sum "$PRIMARY" > "$POST/PRIMARY_SHA256SUM"
"$MM" run -n ogprofiler python "$SRC/scripts/converters/convert_competitor_groups.py" \
 --tool sonicparanoid --input "$PRIMARY" --fasta "$OUT/input" --out "$POST/groups.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/converters/validate_groups.py" \
 "$POST/groups.tsv" --fasta "$OUT/input" > "$POST/validation.log"
printf 'primary_output\t%s\ngroups\t%s\nstatus\tPASS\n' "$PRIMARY" "$POST/groups.tsv" > "$POST/SMOKE_STATUS.tsv"
