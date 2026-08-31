#!/usr/bin/env bash
#SBATCH --job-name=b1-sonic-pfam-index-r2
#SBATCH --cpus-per-task=56
#SBATCH --mem=200G
#SBATCH --time=24:00:00
set -euo pipefail

SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
ENV=/home/mselab/licj/micromamba/envs/sonicparanoid-2.0.9
PROFILE=$ENV/lib/python3.12/site-packages/sonicparanoid/pfam_files/profile_db/pfama.mmseqs
MMSEQS=$ENV/bin/mmseqs
OUT=$ROOT/03_runs/B1_orthobench/environment_preparation/sonicparanoid-2.0.9_pfam_index_split3
FAILED=$OUT/incomplete_index_from_1404342

[[ ! -e "$OUT" ]] || { echo "Refusing existing environment-preparation retry directory: $OUT" >&2; exit 3; }
[[ -x "$MMSEQS" && -f "$PROFILE" && -f "$PROFILE.index" && -f "$PROFILE.dbtype" ]]
mkdir -p "$OUT" "$FAILED"

shopt -s nullglob
old_index=("$PROFILE".idx*)
if ((${#old_index[@]})); then mv "${old_index[@]}" "$FAILED/"; fi
shopt -u nullglob

cat > "$OUT/PREPARATION_MANIFEST.tsv" <<EOF
artifact	SonicParanoid2_bundled_PfamA_MMseqs_index
tool_version	2.0.9
mmseqs_version	$($MMSEQS version)
retry_of_slurm_job_id	1404342
reason	automatic_split_estimator_did_not_partition_profile_index
method	MMseqs_createindex_with_eight_explicit_equal_partitions
explicit_split_count	8
kmer	5
sensitivity	7.0
threads	56
resource_change	56_CPU_200_GB_per_explicit_user_direction
slurm_job_id	${SLURM_JOB_ID:?}
git_commit	${DEV_GIT_COMMIT:?}
source_sha256	${DEV_SOURCE_HASH:?}
EOF

"$SRC/workflows/run_timed.sh" MMseqs2 SonicParanoidPfamIndex 3 "$OUT/timing" -- \
 "$MMSEQS" createindex "$PROFILE" "$OUT/mmseqs_tmp" -k 5 -s 7.0 \
 --threads 56 -v 3 --split 8

[[ -f "$PROFILE.idx" && -f "$PROFILE.idx.index" && -f "$PROFILE.idx.dbtype" ]]
find "$ENV/lib/python3.12/site-packages/sonicparanoid/pfam_files/profile_db" \
 -maxdepth 1 -type f -name 'pfama.mmseqs.idx*' -printf '%s\t%f\n' \
 | LC_ALL=C sort -k2,2 > "$OUT/index_files.tsv"
sha256sum "$PROFILE".idx* > "$OUT/index_files.sha256"
printf 'status\tPASS\n' > "$OUT/PREPARATION_STATUS.tsv"
