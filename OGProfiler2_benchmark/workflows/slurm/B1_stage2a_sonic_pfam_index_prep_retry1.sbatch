#!/usr/bin/env bash
#SBATCH --job-name=b1-sonic-pfam-index-r1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
set -euo pipefail

SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
ENV=/home/mselab/licj/micromamba/envs/sonicparanoid-2.0.9
PROFILE=$ENV/lib/python3.12/site-packages/sonicparanoid/pfam_files/profile_db/pfama.mmseqs
MMSEQS=$ENV/bin/mmseqs
OUT=$ROOT/03_runs/B1_orthobench/environment_preparation/sonicparanoid-2.0.9_pfam_index_split2
FAILED=$OUT/incomplete_index_from_1404331

[[ ! -e "$OUT" ]] || { echo "Refusing existing environment-preparation retry directory: $OUT" >&2; exit 3; }
[[ -x "$MMSEQS" && -f "$PROFILE" && -f "$PROFILE.index" && -f "$PROFILE.dbtype" ]]
mkdir -p "$OUT" "$FAILED"

# Preserve the 80 GB-per-split attempt before recreating the index.
shopt -s nullglob
old_index=("$PROFILE".idx*)
if ((${#old_index[@]})); then mv "${old_index[@]}" "$FAILED/"; fi
shopt -u nullglob

cat > "$OUT/PREPARATION_MANIFEST.tsv" <<EOF
artifact	SonicParanoid2_bundled_PfamA_MMseqs_index
tool_version	2.0.9
mmseqs_version	$($MMSEQS version)
retry_of_slurm_job_id	1404331
reason	80_GB_per_split_left_insufficient_index_construction_headroom
method	MMseqs_createindex_with_automatic_splitting
split_memory_limit	40G
threads	32
slurm_job_id	${SLURM_JOB_ID:?}
git_commit	${DEV_GIT_COMMIT:?}
source_sha256	${DEV_SOURCE_HASH:?}
EOF

"$SRC/workflows/run_timed.sh" MMseqs2 SonicParanoidPfamIndex 2 "$OUT/timing" -- \
 "$MMSEQS" createindex "$PROFILE" "$OUT/mmseqs_tmp" -k 5 -s 7.0 \
 --threads 32 -v 3 --split 0 --split-memory-limit 40G

[[ -f "$PROFILE.idx" && -f "$PROFILE.idx.index" && -f "$PROFILE.idx.dbtype" ]]
find "$ENV/lib/python3.12/site-packages/sonicparanoid/pfam_files/profile_db" \
 -maxdepth 1 -type f -name 'pfama.mmseqs.idx*' -printf '%s\t%f\n' \
 | LC_ALL=C sort -k2,2 > "$OUT/index_files.tsv"
sha256sum "$PROFILE".idx* > "$OUT/index_files.sha256"
printf 'status\tPASS\n' > "$OUT/PREPARATION_STATUS.tsv"
