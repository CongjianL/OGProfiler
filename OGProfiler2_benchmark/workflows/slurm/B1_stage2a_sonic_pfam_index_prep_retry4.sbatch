#!/usr/bin/env bash
#SBATCH --job-name=b1-sonic-pfam-index-r4
#SBATCH --cpus-per-task=56
#SBATCH --mem=200G
#SBATCH --time=24:00:00
set -euo pipefail

SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
MM=/home/mselab/licj/micromamba/bin/micromamba
ENV_NAME=sonicparanoid-2.0.9-mmseqs13
ENV=/home/mselab/licj/micromamba/envs/$ENV_NAME
PROFILE_DIR=$ENV/lib/python3.12/site-packages/sonicparanoid/pfam_files/profile_db
PROFILE=$PROFILE_DIR/pfama.mmseqs
MMSEQS=$ENV/bin/mmseqs
OUT=$ROOT/03_runs/B1_orthobench/environment_preparation/sonicparanoid-2.0.9_mmseqs13_pfam_index_unsplit5
FAILED=$OUT/successful_split_index_from_1404428

[[ ! -e "$OUT" ]] || { echo "Refusing existing environment-preparation retry directory: $OUT" >&2; exit 3; }
[[ -x "$MMSEQS" && -f "$PROFILE" && -f "$PROFILE.index" && -f "$PROFILE.dbtype" ]]
[[ "$($MMSEQS version)" == 13.45111 ]]
mkdir -p "$OUT" "$FAILED"

# Preserve the scientifically valid split index whose disk layout is rejected
# by SonicParanoid's file-existence guard.
shopt -s nullglob
old_index=("$PROFILE".idx*)
((${#old_index[@]} >= 3)) || { echo "Expected split index from 1404428 is absent" >&2; exit 4; }
mv "${old_index[@]}" "$FAILED/"
shopt -u nullglob
find "$FAILED" -maxdepth 1 -type f -printf '%s\t%f\n' | LC_ALL=C sort -k2,2 > "$OUT/preserved_split_index.tsv"
sha256sum "$FAILED"/* > "$OUT/preserved_split_index.sha256"

cat > "$OUT/PREPARATION_MANIFEST.tsv" <<EOF
artifact	SonicParanoid2_bundled_PfamA_MMseqs_index
sonicparanoid_version	2.0.9
sonicparanoid_build	py312hc9302aa_0
mmseqs_version	13.45111
mmseqs_build	pl5321hf1761c0_2
retry_of_slurm_job_id	1404428
reason	SonicParanoid_requires_physical_unsplit_idx_file
method	MMseqs2_13.45111_createindex_with_one_explicit_partition
explicit_split_count	1
kmer	5
sensitivity	7.0
threads	56
resource_change	56_CPU_200_GB_per_explicit_user_direction
slurm_job_id	${SLURM_JOB_ID:?}
git_commit	${DEV_GIT_COMMIT:?}
source_sha256	${DEV_SOURCE_HASH:?}
EOF
"$MM" list -n "$ENV_NAME" --explicit > "$OUT/conda-explicit.txt"

"$SRC/workflows/run_timed.sh" MMseqs2 SonicParanoidPfamIndex 5 "$OUT/timing" -- \
 "$MMSEQS" createindex "$PROFILE" "$OUT/mmseqs_tmp" -k 5 -s 7.0 \
 --threads 56 -v 3 --split 1

[[ -f "$PROFILE.idx" && -f "$PROFILE.idx.index" && -f "$PROFILE.idx.dbtype" ]]
[[ "$($MMSEQS dbtype "$PROFILE.idx")" == Index ]]
find "$PROFILE_DIR" -maxdepth 1 -type f -name 'pfama.mmseqs.idx*' \
 -printf '%s\t%f\n' | LC_ALL=C sort -k2,2 > "$OUT/index_files.tsv"
sha256sum "$PROFILE".idx* > "$OUT/index_files.sha256"
printf 'status\tPASS\n' > "$OUT/PREPARATION_STATUS.tsv"
