#!/usr/bin/env bash
#SBATCH --job-name=b1-sonic-pfam-index-r3
#SBATCH --cpus-per-task=56
#SBATCH --mem=200G
#SBATCH --time=24:00:00
set -euo pipefail

SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
MM=/home/mselab/licj/micromamba/bin/micromamba
ENV_NAME=sonicparanoid-2.0.9-mmseqs13
ENV=/home/mselab/licj/micromamba/envs/$ENV_NAME
PFAM_FILES=$ENV/lib/python3.12/site-packages/sonicparanoid/pfam_files
ARCHIVE=$PFAM_FILES/pfama-mmseqs.tar.gz
PROFILE_DIR=$PFAM_FILES/profile_db
PROFILE=$PROFILE_DIR/pfama.mmseqs
MMSEQS=$ENV/bin/mmseqs
OUT=$ROOT/03_runs/B1_orthobench/environment_preparation/sonicparanoid-2.0.9_mmseqs13_pfam_index_split4

[[ ! -e "$OUT" ]] || { echo "Refusing existing environment-preparation retry directory: $OUT" >&2; exit 3; }
[[ -x "$MMSEQS" && -f "$ARCHIVE" ]]
[[ "$($MMSEQS version)" == 13.45111 ]]
[[ "$($MM run -n "$ENV_NAME" python -c 'import importlib.metadata as m; print(m.version("sonicparanoid"))')" == 2.0.9 ]]
mkdir -p "$OUT"

# This is a new isolated environment. Extract the bundled immutable Pfam archive
# exactly as SonicParanoid does before creating its reusable index.
[[ ! -e "$PROFILE_DIR" ]] || { echo "Refusing pre-existing profile directory in pinned environment" >&2; exit 4; }
mkdir -p "$PROFILE_DIR"
tar -xzf "$ARCHIVE" -C "$PROFILE_DIR"
for suffix in '' .index .dbtype .version _h _h.index _h.dbtype; do
 [[ -f "$PROFILE$suffix" ]] || { echo "Missing extracted Pfam file: $PROFILE$suffix" >&2; exit 5; }
done

cat > "$OUT/PREPARATION_MANIFEST.tsv" <<EOF
artifact	SonicParanoid2_bundled_PfamA_MMseqs_index
sonicparanoid_version	2.0.9
sonicparanoid_build	py312hc9302aa_0
mmseqs_version	13.45111
mmseqs_build	pl5321hf1761c0_2
mmseqs_sha256	$(sha256sum "$MMSEQS" | awk '{print $1}')
retry_of_slurm_job_id	1404345
reason	MMseqs2_18.8cc5c_segfault_during_split_profile_index
method	MMseqs2_13.45111_createindex_with_eight_explicit_equal_partitions
explicit_split_count	8
kmer	5
sensitivity	7.0
threads	56
resource_change	56_CPU_200_GB_per_explicit_user_direction
slurm_job_id	${SLURM_JOB_ID:?}
git_commit	${DEV_GIT_COMMIT:?}
source_sha256	${DEV_SOURCE_HASH:?}
EOF
"$MM" list -n "$ENV_NAME" --explicit > "$OUT/conda-explicit.txt"
sha256sum "$ARCHIVE" "$PROFILE" "$PROFILE.index" "$PROFILE.dbtype" > "$OUT/input_files.sha256"

"$SRC/workflows/run_timed.sh" MMseqs2 SonicParanoidPfamIndex 4 "$OUT/timing" -- \
 "$MMSEQS" createindex "$PROFILE" "$OUT/mmseqs_tmp" -k 5 -s 7.0 \
 --threads 56 -v 3 --split 8

[[ -f "$PROFILE.idx" && -f "$PROFILE.idx.index" && -f "$PROFILE.idx.dbtype" ]]
find "$PROFILE_DIR" -maxdepth 1 -type f -name 'pfama.mmseqs.idx*' \
 -printf '%s\t%f\n' | LC_ALL=C sort -k2,2 > "$OUT/index_files.tsv"
sha256sum "$PROFILE".idx* > "$OUT/index_files.sha256"
printf 'status\tPASS\n' > "$OUT/PREPARATION_STATUS.tsv"
