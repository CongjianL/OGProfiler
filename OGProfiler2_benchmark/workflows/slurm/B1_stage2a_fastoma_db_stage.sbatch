#!/usr/bin/env bash
#SBATCH --job-name=b1-fastoma-db
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
set -euo pipefail
SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
ENVROOT=$ROOT/00_env/fastoma-0.5.1
DB=$ENVROOT/LUCA.h5
PARTS=$ENVROOT/LUCA.h5.parts
URL='https://zenodo.org/records/20814376/files/LUCA.h5?download=1'
EXPECTED_SIZE=9918199250
EXPECTED_MD5=87a322a0a6381a8de0b539d62e277054
N=8
[[ ! -e "$DB" ]] || { echo "Refusing existing database: $DB" >&2; exit 3; }
mkdir -p "$PARTS"
mapfile -t IPS < <(grep -Ev '^[[:space:]]*(#|$)' "$SRC/02_configs/fastoma/ZENODO_ENDPOINT_IPS.txt")
[[ ${#IPS[@]} -gt 0 ]]
IP_COUNT=${#IPS[@]}; if (( IP_COUNT > 3 )); then IP_COUNT=3; fi
chunk=$(( (EXPECTED_SIZE + N - 1) / N ))
pids=()
for ((i=0;i<N;i++)); do
 start=$((i*chunk)); end=$((start+chunk-1)); if (( end >= EXPECTED_SIZE )); then end=$((EXPECTED_SIZE-1)); fi
 expected=$((end-start+1)); printf -v tag '%02d' "$i"
 done_file=$PARTS/part-$tag.done; partial=$PARTS/download-$tag.partial; ip=${IPS[$((i % IP_COUNT))]}
 if [[ -f "$done_file" && $(stat -c %s "$done_file") -eq $expected ]]; then continue; fi
 (
  curl --resolve "zenodo.org:443:$ip" -L --fail --retry 8 --retry-all-errors --connect-timeout 20 \
   --range "$start-$end" -o "$partial" "$URL"
  actual=$(stat -c %s "$partial"); [[ "$actual" -eq "$expected" ]] || { echo "part $tag size $actual != $expected" >&2; exit 5; }
  mv "$partial" "$done_file"
 ) & pids+=("$!")
done
status=0
for pid in "${pids[@]}"; do wait "$pid" || status=1; done
[[ "$status" -eq 0 ]] || exit 6
ASSEMBLY=$ENVROOT/LUCA.h5.assembling
: > "$ASSEMBLY"
for ((i=0;i<N;i++)); do printf -v tag '%02d' "$i"; cat "$PARTS/part-$tag.done" >> "$ASSEMBLY"; done
[[ $(stat -c %s "$ASSEMBLY") -eq $EXPECTED_SIZE ]]
echo "$EXPECTED_MD5  $ASSEMBLY" | md5sum -c -
mv "$ASSEMBLY" "$DB"
md5sum "$DB" > "$ENVROOT/LUCA.h5.md5"
sha256sum "$DB" > "$ENVROOT/LUCA.h5.sha256"
cat > "$ENVROOT/LUCA.h5.provenance.tsv" <<EOF
source_url\t$URL
zenodo_record\t20814376
published_size_bytes\t$EXPECTED_SIZE
published_md5\t$EXPECTED_MD5
staging_slurm_job_id\t${SLURM_JOB_ID:?}
git_commit\t${DEV_GIT_COMMIT:?}
source_sha256\t${DEV_SOURCE_HASH:?}
EOF
