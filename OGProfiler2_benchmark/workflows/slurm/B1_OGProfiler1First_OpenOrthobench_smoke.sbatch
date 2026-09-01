#!/usr/bin/env bash
#SBATCH --job-name=b1-ogp1-smoke
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
set -euo pipefail
SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
DATA=/home/mselab/licj/project_data/OGProfiler/Open_Orthobench-master/BENCHMARKS/Input
MM=/home/mselab/licj/micromamba/bin/micromamba
ENV=ogprofiler-v1-first
OUT=$ROOT/03_runs/B1_orthobench/OGProfiler1First/smoke_retry2
[[ ! -e "$OUT" ]] || { echo "Refusing existing smoke directory: $OUT" >&2; exit 3; }
mkdir -p "$OUT"
"$MM" run -n ogprofiler python "$SRC/scripts/preprocess/prepare_ogprofiler_v1_input.py" \
 --input "$DATA" --out "$OUT/input" --map "$OUT/id_map.tsv" --limit-per-file 10 --min-length 20
"$MM" list -n "$ENV" --explicit > "$OUT/environment.explicit.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BENCHMARK_THREADS=4
cd "$OUT"
"$SRC/workflows/run_timed.sh" OGProfiler1First OpenOrthobenchSmoke 1 "$OUT/raw_tool" -- \
 "$MM" run -n "$ENV" python "$SRC/00_env/ogprofiler_v1_first/OGProfiler.py" \
 -i "$OUT/input" -o "$OUT/ogprofiler_v1_first.gml" -s diamond -e 1e-5 -t 4 -d lrb -w NBS -m rber -a 4 -g 1.0 --so 0
PRIMARY="$OUT/input/WorkingDirectory/OGFile_coalescence_SameGenome.txt"
[[ -s "$PRIMARY" ]]
"$MM" run -n ogprofiler python "$SRC/scripts/converters/convert_competitor_groups.py" \
 --tool ogprofiler_v1 --input "$PRIMARY" --id-map "$OUT/id_map.tsv" --fasta "$DATA" --out "$OUT/groups.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/converters/validate_groups.py" "$OUT/groups.tsv" --fasta "$DATA" > "$OUT/validation.log"
