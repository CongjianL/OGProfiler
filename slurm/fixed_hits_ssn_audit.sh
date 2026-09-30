#!/usr/bin/env bash
# Read existing immutable hits only; no search, Leiden, MCL, or OG changes.
# Same project scheduler settings as S0; single Python worker, conservative RAM.
#SBATCH --job-name=ogp-fixed-ssn-audit
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=12:00:00
set -euo pipefail
: "${DEV_SOURCE_DIR:?Missing source snapshot}"
: "${DEV_RUN_DIR:?Missing run directory}"
: "${DEV_CONDA:?Missing configured environment manager}"
: "${DEV_CONDA_ENV:?Missing configured execution environment}"
INPUT_RUN=${1:?Usage: fixed_hits_ssn_audit.sh INPUT_RUN OF_SOURCE_ROOT}
OF_SOURCE_ROOT=${2:?Provide the installed orthofinder source directory}
[[ -f "$INPUT_RUN/search/hits.parquet" && -f "$INPUT_RUN/input/proteins.parquet" ]]
[[ -f "$OF_SOURCE_ROOT/tools/waterfall.py" && -f "$OF_SOURCE_ROOT/orthogroups/gathering.py" ]]
if [[ ${OGP_AUDIT_ENV_READY:-0} != 1 ]]; then
    exec "$DEV_CONDA" run -n "$DEV_CONDA_ENV" env OGP_AUDIT_ENV_READY=1 bash "$0" "$@"
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONPATH="$DEV_SOURCE_DIR/src"
"$DEV_CONDA" list -n "$DEV_CONDA_ENV" > "$DEV_RUN_DIR/environment.txt"
/usr/bin/time -v -o "$DEV_RUN_DIR/audit.time" python \
    "$DEV_SOURCE_DIR/benchmarks/compare_fixed_hits_ssn.py" \
    --run "$INPUT_RUN" --of-source "$OF_SOURCE_ROOT" --out "$DEV_RUN_DIR/ssn-audit"
