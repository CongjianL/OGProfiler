#!/usr/bin/env bash
#SBATCH --job-name=ogp-v1-mean-targets
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
# Four full small components (134 proteins), plus two isolated proteins.
# New diagnostic envelope, not a change to H5 resources or scientific budgets.
set -euo pipefail
: "${DEV_SOURCE_DIR:?Missing immutable source}"
: "${DEV_RUN_DIR:?Missing independent run directory}"
: "${DEV_CONDA:?Missing configured environment manager}"
: "${DEV_CONDA_ENV:?Missing configured environment}"
ORIGIN=${1:?Provide frozen H5 job1410868 root}
if [[ ${OGP_V1_TARGET_READY:-0} != 1 ]]; then
    exec "$DEV_CONDA" run -n "$DEV_CONDA_ENV" env OGP_V1_TARGET_READY=1 bash "$0" "$@"
fi
export PYTHONPATH="$DEV_SOURCE_DIR/src:$DEV_SOURCE_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONHASHSEED=0
[[ -f "$ORIGIN/h5-hierarchy-freeze.json" && -f "$ORIGIN/h5-metrics/report.json" ]]
python - "$ORIGIN" <<'PY'
import json, sys
from pathlib import Path
root = Path(sys.argv[1])
f = json.loads((root/'h5-hierarchy-freeze.json').read_text())
assert f['resolved_gate_passed'] and f['component0_unchanged']
assert (root/'provenance/job_id.txt').read_text().strip() == '1410868'
PY
"$DEV_CONDA" list -n "$DEV_CONDA_ENV" > "$DEV_RUN_DIR/environment.txt"
/usr/bin/time -v -o "$DEV_RUN_DIR/v1-original.time" \
    python -m benchmarks.og_extraction.v1_original_hierarchy \
    --origin "$ORIGIN" --out "$DEV_RUN_DIR/v1-original" \
    > "$DEV_RUN_DIR/v1-original.stdout" 2> "$DEV_RUN_DIR/v1-original.stderr"
