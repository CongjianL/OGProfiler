#!/usr/bin/env bash
# Existing P6 resource envelope; explicit controls, no job array.
#SBATCH --job-name=ogp-h4-h5-hierarchy
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=56
#SBATCH --mem=250G
#SBATCH --time=72:00:00
set -euo pipefail
: "${DEV_SOURCE_DIR:?Missing immutable source}"
: "${DEV_RUN_DIR:?Missing independent run directory}"
: "${DEV_CONDA:?Missing configured environment manager}"
: "${DEV_CONDA_ENV:?Missing configured environment}"
ORIGIN=${1:?Usage: h4_h5_hierarchy_regression.sh P6_RUN_ROOT PRIOR_H4_RUN_ROOT [kway_v1|soft_binary_24_v2] [full|depth42-h4-only]}
BUDGET_BASELINE=${2:?Provide the previous frozen H4 run for controlled comparison}
TOPOLOGY_POLICY=${3:-kway_v1}
CAMPAIGN_MODE=${4:-full}
case "$CAMPAIGN_MODE" in
    full) ;;
    depth42-h4-only) [[ "$TOPOLOGY_POLICY" == soft_binary_24_v2 ]] || exit 2 ;;
    *) echo "Unknown campaign mode" >&2; exit 2 ;;
esac
case "$TOPOLOGY_POLICY" in kway_v1|soft_binary_24_v2) ;; *) echo "Unknown topology policy" >&2; exit 2 ;; esac
if [[ ${OGP_H4_H5_READY:-0} != 1 ]]; then
    exec "$DEV_CONDA" run -n "$DEV_CONDA_ENV" env OGP_H4_H5_READY=1 bash "$0" "$@"
fi
export PYTHONPATH="$DEV_SOURCE_DIR/src:$DEV_SOURCE_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OGP_REGRESSION_ORIGIN="$ORIGIN"
export OGP_BUDGET_BASELINE="$BUDGET_BASELINE"
export OGP_TOPOLOGY_POLICY="$TOPOLOGY_POLICY"
export OGP_CAMPAIGN_MODE="$CAMPAIGN_MODE"
python - <<'PY'
import json, os, shutil
from pathlib import Path
from ogprofiler.core.manifest import sha256_file
from benchmarks.og_extraction.orthobench import OFFICIAL_SCORER_SHA256
root, origin = Path(os.environ['DEV_RUN_DIR']), Path(os.environ['OGP_REGRESSION_ORIGIN'])
previous = origin/'repaired-mean'
completion = json.loads((origin/'p6-completion.json').read_text())
assert completion['job_id'] == '1410751' and completion['benchmark_unchanged']
assert sha256_file(previous/'input/proteins.parquet') == 'e5323530d8f6acbdf6cab152e46869712e1965ff2cdda1b69300f6f245db106f'
run = root/'new-hierarchy'
run.mkdir()
fixed = {}
for name in ('input', 'edges', 'components'):
    shutil.copytree(previous/name, run/name)
    for p in (previous/name).rglob('*'):
        if p.is_file():
            relative = str(p.relative_to(previous))
            fixed[relative] = sha256_file(p)
            assert sha256_file(run/relative) == fixed[relative]
(root/'fixed-inputs.json').write_text(json.dumps(fixed, indent=2))
shutil.copytree(origin/'benchmark', root/'benchmark')
assert sha256_file(root/'benchmark/benchmark.py') == OFFICIAL_SCORER_SHA256
benchmark = {str(p.relative_to(root/'benchmark')): sha256_file(p)
             for p in (root/'benchmark').rglob('*') if p.is_file()}
(root/'benchmark-inputs.json').write_text(json.dumps(benchmark, indent=2))
(root/'h4-h5-input-origin.json').write_text(json.dumps(dict(origin=str(origin),
    origin_job='1410751', fixed_ssn='repaired_mean', component0_proteins=69642,
    new_hits=False, new_ssn=False, seed=42, robust_seed_count=3, stability_threshold=.9,
    git_commit=os.environ['DEV_GIT_COMMIT'], dirty=os.environ['DEV_GIT_DIRTY'],
    source_sha256=os.environ['DEV_SOURCE_HASH']), indent=2))
PY
if [[ "$CAMPAIGN_MODE" == depth42-h4-only ]]; then
    python -m benchmarks.og_extraction.depth_budget_regression configure \
        --baseline "$BUDGET_BASELINE" --current "$DEV_RUN_DIR" \
        --out "$DEV_RUN_DIR/depth-only-control.json"
else
python -m benchmarks.og_extraction.hierarchy_regression migrate \
    --run "$DEV_RUN_DIR/new-hierarchy" --origin "$ORIGIN" --out "$DEV_RUN_DIR/config-migration.json" \
    --topology-policy "$TOPOLOGY_POLICY"
fi
python - <<'PY'
import json, os, yaml
from pathlib import Path
from ogprofiler.core.manifest import sha256_file
root, baseline = Path(os.environ['DEV_RUN_DIR']), Path(os.environ['OGP_BUDGET_BASELINE'])
old = yaml.safe_load((baseline/'new-hierarchy/run.yaml').read_text())
new = yaml.safe_load((root/'new-hierarchy/run.yaml').read_text())
old['hierarchy'].setdefault('topology_policy', 'kway_v1')
policy = os.environ['OGP_TOPOLOGY_POLICY']
assert new['hierarchy']['leiden_iterations'] == 10
if os.environ['OGP_CAMPAIGN_MODE'] == 'depth42-h4-only':
    assert old['hierarchy']['topology_policy'] == policy == 'soft_binary_24_v2'
    assert old['hierarchy']['max_depth'] == 20 and new['hierarchy']['max_depth'] == 42
    changed, previous, current = 'hierarchy.max_depth', 20, 42
    old['hierarchy']['max_depth'] = 42
elif policy == 'soft_binary_24_v2':
    assert old['hierarchy']['leiden_iterations'] == 10
    assert old['hierarchy']['topology_policy'] == 'kway_v1'
    changed, previous, current = 'hierarchy.topology_policy', 'kway_v1', policy
    old['hierarchy']['topology_policy'] = policy
else:
    assert old['hierarchy']['leiden_iterations'] == 2
    changed, previous, current = 'hierarchy.leiden_iterations', 2, 10
    old['hierarchy']['leiden_iterations'] = 10
assert old == new, 'Unexpected parameter change beyond explicit controlled parameter'
fixed = json.loads((root/'fixed-inputs.json').read_text())
assert all(sha256_file(baseline/'new-hierarchy'/p) == h for p, h in fixed.items())
name = ('depth-parameter-control.json' if os.environ['OGP_CAMPAIGN_MODE'] == 'depth42-h4-only'
        else 'topology-only-control.json' if policy == 'soft_binary_24_v2'
        else 'iteration-only-control.json')
(root/name).write_text(json.dumps(dict(
    baseline=str(baseline), fixed_inputs_equal=True,
    changed_parameter=changed, previous=previous, current=current), indent=2))
PY
"$DEV_CONDA" list -n "$DEV_CONDA_ENV" > "$DEV_RUN_DIR/environment.txt"
echo '==> H4: fixed component 0, root and actual descendants'
set +e
/usr/bin/time -v -o "$DEV_RUN_DIR/h4.time" \
    python -m ogprofiler hierarchy --run "$DEV_RUN_DIR/new-hierarchy" --component-id 0 \
    > "$DEV_RUN_DIR/h4.stdout" 2> "$DEV_RUN_DIR/h4.stderr"
H4_EXIT=$?
set -e
printf '%s\n' "$H4_EXIT" > "$DEV_RUN_DIR/h4-exit-code.txt"
python -m benchmarks.og_extraction.hierarchy_regression summary \
    --run "$DEV_RUN_DIR/new-hierarchy" --component 0 --out "$DEV_RUN_DIR/h4-summary.json"
if [[ "$CAMPAIGN_MODE" == depth42-h4-only ]]; then
    python -m benchmarks.og_extraction.depth_budget_regression validate \
        --baseline "$BUDGET_BASELINE" --current "$DEV_RUN_DIR" \
        --out "$DEV_RUN_DIR/depth-prefix-acceptance.json"
fi
if [[ "$H4_EXIT" != 0 ]]; then
    printf '{"h4_exit":%s,"h5_status":"BLOCKED_BY_H4","scoring_started":false}\n' "$H4_EXIT" > "$DEV_RUN_DIR/h4-h5-completion.json"
    exit "$H4_EXIT"
fi
echo '==> H4 acceptance: independent two-worker replay of complete component 0'
python - <<'PY'
import os, shutil
from pathlib import Path
root = Path(os.environ['DEV_RUN_DIR'])
parallel = root/'h4-parallel'
parallel.mkdir()
for name in ('input', 'edges', 'components'):
    (parallel/name).symlink_to(root/'new-hierarchy'/name, target_is_directory=True)
shutil.copy2(root/'new-hierarchy/run.yaml', parallel/'run.yaml')
PY
/usr/bin/time -v -o "$DEV_RUN_DIR/h4-parallel.time" \
    python -m ogprofiler hierarchy --run "$DEV_RUN_DIR/h4-parallel" --component-id 0 \
    --set hierarchy.subtree_workers=2 \
    > "$DEV_RUN_DIR/h4-parallel.stdout" 2> "$DEV_RUN_DIR/h4-parallel.stderr"
python -m benchmarks.og_extraction.hierarchy_regression validate-h4 \
    --run "$DEV_RUN_DIR/new-hierarchy" --parallel "$DEV_RUN_DIR/h4-parallel" \
    --out "$DEV_RUN_DIR/h4-acceptance.json"
if [[ "$CAMPAIGN_MODE" == depth42-h4-only ]]; then
    printf '{"h4_completed":true,"h4_accepted":true,"h5_status":"NOT_REQUESTED_H4_ONLY","scoring_started":false}\n' > "$DEV_RUN_DIR/h4-h5-completion.json"
    exit 0
fi
echo '==> H5: all remaining components, same policy/SSN; reuse verified component 0'
set +e
/usr/bin/time -v -o "$DEV_RUN_DIR/h5-hierarchy.time" \
    python -m ogprofiler hierarchy-all --run "$DEV_RUN_DIR/new-hierarchy" \
    > "$DEV_RUN_DIR/h5-hierarchy.stdout" 2> "$DEV_RUN_DIR/h5-hierarchy.stderr"
H5_EXIT=$?
set -e
printf '%s\n' "$H5_EXIT" > "$DEV_RUN_DIR/h5-hierarchy-exit-code.txt"
python -m benchmarks.og_extraction.hierarchy_regression summary \
    --run "$DEV_RUN_DIR/new-hierarchy" --out "$DEV_RUN_DIR/h5-hierarchy-summary.json"
if [[ "$H5_EXIT" != 0 ]]; then
    printf '{"h5_hierarchy_exit":%s,"h5_status":"BLOCKED_BY_HIERARCHY","scoring_started":false}\n' "$H5_EXIT" > "$DEV_RUN_DIR/h4-h5-completion.json"
    exit "$H5_EXIT"
fi
python -m ogprofiler annotate-network --run "$DEV_RUN_DIR/new-hierarchy" \
    > "$DEV_RUN_DIR/h5-events.stdout" 2> "$DEV_RUN_DIR/h5-events.stderr"
/usr/bin/time -v -o "$DEV_RUN_DIR/h5-parity.time" \
    python -m benchmarks.og_extraction.fixed_hierarchy --run "$DEV_RUN_DIR/new-hierarchy" \
    --out "$DEV_RUN_DIR/h5-parity" > "$DEV_RUN_DIR/h5-parity.stdout" 2> "$DEV_RUN_DIR/h5-parity.stderr"
python -m ogprofiler export --run "$DEV_RUN_DIR/new-hierarchy" --strategy terminal \
    > "$DEV_RUN_DIR/h5-terminal.stdout" 2> "$DEV_RUN_DIR/h5-terminal.stderr"
/usr/bin/time -v -o "$DEV_RUN_DIR/h5-scoring.time" \
    python -m benchmarks.og_extraction.hierarchy_regression score --run "$DEV_RUN_DIR" --origin "$ORIGIN" \
    > "$DEV_RUN_DIR/h5-scoring.stdout" 2> "$DEV_RUN_DIR/h5-scoring.stderr"
printf '{"h4_completed":true,"h5_status":"EVALUATED","accuracy_improvement_claimed":false}\n' > "$DEV_RUN_DIR/h4-h5-completion.json"
