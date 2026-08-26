#!/usr/bin/env bash
set -euo pipefail
OUT=$1; EXPECTED_COMMIT=$2; EXPECTED_SOURCE_SHA256=$3; ROOT=${4:-$PWD}
mkdir -p "$(dirname "$OUT")"
py=$(command -v python); og=$(command -v ogprofiler || true)
deploy_commit=$(awk -F= '$1=="git_commit" {print $2}' "$ROOT/.dev-deploy-meta" 2>/dev/null || true)
source_digest=$(cd "$ROOT" && find src/ogprofiler -type f -name '*.py' -print0 | LC_ALL=C sort -z | xargs -0 shasum -a 256 | shasum -a 256 | awk '{print $1}')
version=$(grep -A3 '^\[project\]' "$ROOT/pyproject.toml" | awk -F'"' '/version/{print $2}')
{
 echo '# Remote execution audit'; echo; echo "- hostname: $(hostname)"; echo "- repository: $ROOT"; echo "- deployed source commit: ${deploy_commit:-UNKNOWN}"; echo "- algorithm baseline commit: $EXPECTED_COMMIT"; echo "- source-tree SHA256: $source_digest"; echo "- pyproject version: $version"; echo "- python: $py ($($py --version))"; echo "- ogprofiler: ${og:-NOT_ON_PATH}"; [[ -n "$og" ]] && "$og" --version || true
 PYTHONPATH="$ROOT/src" "$py" - <<'PY'
import ogprofiler
print(f'- package path: {ogprofiler.__file__}')
print(f'- package version: {ogprofiler.__version__}')
PY
 echo; [[ "$source_digest" == "$EXPECTED_SOURCE_SHA256" ]] && echo 'REMOTE_SOURCE_MATCH=YES' || echo 'REMOTE_SOURCE_MATCH=NO'
} > "$OUT"
[[ "$source_digest" == "$EXPECTED_SOURCE_SHA256" ]]
