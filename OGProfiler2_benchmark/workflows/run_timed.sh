#!/usr/bin/env bash
# Run one benchmark command while preserving logs and resource/provenance records.
set -uo pipefail
if [[ $# -lt 6 || "$5" != "--" ]]; then
  echo "Usage: $0 TOOL DATASET REPLICATE OUTPUT_DIR -- COMMAND [ARG ...]" >&2
  exit 2
fi
TOOL=$1; DATASET=$2; REPLICATE=$3; OUT=$4; shift 5
mkdir -p "$OUT"
printf '%q ' "$@" > "$OUT/command.txt"; printf '\n' >> "$OUT/command.txt"
start=$(date -u +%Y-%m-%dT%H:%M:%SZ)
host=$(hostname)
source_root="$(cd "$(dirname "$0")/../.." && pwd -P)"
git_commit=${DEV_GIT_COMMIT:-${BENCHMARK_GIT_COMMIT:-}}
if [[ -z "$git_commit" ]]; then
  git_commit=$(git -C "$source_root" rev-parse HEAD 2>/dev/null || awk -F= '$1 == "git_commit" {print $2}' "$source_root/.dev-deploy-meta" 2>/dev/null || true)
fi
git_commit=${git_commit:-UNKNOWN}
cpu=$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || printf UNKNOWN)
threads=${BENCHMARK_THREADS:-UNKNOWN}; seed=${BENCHMARK_SEED:-UNKNOWN}
set +e
/usr/bin/time -v -o "$OUT/resources.time.txt" "$@" >"$OUT/stdout.log" 2>"$OUT/stderr.log"
code=$?
set -e
end=$(date -u +%Y-%m-%dT%H:%M:%SZ)
status=OK
[[ $code -ne 0 ]] && status=FAILED
grep -qiE 'out of memory|oom|killed' "$OUT/stderr.log" "$OUT/resources.time.txt" 2>/dev/null && status=OOM || true
grep -qiE 'timed out|timeout' "$OUT/stderr.log" "$OUT/resources.time.txt" 2>/dev/null && status=TIMEOUT || true
peak_kb=$(awk -F: '/Maximum resident set size/{gsub(/ /,"",$2); print $2}' "$OUT/resources.time.txt")
wall=$(sed -n 's/.*Elapsed (wall clock) time (h:mm:ss or m:ss): *//p' "$OUT/resources.time.txt")
user_cpu=$(awk -F: '/User time \(seconds\)/{gsub(/ /,"",$2); print $2}' "$OUT/resources.time.txt")
system_cpu=$(awk -F: '/System time \(seconds\)/{gsub(/ /,"",$2); print $2}' "$OUT/resources.time.txt")
disk=$(du -sk "$OUT" | awk '{print $1 * 1024}')
{
 echo -e "key\tvalue"
 echo -e "tool\t$TOOL"; echo -e "dataset\t$DATASET"; echo -e "replicate\t$REPLICATE"
 echo -e "start_time_utc\t$start"; echo -e "end_time_utc\t$end"; echo -e "hostname\t$host"
 echo -e "cpu_count\t$cpu"; echo -e "threads\t$threads"; echo -e "seed\t$seed"; echo -e "git_commit\t$git_commit"
 echo -e "exit_code\t$code"; echo -e "status\t$status"; echo -e "peak_rss_kb\t${peak_kb:-}"
 echo -e "wall_clock\t${wall:-}"; echo -e "user_cpu_seconds\t${user_cpu:-}"; echo -e "system_cpu_seconds\t${system_cpu:-}"; echo -e "disk_bytes\t$disk"
} > "$OUT/metadata.tsv"
exit "$code"
