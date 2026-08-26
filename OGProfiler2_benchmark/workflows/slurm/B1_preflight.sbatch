#!/usr/bin/env bash
#SBATCH --job-name=ogp2-ob-preflight
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
set -euo pipefail
# Slurm compute nodes start with a minimal PATH; re-exec in the configured immutable environment.
if [[ "${OGP2_ENV_READY:-0}" != 1 ]]; then
  exec "${DEV_CONDA:?}" run -n "${DEV_CONDA_ENV:?}" env OGP2_ENV_READY=1 bash "$0"
fi
B=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
SOURCE_ROOT=${DEV_SOURCE_DIR:?}
mkdir -p "$B/00_env" "$B/manifest"
{
 echo -e 'software\tversion\texecutable\tsource'
 echo -e "OGProfiler\t$(PYTHONPATH="$SOURCE_ROOT/src" python -c 'import ogprofiler; print(ogprofiler.__version__)')\t$SOURCE_ROOT/src/ogprofiler\timport ogprofiler"
 echo -e "Python\t$(python --version 2>&1 | awk '{print $2}')\t$(command -v python)\tpython --version"
 echo -e "DIAMOND\t$(diamond version | awk '{print $3}')\t$(command -v diamond)\tdiamond version"
 PYTHONPATH="$SOURCE_ROOT/src" python - <<'PY'
import importlib.metadata as m
for name in ('igraph','leidenalg','numpy','pyarrow'):
 print(f"{name}\t{m.version(name)}\tpython-package\timportlib.metadata")
PY
 echo -e "Slurm\t$(scontrol --version | awk '{print $2}')\t$(command -v scontrol)\tscontrol --version"
 echo -e "/usr/bin/time\t$(/usr/bin/time --version 2>&1 | head -1)\t/usr/bin/time\t/usr/bin/time --version"
} > "$B/manifest/software_versions.tsv"
{
 echo -e 'key\tvalue'
 echo -e "hostname\t$(hostname)"; echo -e "cluster\t${SLURM_CLUSTER_NAME:-unknown}"; echo -e "slurm_partition\t${SLURM_JOB_PARTITION:-unknown}"
 echo -e "CPU_model\t$(lscpu | awk -F: '/Model name/{sub(/^ +/,"",$2);print $2}')"; echo -e "physical_cores\t$(lscpu -p=Core | grep -v '^#' | sort -u | wc -l)"; echo -e "logical_CPUs\t$(nproc)"
 echo -e "RAM_total_bytes\t$(awk '/MemTotal/{print $2*1024}' /proc/meminfo)"; echo -e "OS\t$(. /etc/os-release; echo "$PRETTY_NAME")"; echo -e "kernel\t$(uname -r)"
 echo -e "filesystem_type_of_run_directory\t$(stat -f -c %T "$B")"; echo -e "storage_mount\t$(df -P "$B" | awk 'NR==2{print $6}')"
 echo -e "slurm_job_cpus\t${SLURM_CPUS_PER_TASK:-}"; echo -e "slurm_job_memory\t${SLURM_MEM_PER_NODE:-}"; echo -e "slurm_time_limit\t${SLURM_TIMELIMIT:-}"
} > "$B/manifest/hardware.tsv"
TMP=${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/ogp2_diamond_preflight_$SLURM_JOB_ID; mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT
printf '>a\nMKTIIALSYIFCLVFADYKDDDDK\n>b\nMKTIIALSYIFCLVFADYKDDDE\n' > "$TMP/tiny.fa"
diamond makedb --in "$TMP/tiny.fa" --db "$TMP/tiny" --threads 1
diamond blastp --query "$TMP/tiny.fa" --db "$TMP/tiny" --out "$TMP/out.tsv" --outfmt 6 --threads 1
[[ -s "$TMP/out.tsv" ]]
printf 'DIAMOND_PREFLIGHT=PASS\nPATH=%s\nTMP=%s\n' "$PATH" "$TMP" > "$B/00_env/diamond_preflight_${SLURM_JOB_ID}.txt"
