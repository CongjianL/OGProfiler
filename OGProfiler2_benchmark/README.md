# OGProfiler 2 benchmark engineering

This repository-contained scaffold is deployed with the source. On the compute host, run the audit into `/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark`; raw inputs remain symlink-only and read-only. Full QfO is intentionally blocked by `manifest/QFO_READINESS.md`.

`workflows/run_ogprofiler.sh INPUT OUTPUT 32 42` is the canonical OGProfiler invocation. Run Open Orthobench only after the small subset smoke test and with the official `BENCHMARKS/benchmark.py` interface.
