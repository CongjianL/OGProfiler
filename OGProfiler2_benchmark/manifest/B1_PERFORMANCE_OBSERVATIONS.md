# B1 performance observations

Formal run: `B1_orthobench_OGProfiler2_seed42_rep2`  
Slurm job: `1404237`

| Quantity | Value |
|---|---:|
| Wall time | 14,189 s |
| User CPU time | 18,030.04 s |
| System CPU time | 830.23 s |
| Total CPU time | 18,860.27 s |
| Allocated CPUs | 32 |
| Effective mean CPU cores (`CPU / wall`) | 1.3292 |
| Approximate allocation efficiency (`1.3292 / 32`) | 4.1538% |
| Peak RSS | 72,819,420 KiB (69.446 GiB) |

GNU `/usr/bin/time -v` measured the wrapped OGProfiler command. Its user/system values
normally accumulate CPU for the measured process and waited-for descendants, including
the pipeline child processes that terminate under that command, while multithreaded CPU
time is summed across threads. They do not constitute per-stage profiling and may omit
detached work or work not waited for by the measured process. Consequently, the ratio is
an allocation-level observation, not a diagnosis of which stage underused CPUs.

This observation is deferred to B4 scalability. Stage 1C makes no parallelism change.
