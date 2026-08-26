# Phase 7 scheduler and resume benchmark

## Provenance

- Slurm array job: `1404182` (`0=C`, `1=E`)
- Run ID: `20260826T060816Z_e46af3ed522d_55098951_14986`
- Immutable source SHA-256: `550989519c3a648d7fdcb9a47d4ba09485121df6e848e2afc0939af7abde8f45`
- Git commit: `e46af3ed522d4d8b21870a23409791e4d08e1953` (`dirty=1`, captured snapshot)
- Resources per task: 4 CPUs, 16 GiB, 30-minute limit
- Both array tasks: `COMPLETED`, exit `0:0`, elapsed 8 seconds
- Scheduler workers: 4; hierarchy stability mode: robust

## Results

| Dataset | Components scheduled | Start order | Terminal members | Hierarchy nodes | Initial wall s | Resume wall s | Initial RSS MiB | Resume RSS MiB |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| C | 6 | 0,1,2,3,4,5 | 80 | 9 | 2.91 | 0.76 | 74.7 | 69.9 |
| E | 1 | 0 | 240 | 7 | 2.29 | 0.72 | 76.4 | 69.8 |

C component 0 contains 60 proteins; components 1–5 contain four proteins each.
The recorded dispatch order is therefore largest-first, followed by the stable
component-ID tie-break. All seven C/E tasks finished `DONE` with one attempt.

The immediate resume invocation scheduled and completed zero tasks and skipped
all six C or one E verified artifacts. Resume wall time includes Python startup,
configuration, SQLite inspection, and output checksum verification.

## Resume and fault-isolation evidence

The automated regression suite additionally covers:

- actual subprocess termination after a task reaches `RUNNING`, followed by
  stale-task recovery;
- a worker exception caused by a damaged component partition while another
  component completes successfully;
- `--failed-only` completion after the damaged input is restored;
- partial hierarchy-output corruption and component-local reconstruction;
- hierarchy seed/config changes producing audited `INVALID` transitions;
- migration of pre-Phase-7 `run.db` task schemas;
- atomic artifact staging and manifest-last publication.

The archived `run.db`, initial manifest, resume manifest, and summary were
cross-checked locally. Task counts, states, attempts, start order, hierarchy
member counts, and skip counts agree across all artifacts.
