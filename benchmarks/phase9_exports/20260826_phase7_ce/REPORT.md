# Phase 9 terminal-family final export regression

## Provenance

- Source hierarchy: Phase 7 Slurm run `1404182`
- Source network annotations: Phase 8 `network-overlap-v1`
- Source run ID: `20260826T060816Z_e46af3ed522d_55098951_14986`
- Source Git commit: `e46af3ed522d4d8b21870a23409791e4d08e1953`
- Export algorithm: `terminal-family-export-v1`
- Family-ID algorithm: `canonical-terminal-membership-rank-v1`
- Remote environment: project `ogprofiler` Conda environment

The short export regression reused completed C/E hierarchy and evolution
artifacts, so no additional Slurm workload was scheduled. Phase 9 was deployed
from a captured dirty development state; the Git commit identifies its base
revision rather than the uncommitted Phase 9 source.

## Results

| Dataset | Families | Members | Hierarchy nodes | Events | OG ID range | Resume |
|---|---:|---:|---:|---:|---|---|
| C | 8 | 80 | 9 | 9 | `OG000000000`–`OG000000007` | verified reuse |
| E | 6 | 240 | 7 | 7 | `OG000000000`–`OG000000005` | verified reuse |

Every prepared protein occurs exactly once in `members.tsv`; every referenced
family occurs in `families.tsv`; hierarchy and event row counts agree; and all
four archived output hashes match `export-manifest.json`. No FASTA files were
created by default.

The immediate second invocation preserved all result tables and returned the
verified export manifest. Automated regression also permutes terminal
membership row order and confirms byte-identical family and member tables.

## Validation

- local suite: 54 tests passed;
- remote Conda suite: 54 tests passed;
- Ruff and mypy (50 source files) passed;
- Phase 0 fixture checksum verification passed;
- selected-family FASTA, corruption reconstruction, singleton handoff, and
  opt-in component GraphML are covered by automated tests.
