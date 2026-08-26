# OGProfiler 1 reference implementation

`OGProfiler_v1.py` is a frozen copy of the legacy implementation from commit
`e46af3e`. It exists only as a reference implementation for regression and
performance comparisons with OGProfiler 2.

## Freeze policy

- Do not add features to this file.
- Record critical defects as annotations or tests before changing behavior.
- Apply product development to `src/ogprofiler/`.
- Preserve the exact command, input manifest, dependency versions, and outputs
  for every baseline run.

## Captured development environment (2026-08-25)

- Python available on the local development machine: 3.12.6
- Python packages required by V1: `ete3`, `igraph`, `leidenalg`, `numpy`,
  `progressbar2`, and `scipy`
- NCBI BLAST+: 2.16.0+
- FastTree: 2.2.0 (double precision)
- MAFFT: 7.526
- DIAMOND: not present in the inspected local environment
- MMseqs2: not present in the inspected local environment

Package versions for a baseline run must be captured from the execution
environment rather than inferred from this development-machine inventory.

## Complete command template

```bash
python legacy/OGProfiler_v1.py \
  --in INPUT_PROTEOME_DIRECTORY \
  --out BASELINE_OUTPUT_DIRECTORY \
  --extension faa \
  --search_method diamond \
  --evalue 0.001 \
  --threads 8 \
  --distance lrb \
  --weight NBS \
  --community_method rber \
  --network_threads 8 \
  --gamma_coefficient 1.0 \
  --species_overlap 0
```

Use `--continue` only when documenting a resume test. Use `--refined` only
when the external alignment and tree-building tools are part of the recorded
baseline environment.

