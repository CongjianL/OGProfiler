# Reproducibility guide

## Source and environment

Use a clean, tagged source revision for formal analysis. Record `git rev-parse
HEAD`, the tag, dirty status, and any patch. Build the environment from the
version-controlled `environment.yml`; preserve an explicit package export when
long-term recreation matters. OGProfiler's package and algorithm versions are
recorded independently because a package release may contain several stage
algorithms.

## Seed

`hierarchy.seed` controls deterministic Leiden calls and defaults to `42`.
Scientific comparisons must use the same seed policy. Robust/publication modes
may evaluate additional deterministic seeds; these parameters are included in
hierarchy manifests. Synthetic benchmark scenarios also record their own seed.

## Resolved configuration

`prepare` writes the complete resolved configuration to `run.yaml`. Later stages
load it automatically unless another `--config` is supplied. Preserve this file
rather than reconstructing settings from shell history. `--set` overrides are
part of each stage's recorded command and parameter identity.

## Manifests and checksums

Keep the complete run workspace. Important manifests include:

- `manifest.json`: prepared dataset identity, counts, seed, and config checksum;
- stage manifests under `search/`, `edges/`, `components/`, `hierarchy/`, and
  `evolution/`;
- `results/export-manifest.json`: final table checksums and counts;
- `provenance.json`: runtime, platform, configuration, and external versions.

A resume is accepted only when the relevant inputs, parameters, algorithm
version, and output checksums match. Editing an artifact makes that stage
ineligible for verified reuse.

## External tools

`provenance.json` probes the configured search, alignment, and tree executables.
For the default workflow this includes DIAMOND; optional phylogenetic refinement
also uses MAFFT and FastTree. Preserve executable paths, version output, and the
Conda/container identity. A missing optional tool is recorded as unavailable
and does not affect stages that do not use it.

## Formal execution record

For a local formal run, archive:

```text
git commit/tag and dirty status
environment lock/export
complete run workspace
stdout/stderr and ogprofiler.log
input manifest and checksums
run.yaml and provenance.json
random seed policy
```

For Slurm, additionally archive job/array IDs, immutable source hash, submission
script and arguments, requested resources, `sacct`, stdout/stderr, run directory,
and retry lineage. Queue completion and scientific acceptance are separate
claims.
