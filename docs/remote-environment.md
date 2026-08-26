# Remote OGProfiler environment

The remote development and Slurm environment is declared by `environment.yml`
and resolved as the micromamba environment named `ogprofiler`.

## Provisioning

```bash
/home/mselab/licj/micromamba/bin/micromamba create -y -f environment.yml
```

Project wrappers resolve it through `REMOTE_CONDA_ENV="ogprofiler"` in the
untracked/deployment-specific `dev/remote.env` configuration.

## Validation witness (2026-08-25)

Validated on `campus-server`:

- Python 3.12.14
- PyArrow, NumPy, SciPy, igraph, leidenalg, and ETE imports
- DIAMOND 2.2.5
- MMseqs2 18.8cc5c
- BLAST+ 2.17.0
- FastTree 2.2.0
- MAFFT 7.526
- complete project test suite

Recreate from `environment.yml` after dependency changes, then run
`./dev/remote-test` to verify the documented project invocation.

