# OGProfiler release checklist

## Code and compatibility

- [ ] Working tree contains only intended release changes.
- [x] `ruff check .` passes.
- [x] `mypy src/ogprofiler` passes.
- [x] Unit, integration, and regression tests pass with `pytest -q`.
- [x] Frozen Phase 0 datasets pass `generate_phase0_datasets.py --check`.
- [x] Remote smoke test passes in the supported `ogprofiler` environment.

## Scientific and performance evidence

- [x] Phase 12 scientific parameter benchmark archived.
- [x] Phase 13 synthetic evolution benchmark archived.
- [x] Phase 14 scheduler, resume, memory, and I/O benchmark archived.
- [ ] Any algorithm change since those runs has a new targeted regression.

## Documentation and artifacts

- [x] README includes concept, installation, quick start, outputs, and caveats.
- [x] Reproducibility guide covers seeds, versions, manifests, `run.yaml`, and tools.
- [x] CLI help and examples match the release build.
- [x] Changelog contains the release version and date.
- [x] Wheel and source distribution build successfully.
- [x] Installed wheel passes `ogprofiler --version` and CLI smoke tests.

## Release publication

- [x] Package versions in `pyproject.toml` and `ogprofiler/__init__.py` agree.
- [x] Release commit is clean and pushed.
- [x] Annotated version tag is created from the release commit and pushed.
- [ ] Release notes link benchmark and reproducibility evidence.
