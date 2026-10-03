"""H5-only continuation from a verified immutable H4, without mutating its tree."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from ogprofiler.config import hierarchy_config, load_config
from ogprofiler.core.manifest import sha256_file
from ogprofiler.hierarchy.gate import require_resolved_hierarchy
from ogprofiler.hierarchy.stage import hierarchy_component_is_verified


def prepare(root: Path, baseline: Path, origin: Path):
    from benchmarks.og_extraction.orthobench import OFFICIAL_SCORER_SHA256

    acceptance = json.loads((baseline / "h4-acceptance.json").read_text())
    prefix = json.loads((baseline / "depth-prefix-acceptance.json").read_text())
    assert acceptance["passed"] and all(acceptance["checks"].values())
    assert prefix["passed"] and prefix["resolved"]
    source = baseline / "new-hierarchy"
    options = load_config(str(source / "run.yaml"))
    config = hierarchy_config(options["hierarchy"])
    assert config.max_depth == 42 and config.resolution.topology_policy == "soft_binary_24_v2"
    assert config.recursion_stop_size == 1 and config.resolution.stability_threshold == 0.9
    assert hierarchy_component_is_verified(source, 0, config)
    run = root / "new-hierarchy"
    run.mkdir(exist_ok=False)
    for name in ("input", "edges", "components"):
        (run / name).symlink_to(source / name, target_is_directory=True)
    shutil.copy2(source / "run.yaml", run / "run.yaml")
    component = "hierarchy/components/component=00000000"
    shutil.copytree(source / component, run / component)
    assert hierarchy_component_is_verified(run, 0, config)
    fixed = json.loads((baseline / "fixed-inputs.json").read_text())
    assert all(sha256_file(run / p) == h for p, h in fixed.items())
    (root / "fixed-inputs.json").write_text(json.dumps(fixed, indent=2))
    shutil.copytree(origin / "benchmark", root / "benchmark")
    assert sha256_file(root / "benchmark/benchmark.py") == OFFICIAL_SCORER_SHA256
    hashes = {
        str(p.relative_to(root / "benchmark")): sha256_file(p)
        for p in (root / "benchmark").rglob("*")
        if p.is_file()
    }
    (root / "benchmark-inputs.json").write_text(json.dumps(hashes, indent=2))
    manifest = json.loads((source / component / "hierarchy-manifest.json").read_text())
    record = dict(
        baseline=str(baseline),
        baseline_job="1410845",
        origin=str(origin),
        component0_reused=True,
        configuration_unchanged=True,
        baseline_acceptance=acceptance,
        baseline_prefix=prefix,
        component0_output_checksums=manifest["output_checksums"],
    )
    (root / "h5-reuse-provenance.json").write_text(json.dumps(record, indent=2))


def freeze(root: Path):
    run = root / "new-hierarchy"
    require_resolved_hierarchy(run)
    reused = json.loads((root / "h5-reuse-provenance.json").read_text())
    folder = run / "hierarchy/components/component=00000000"
    assert all(
        sha256_file(folder / p) == h for p, h in reused["component0_output_checksums"].items()
    )
    manifests = {
        str(p.relative_to(run)): sha256_file(p)
        for p in (run / "hierarchy/components").glob("*/hierarchy-manifest.json")
    }
    (root / "h5-hierarchy-freeze.json").write_text(
        json.dumps(
            dict(
                resolved_gate_passed=True,
                component0_unchanged=True,
                configuration_sha256=sha256_file(run / "run.yaml"),
                manifests=manifests,
            ),
            indent=2,
        )
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("command", choices=("prepare", "freeze"))
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--baseline", type=Path)
    p.add_argument("--origin", type=Path)
    args = p.parse_args()
    if args.command == "prepare":
        assert args.baseline and args.origin
        prepare(args.run, args.baseline, args.origin)
    else:
        freeze(args.run)


if __name__ == "__main__":
    main()
