"""Fresh soft42 hierarchy on frozen Open Orthobench SSN; no parameter fitting."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import pyarrow.parquet as pq
import yaml

from benchmarks.og_extraction.orthobench import OFFICIAL_SCORER_SHA256
from ogprofiler.config import load_config, validate_config
from ogprofiler.core.manifest import sha256_file
from ogprofiler.hierarchy.gate import require_resolved_hierarchy


def execution_config(preset: Path, workers: int):
    if workers < 1:
        raise ValueError("workers must be positive")
    config = load_config(str(preset))
    if (config["hierarchy"]["topology_policy"], config["hierarchy"]["max_depth"]) != (
        "soft_binary_24_v2",
        42,
    ):
        raise ValueError("Expected frozen soft42 preset")
    config["runtime"]["workers"] = workers
    validate_config(config)
    return config


def hashes(root: Path):
    return {
        str(p.relative_to(root)): sha256_file(p) for p in sorted(root.rglob("*")) if p.is_file()
    }


def prepare(root: Path, origin: Path, preset: Path, workers: int):
    config = execution_config(preset, workers)
    previous = origin / "repaired-mean"
    completion = json.loads((origin / "p6-completion.json").read_text())
    assert completion["job_id"] == "1410751" and completion["benchmark_unchanged"]
    assert completion["evaluation_completed"] and completion["strategy_parity_passed"]
    assert sha256_file(previous / "input/proteins.parquet") == (
        "e5323530d8f6acbdf6cab152e46869712e1965ff2cdda1b69300f6f245db106f"
    )
    old = yaml.safe_load((previous / "run.yaml").read_text())
    for section in ("input", "similarity", "edges", "components"):
        assert config[section] == old[section], f"Frozen upstream mismatch: {section}"
    a, b = dict(config["search"]), dict(old["search"])
    a.pop("threads")
    b.pop("threads")
    assert a == b, "Frozen search settings mismatch"
    proteins = pq.read_table(previous / "input/proteins.parquet")
    assert proteins.num_rows == 251378
    assert len(set(proteins["species_id"].to_pylist())) == 12
    assert len(set(proteins["original_id"].to_pylist())) == 251378
    run = root / "new-hierarchy"
    run.mkdir(exist_ok=False)
    fixed = {}
    for name in ("input", "edges", "components"):
        before = hashes(previous / name)
        shutil.copytree(previous / name, run / name)
        assert hashes(run / name) == before
        fixed.update({f"{name}/{p}": h for p, h in before.items()})
    (run / "run.yaml").write_text(yaml.safe_dump(config, sort_keys=True))
    shutil.copy2(preset, root / "selected-soft42.yaml")
    (root / "fixed-inputs.json").write_text(json.dumps(fixed, indent=2))
    shutil.copytree(origin / "benchmark", root / "benchmark")
    assert sha256_file(root / "benchmark/benchmark.py") == OFFICIAL_SCORER_SHA256
    (root / "benchmark-inputs.json").write_text(json.dumps(hashes(root / "benchmark"), indent=2))
    (root / "soft42-provenance.json").write_text(
        json.dumps(
            dict(
                origin=str(origin),
                origin_job="1410751",
                preset_sha256=sha256_file(preset),
                runtime_override={"runtime.workers": workers},
                original_workers=load_config(str(preset))["runtime"]["workers"],
                reused_stages=["input", "search", "edges", "components"],
                reused_hierarchy_components=0,
                fitting_to_reference=False,
                scope="fresh full hierarchy and OG extraction on frozen mean SSN",
            ),
            indent=2,
        )
    )


def freeze(root: Path):
    run = root / "new-hierarchy"
    require_resolved_hierarchy(run)
    fixed = json.loads((root / "fixed-inputs.json").read_text())
    assert all(sha256_file(run / p) == h for p, h in fixed.items())
    (root / "soft42-hierarchy-freeze.json").write_text(
        json.dumps(
            dict(
                resolved_gate_passed=True,
                configuration_sha256=sha256_file(run / "run.yaml"),
                manifests={
                    str(p.relative_to(run)): sha256_file(p)
                    for p in (run / "hierarchy/components").glob("*/hierarchy-manifest.json")
                },
            ),
            indent=2,
        )
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=("prepare", "freeze"))
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--origin", type=Path)
    p.add_argument("--preset", type=Path)
    p.add_argument("--workers", type=int, default=56)
    args = p.parse_args()
    if args.command == "prepare":
        if args.origin is None or args.preset is None:
            p.error("prepare requires --origin and --preset")
        prepare(args.run, args.origin, args.preset, args.workers)
    else:
        freeze(args.run)


if __name__ == "__main__":
    main()
