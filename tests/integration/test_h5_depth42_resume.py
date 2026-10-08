"""H5 completion gate and reused component integrity on small public fixtures."""

import json
from pathlib import Path

import pytest
from test_hierarchy_scheduler import _prepare_multicomponent_run

from benchmarks.og_extraction.h5_depth42_resume import freeze, prepare
from ogprofiler.config import hierarchy_config, load_config
from ogprofiler.hierarchy.stage import run_hierarchy_component_stage


def test_h5_freeze_requires_complete_hierarchy_and_preserves_reused_component(tmp_path):
    run = _prepare_multicomponent_run(tmp_path)
    root = tmp_path / "campaign"
    root.mkdir()
    run.rename(root / "new-hierarchy")
    run = root / "new-hierarchy"
    config = hierarchy_config(
        load_config(
            overrides=["hierarchy.topology_policy=soft_binary_24_v2", "hierarchy.max_depth=42"]
        )["hierarchy"]
    )
    folder, _ = run_hierarchy_component_stage(run, 0, config, [])
    manifest = json.loads((folder / "hierarchy-manifest.json").read_text())
    (root / "h5-reuse-provenance.json").write_text(
        json.dumps(dict(component0_output_checksums=manifest["output_checksums"]))
    )
    with pytest.raises(Exception, match="Missing hierarchy component 1"):
        freeze(root)
    run_hierarchy_component_stage(run, 1, config, [])
    freeze(root)
    assert json.loads((root / "h5-hierarchy-freeze.json").read_text())["component0_unchanged"]
    (folder / "metrics.json").write_text("{}")
    with pytest.raises(Exception, match="checksum mismatch"):
        freeze(root)


def test_h5_prepare_rejects_unaccepted_h4_before_creating_run(tmp_path):
    baseline = tmp_path / "baseline"
    baseline.mkdir()
    (baseline / "h4-acceptance.json").write_text('{"passed": false}')
    (baseline / "depth-prefix-acceptance.json").write_text('{"passed": false}')
    root = tmp_path / "campaign"
    root.mkdir()
    with pytest.raises(AssertionError):
        prepare(root, baseline, Path("/unused"))
    assert not (root / "new-hierarchy").exists()
