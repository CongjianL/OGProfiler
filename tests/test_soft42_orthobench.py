from pathlib import Path

import pytest
import yaml

from benchmarks.og_extraction.soft42_orthobench import execution_config
from ogprofiler.config import load_config

PRESET = Path(__file__).resolve().parents[1] / "presets/embleya-soft42.yaml"


def test_only_execution_workers_change():
    baseline = load_config(str(PRESET))
    current = execution_config(PRESET, 56)
    assert baseline["runtime"]["workers"] == 1
    assert current["runtime"]["workers"] == 56
    current["runtime"]["workers"] = 1
    assert current == baseline


def test_exact_preset_single_worker():
    assert execution_config(PRESET, 1) == load_config(str(PRESET))


def test_reject_invalid_workers():
    with pytest.raises(ValueError, match="positive"):
        execution_config(PRESET, 0)


def test_reject_wrong_policy(tmp_path):
    config = load_config(str(PRESET))
    config["hierarchy"]["max_depth"] = 20
    path = tmp_path / "wrong.yaml"
    path.write_text(yaml.safe_dump(config))
    with pytest.raises(ValueError, match="soft42"):
        execution_config(path, 56)
