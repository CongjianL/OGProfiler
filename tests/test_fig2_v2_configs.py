import runpy
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MOD = runpy.run_path(
    str(ROOT / "OGProfiler2_benchmark/scripts/statistics/update_fig2_v2_configs.py")
)


def test_raw_refog_diagnostics_and_official_exclusion_are_distinct():
    groups = {"a": {"x", "y", "low", "extra"}, "b": {"z"}}
    refs, fp, pairwise = MOD["diagnostics"](
        groups, {"RefOG001": {"x", "y", "z", "low"}}, {"RefOG001": {"low"}}
    )
    row = refs[0]
    assert row["F1"] == 0.75
    assert row["split_count"] == 1
    assert row["contamination"] == 0.25
    assert row["missing_fraction"] == 0.25
    assert not row["exact"] == "true"
    assert pairwise["precision"] == pytest.approx(1 / 3)
    assert pairwise["recall"] == pytest.approx(1 / 3)
    assert fp == {1: 1.0, 5: 1.0, 10: 1.0}


def test_raw_tie_breaking_matches_existing_metric_script():
    refs, _, _ = MOD["diagnostics"]({"a": {"x"}, "b": {"y"}}, {"RefOG001": {"x", "y"}}, {})
    assert refs[0]["best_predicted_group"] == "b"


def test_reject_duplicate_assignment(tmp_path):
    path = tmp_path / "prediction.txt"
    path.write_text("a: x y\nb: y z\n")
    with pytest.raises(ValueError, match="Duplicate"):
        MOD["predictions"](path)
