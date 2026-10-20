from __future__ import annotations

import os
import runpy
from pathlib import Path

import pytest

from benchmarks.og_extraction.orthobench import (
    OFFICIAL_SCORER_SHA256,
    refog_diagnostics,
    validate_mapping,
)
from ogprofiler.core.manifest import sha256_file


def test_mapping_retains_missing_predictions():
    assert validate_mapping({"a": ["p"]}, {"p", "q"}) == {"p"}


@pytest.mark.parametrize("groups", [{"a": ["p"], "b": ["p"]}, {"a": ["x"]}, {}])
def test_mapping_rejects_corruption(groups):
    with pytest.raises(ValueError):
        validate_mapping(groups, {"p"})


def test_diagnostics_keep_uncertainty_separate():
    result = refog_diagnostics(
        {"a": ["p", "u", "extra"], "b": ["q"]}, {"R": {"p", "q", "u", "missing"}}, {"R": {"u"}}
    )[0]
    assert result["confident_genes"] == 3
    assert result["fragments"] == 2 and result["split_excess"] == 1
    assert result["best_group_extra_genes"] == 0
    assert result["contaminated_fragments"] == 1
    assert result["missing_confident"] == 1
    assert result["covered_raw"] == 3
    assert result["classification"] == "MERGE_AND_SPLIT"


def test_actual_official_function_uncertainty_and_missing():
    source = os.environ.get("ORTHOBENCH_SOURCE")
    if not source:
        pytest.skip("Set ORTHOBENCH_SOURCE to the supplied benchmark.py")
    path = Path(source)
    assert sha256_file(path) == OFFICIAL_SCORER_SHA256
    official = runpy.run_path(str(path))
    # After removing uncertain u: reference p,q,r and prediction p,q,x.
    # TP=1, FP=2, FN=2 => P=R=F1=1/3 (per-RefOG normalization cancels).
    values = official["calculate_benchmarks_pairwise"](
        [{"p", "q", "r", "u"}], [{"u"}], [{"p", "q", "u", "x"}]
    )
    assert values == pytest.approx([100 / 3] * 3)
