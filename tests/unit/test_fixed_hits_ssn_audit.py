from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import pytest

pytest.importorskip("scipy")

from ogprofiler.similarity.models import DirectionalHit

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "fixed_hits_audit", ROOT / "benchmarks/compare_fixed_hits_ssn.py"
)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


@pytest.fixture
def reference_root():
    value = os.environ.get("OF_SOURCE_ROOT")
    if not value:
        pytest.skip("Set OF_SOURCE_ROOT to an installed OF3 orthofinder source directory")
    return Path(value)


def h(q, t, qs, ts, score):
    return DirectionalHit(q, t, qs, ts, score, 50.0, 100.0, 100.0, 1e-20)


def rows(species, lengths=None):
    lengths = lengths or {}
    return [
        dict(
            protein_id=p,
            species_id=s,
            original_id=f"protein_{p}",
            length=lengths.get(p, 100 + p * 10),
        )
        for p, s in species.items()
    ]


def test_fixed_hits_traces_all_layers_without_claiming_projection_equivalence(
    tmp_path, reference_root
):
    hits = [
        h(0, 2, 0, 1, 4.0),
        h(2, 0, 1, 0, 6.0),
        h(0, 3, 0, 1, 8.0),
        h(3, 0, 1, 0, 8.0),
        h(2, 4, 1, 2, 5.0),
        h(4, 2, 2, 1, 5.0),
        h(0, 2, 0, 1, 2.0),
    ]
    report = audit.run(
        hits, rows({0: 0, 1: 0, 2: 1, 3: 1, 4: 2}), reference_root, tmp_path / "audit"
    )
    assert report["completed"]
    assert report["hits"]["duplicate_usable"] == 1
    assert report["projection_contract"]["selected_projection"] == "mean"
    assert len(report["stages_complete"]) == 8
    assert all(
        x["outside_tolerance_count"] == 0
        for x in report["stage_comparisons"]
        if x["stage"].startswith("max_bitscore_input_diagnostic")
    )
    assert all(
        x["outside_tolerance_count"] == 0
        for x in report["stage_comparisons"]
        if x["stage"].startswith("B/")
    )
    assert all(
        x["outside_tolerance_count"] == 0
        for x in report["stage_comparisons"]
        if x["stage"].startswith(("BH_sameB_control/", "RBH_sameB_control/"))
    )
    assert (tmp_path / "audit/reference/orthofinder/tools/waterfall.py").read_bytes() == (
        reference_root / "tools/waterfall.py"
    ).read_bytes()


def test_true_identity_zero_diff_and_symmetric_lift_double_direction(tmp_path, reference_root):
    # Two hits per directed pair avoid the too-few-fit-points branch. Equal
    # scores/length products make normalized forward/reverse weights equal.
    hits = [h(0, 2, 0, 1, 8.0), h(1, 3, 0, 1, 8.0), h(2, 0, 1, 0, 8.0), h(3, 1, 1, 0, 8.0)]
    report = audit.run(
        hits,
        rows({0: 0, 1: 0, 2: 1, 3: 1}, {0: 100, 1: 200, 2: 100, 3: 200}),
        reference_root,
        tmp_path / "audit",
    )
    assert report["graph_objects"]["V2_undirected"]["undirected_edges"] == 2
    assert report["graph_objects"]["V2_undirected"]["directed_nonzero"] == 4
    diff = next(
        x
        for x in report["stage_comparisons"]
        if x["stage"] == "OF_directional_W_vs_V2_symmetric_lift"
    )
    assert diff["outside_tolerance_count"] == 0
    assert diff["max_abs_error"] == 0


def test_preserves_duplicate_original_ids_across_species_and_isolates(tmp_path, reference_root):
    proteins = rows({10: 7, 20: 12, 30: 12})
    for row in proteins:
        row["original_id"] = "same_id"
    report = audit.run([], proteins, reference_root, tmp_path / "audit")
    assert report["graph_objects"]["V2_undirected"]["singletons"] == 3
    assert report["graph_objects"]["OF_directional_full_precision"]["singletons"] == 3


def test_fresh_output_is_required(tmp_path, reference_root):
    out = tmp_path / "existing"
    out.mkdir()
    with pytest.raises(FileExistsError):
        audit.run([], rows({0: 0}), reference_root, out)


def test_repaired_pipeline_matches_of_and_mean_projection(tmp_path, reference_root):
    # Equal length products, unsorted hits, incomplete bins, duplicate HSPs,
    # within-species hits and asymmetric directions exercise the full path.
    import random

    rng = random.Random(17)
    species = {p: p // 23 for p in range(46)}
    proteins = rows(species, {p: 100 + (p % 5) * 20 for p in species})
    hits = [
        h(q, t, species[q], species[t], float(rng.randint(20, 500)))
        for q in species
        for t in species
        if q != t
    ]
    hits += [h(0, 24, 0, 1, 1.0), h(1, 25, 0, 1, 1.0)]
    rng.shuffle(hits)
    report = audit.run(hits, proteins, reference_root, tmp_path / "audit")
    expected_equal = {
        "max_bitscore_input_diagnostic",
        "B",
        "BH",
        "RBH_cross_species",
        "BH_sameB_control",
        "RBH_sameB_control",
        "cutoff",
        "cutoff_sameB_control",
        "connect",
        "connect_sameB_control",
        "directional_W_same_assembly_DIAGNOSTIC",
        "OF_mean_projection_vs_V2_undirected",
    }
    for comparison in report["stage_comparisons"]:
        if comparison["stage"].split("/")[0] in expected_equal:
            assert comparison["outside_tolerance_count"] == 0, comparison
    for pair in report["normalization_pairs"]:
        assert pair["fit_sample_only_v2"] == pair["fit_sample_only_of"] == 0
        assert pair["v2_parameters"] == pytest.approx(pair["of_parameters"])
