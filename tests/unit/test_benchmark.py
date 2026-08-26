from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa

from ogprofiler.benchmark.comparison import compare_methods
from ogprofiler.benchmark.matrix import generate_ofat_matrix, write_matrix
from ogprofiler.benchmark.metrics import (
    evaluate_legacy_membership,
    evaluate_run,
    write_evaluation,
)


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _perfect_run(root: Path) -> Path:
    results = root / "results"
    _write(
        results / "families.tsv",
        "family_id\tcomponent_id\tcluster_id\tn_genes\tn_species\tterminal_reason\t"
        "network_event\nOG000000000\t0\t1\t4\t4\tNO_SPLIT\tAMBIGUOUS\n",
    )
    _write(
        results / "members.tsv",
        "family_id\tprotein_id\tspecies_id\toriginal_id\n"
        "OG000000000\t0\t0\ts0|a\nOG000000000\t1\t1\ts1|a\n"
        "OG000000000\t2\t2\ts2|a\nOG000000000\t3\t3\ts3|a\n",
    )
    _write(
        results / "hierarchy.tsv",
        "cluster_id\tparent_id\tcomponent_id\tdepth\tn_genes\tn_species\tresolution\t"
        "quality\tchild_count\tterminal_reason\n"
        "0\t\t0\t0\t4\t4\t1\t1\t4\t\n"
        "1\t0\t0\t1\t1\t1\t\t\t0\tONE_SPECIES\n"
        "2\t0\t0\t1\t1\t1\t\t\t0\tONE_SPECIES\n"
        "3\t0\t0\t1\t1\t1\t\t\t0\tONE_SPECIES\n"
        "4\t0\t0\t1\t1\t1\t\t\t0\tONE_SPECIES\n",
    )
    _write(
        results / "events.tsv",
        "component_id\tcluster_id\tnetwork_event\toverlap_score\tconfidence\n"
        "0\t0\tPOLYTOMY\t0\t1\n0\t1\tSPECIES_SPECIFIC\t0\t1\n"
        "0\t2\tSPECIES_SPECIFIC\t0\t1\n0\t3\tSPECIES_SPECIFIC\t0\t1\n"
        "0\t4\tSPECIES_SPECIFIC\t0\t1\n",
    )
    truth = root / "ground_truth.tsv"
    _write(
        truth,
        "species\tprotein_id\tfamily\trole\tlength\n"
        "s0\ts0|a\tF0\tone_to_one\t100\n"
        "s1\ts1|a\tF0\tone_to_one\t100\n"
        "s2\ts2|a\tF0\tone_to_one\t100\n"
        "s3\ts3|a\tF0\tone_to_one\t100\n",
    )
    rows = [
        "protein_a_id\tprotein_b_id\tspecies_a_id\tspecies_b_id\tcomponent_id\t"
        "supporting_cluster_id\trelationship\n"
    ]
    for left in range(4):
        for right in range(left + 1, 4):
            rows.append(f"{left}\t{right}\t{left}\t{right}\t0\t0\tCO_ORTHOLOG_CANDIDATE\n")
    with pa.output_stream(str(results / "ortholog_pairs.tsv.zst"), compression="zstd") as out:
        out.write("".join(rows).encode())
    return truth


def test_parameter_matrix_covers_every_required_axis(tmp_path: Path) -> None:
    runs = generate_ofat_matrix()
    assert len(runs) == 16
    assert runs[0].run_id == "baseline"
    assert {run.axis for run in runs} == {
        "baseline",
        "normalization",
        "coverage",
        "symmetrization",
        "leiden_method",
        "gamma_strategy",
        "split_acceptance",
        "seed",
    }
    manifest, table = write_matrix(tmp_path, runs)
    assert json.loads(manifest.read_text())["run_count"] == 16
    assert len(table.read_text().splitlines()) == 17


def test_scientific_metrics_cover_family_evolution_and_orthology(tmp_path: Path) -> None:
    truth = _perfect_run(tmp_path)
    result = evaluate_run(tmp_path, truth, method="ogprofiler2", dataset="fixture")
    assert result["family"]["pairwise_clustering"]["f1"] == 1.0
    assert result["family"]["single_copy_family_recovery"]["rate"] == 1.0
    assert result["evolution"]["species_overlap_consistency"] == 1.0
    assert result["orthology"]["precision"] == 1.0
    assert result["orthology"]["recall"] == 1.0


def test_method_comparison_uses_closed_schema_and_ranks_scores(tmp_path: Path) -> None:
    truth = _perfect_run(tmp_path / "run")
    strong = evaluate_run(
        tmp_path / "run", truth, method="ogprofiler2", dataset="fixture"
    )
    weak = json.loads(json.dumps(strong))
    weak["method"] = "ogprofiler1"
    weak["family"]["pairwise_clustering"]["f1"] = 0.5
    strong_path, weak_path = tmp_path / "strong.json", tmp_path / "weak.json"
    write_evaluation(strong_path, strong)
    write_evaluation(weak_path, weak)
    comparison, table = compare_methods([weak_path, strong_path], tmp_path / "comparison")
    value = json.loads(comparison.read_text())
    assert value["rows"][0]["method"] == "ogprofiler2"
    assert value["method_registry"]["orthofinder"]["status"] == "planned-import"
    assert len(table.read_text().splitlines()) == 3


def test_legacy_membership_uses_same_family_metric_schema(tmp_path: Path) -> None:
    truth = _perfect_run(tmp_path / "run")
    normalized = tmp_path / "legacy.json"
    normalized.write_text(
        json.dumps(
            {
                "membership": {
                    "s0|a": "legacy_family",
                    "s1|a": "legacy_family",
                    "s2|a": "legacy_family",
                    "s3|a": "legacy_family",
                }
            }
        ),
        encoding="utf-8",
    )
    result = evaluate_legacy_membership(
        normalized, truth, method="ogprofiler1", dataset="fixture"
    )
    assert result["family"]["pairwise_clustering"]["f1"] == 1.0
    assert result["family"]["exact_family_recovery"]["rate"] == 1.0
    assert result["evolution"]["species_overlap_consistency"] is None
    assert result["orthology"] is None
