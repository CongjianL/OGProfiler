from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa

from ogprofiler.benchmark.synthetic import (
    SyntheticScenario,
    generate_scenario_matrix,
    generate_synthetic_dataset,
    write_scenario_matrix,
)
from ogprofiler.benchmark.synthetic_metrics import (
    aggregate_applicability,
    evaluate_synthetic_run,
    write_synthetic_evaluation,
)


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_synthetic_matrix_is_full_factorial_and_stable(tmp_path: Path) -> None:
    scenarios = generate_scenario_matrix()
    assert len(scenarios) == 243
    assert len({item.scenario_id for item in scenarios}) == 243
    assert scenarios == generate_scenario_matrix()
    manifest, table = write_scenario_matrix(tmp_path, scenarios)
    assert json.loads(manifest.read_text())["design"] == "full-factorial"
    assert len(table.read_text().splitlines()) == 244


def test_synthetic_generator_records_complete_deterministic_truth(tmp_path: Path) -> None:
    scenario = SyntheticScenario(
        divergence=0.2,
        duplication_rate=0.5,
        loss_rate=0.1,
        expansion=2,
        fusion_rate=0.2,
    )
    left = tmp_path / "left"
    right = tmp_path / "right"
    generate_synthetic_dataset(left, scenario, ancestral_families=3, sequence_length=60)
    generate_synthetic_dataset(right, scenario, ancestral_families=3, sequence_length=60)
    left_manifest = json.loads((left / "manifest.json").read_text())
    right_manifest = json.loads((right / "manifest.json").read_text())
    assert left_manifest["artifact_checksums"] == right_manifest["artifact_checksums"]
    assert left_manifest["terminal_family_count"] == 6
    assert left_manifest["true_event_count"] > 0
    assert left_manifest["true_ortholog_pair_count"] > 0
    assert left_manifest["fusion_count"] > 0
    assert (left / "genealogy.tsv").is_file()
    assert (left / "true_events.tsv").is_file()
    assert (left / "true_orthologs.tsv").is_file()
    assert (left / "domain_architecture.tsv").is_file()


def _perfect_synthetic_fixture(root: Path) -> tuple[Path, Path]:
    dataset = root / "dataset"
    run = root / "run"
    scenario = SyntheticScenario()
    _write(
        dataset / "manifest.json",
        json.dumps(
            {
                "schema_version": "synthetic-evolution-dataset-v1",
                "scenario_id": scenario.scenario_id,
                "scenario": {
                    "divergence": 0.1,
                    "duplication_rate": 0.1,
                    "loss_rate": 0.05,
                    "expansion": 1,
                    "fusion_rate": 0.0,
                    "replicate": 0,
                    "seed": 42,
                },
            }
        ),
    )
    _write(
        dataset / "ground_truth.tsv",
        "species\tprotein_id\tfamily\trole\tlength\n"
        "s0\ta0\tF0\tsingle_copy\t60\ns1\ta1\tF0\tsingle_copy\t60\n"
        "s0\tb0\tF1\tsingle_copy\t60\ns1\tb1\tF1\tsingle_copy\t60\n",
    )
    _write(
        dataset / "genealogy.tsv",
        "family\tnode_id\tparent_id\tevent\tspecies\tprotein_id\tlineage\n"
        "AF\troot\t\tDUPLICATION\t\t\t-1\n"
        "F0\tleft\troot\tLINEAGE\t\t\t0\n"
        "F0\ta0\tleft\tGENE\ts0\ta0\t0\nF0\ta1\tleft\tGENE\ts1\ta1\t0\n"
        "F1\tright\troot\tLINEAGE\t\t\t1\n"
        "F1\tb0\tright\tGENE\ts0\tb0\t1\nF1\tb1\tright\tGENE\ts1\tb1\t1\n",
    )
    _write(
        dataset / "true_orthologs.tsv",
        "protein_a\tprotein_b\tfamily\na0\ta1\tF0\nb0\tb1\tF1\n",
    )
    _write(
        run / "results/families.tsv",
        "family_id\tcomponent_id\tcluster_id\tn_genes\tn_species\tterminal_reason\t"
        "network_event\nOG0\t0\t1\t2\t2\tNO_SPLIT\tAMBIGUOUS\n"
        "OG1\t0\t2\t2\t2\tNO_SPLIT\tAMBIGUOUS\n",
    )
    _write(
        run / "results/members.tsv",
        "family_id\tprotein_id\tspecies_id\toriginal_id\n"
        "OG0\t0\t0\ta0\nOG0\t1\t1\ta1\nOG1\t2\t0\tb0\nOG1\t3\t1\tb1\n",
    )
    _write(
        run / "results/hierarchy.tsv",
        "cluster_id\tparent_id\tcomponent_id\tdepth\tn_genes\tn_species\tresolution\t"
        "quality\tchild_count\tterminal_reason\n"
        "0\t\t0\t0\t4\t2\t1\t1\t2\t\n"
        "1\t0\t0\t1\t2\t2\t\t\t0\tNO_SPLIT\n"
        "2\t0\t0\t1\t2\t2\t\t\t0\tNO_SPLIT\n",
    )
    _write(
        run / "results/events.tsv",
        "component_id\tcluster_id\tnetwork_event\toverlap_score\tconfidence\n"
        "0\t0\tDUPLICATION_LIKE\t1\t1\n"
        "0\t1\tAMBIGUOUS\t0\t0\n0\t2\tAMBIGUOUS\t0\t0\n",
    )
    pairs = (
        "protein_a_id\tprotein_b_id\tspecies_a_id\tspecies_b_id\tcomponent_id\t"
        "supporting_cluster_id\trelationship\n"
        "0\t1\t0\t1\t0\t1\tCO_ORTHOLOG_CANDIDATE\n"
        "2\t3\t0\t1\t0\t2\tCO_ORTHOLOG_CANDIDATE\n"
    )
    with pa.output_stream(str(run / "results/ortholog_pairs.tsv.zst"), compression="zstd") as out:
        out.write(pairs.encode())
    return run, dataset


def test_synthetic_recovery_and_applicability_schema(tmp_path: Path) -> None:
    run, dataset = _perfect_synthetic_fixture(tmp_path)
    result = evaluate_synthetic_run(run, dataset)
    assert result["family"]["pairwise_clustering"]["f1"] == 1.0
    assert result["hierarchy"]["f1"] == 1.0
    assert result["events"]["end_to_end_accuracy"] == 1.0
    assert result["orthology"]["f1"] == 1.0
    metrics = tmp_path / "synthetic-metrics.json"
    write_synthetic_evaluation(metrics, result)
    applicability, table = aggregate_applicability([metrics], tmp_path / "map")
    value = json.loads(applicability.read_text())
    assert value["overall_applicable_scenarios"] == 1
    assert value["family_applicable_scenarios"] == 1
    assert value["rows"][0]["overall_applicable"] is True
    assert len(table.read_text().splitlines()) == 2
    axis = json.loads((tmp_path / "map/axis-summary.json").read_text())
    assert len(axis["rows"]) == 5
    assert (tmp_path / "map/axis-summary.tsv").is_file()
