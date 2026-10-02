from __future__ import annotations

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from test_og_v1_events import FIXTURES, component_inputs

from benchmarks.og_extraction.fixed_hierarchy import audit, compare_component
from ogprofiler.evolution.stage import run_network_annotation_stage


@pytest.mark.parametrize("spec", FIXTURES, ids=lambda spec: spec["name"])
def test_real_hierarchy_adapter_against_frozen_fixtures(spec):
    genes = sorted({g for v in spec["vertices"] for g in v["genes"]})
    for nodes, members, species, _ in component_inputs(spec):
        originals = {p: genes[p].split("|", 1)[1] for p in species}
        isolates = {p for p in species if genes[p] in spec.get("isolates", [])}
        _, report = compare_component(
            nodes,
            members,
            species,
            originals,
            spec["n_species"],
            isolates,
            spec.get("overlap_count", 0),
        )
        assert all(report["checks"].values()), report
        assert report["passed"] == (not spec["expected"]["duplicate_members"])


def test_p5_audits_actual_parquet_and_tsv(tmp_path):
    spec = next(s for s in FIXTURES if s["name"] == "binary_disjoint")
    nodes, members, species, _ = next(component_inputs(spec))
    run = tmp_path / "run"

    def put(name, records, schema=None):
        path = run / name
        path.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.Table.from_pylist(records, schema=schema), path)

    component = nodes[0]["component_id"]
    folder = f"hierarchy/components/component={component:08d}"
    for node in nodes:
        node["child_count"] = sum(n["parent_id"] == node["cluster_id"] for n in nodes)
        node["terminal_reason"] = None if node["child_count"] else "TEST"
    put(folder + "/nodes.parquet", nodes)
    put(folder + "/members.parquet", members)
    put(
        "input/proteins.parquet",
        [dict(protein_id=p, species_id=s, original_id=f"g{p}") for p, s in species.items()],
    )
    put("input/species.parquet", [dict(species_id=s) for s in range(spec["n_species"])])
    (run / "input/proteins.faa").write_text("".join(f">OGP2P{p:012d}\nAAAA\n" for p in species))
    put("components/index.parquet", [dict(component_id=component, protein_id=p) for p in species])
    put(
        "components/singleton_terminal_families.parquet",
        [],
        pa.schema([("component_id", pa.int64()), ("protein_id", pa.int64())]),
    )
    run_network_annotation_stage(run, 0.0, [])
    report = audit(run, tmp_path / "audit")
    assert report["passed"], report
    assert report["artifacts_passed"] and report["frozen_inputs_unchanged"]
