from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from test_og_v1_events import component_inputs

from ogprofiler.evolution.stage import run_network_annotation_stage
from ogprofiler.exceptions import ExportError, HierarchyError
from ogprofiler.orthogroups.stage import OrthogroupConfig, run_orthogroup_stage
from ogprofiler.output.orthogroups import assemble_orthogroup_export_tables
from ogprofiler.output.stage import run_export_stage

CASES = json.loads(
    (Path(__file__).parents[1] / "fixtures/og_extraction/v1_contract.json").read_text()
)["cases"]


def put(path, rows, schema=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows, schema=schema), path)


@pytest.mark.parametrize("spec", CASES, ids=lambda spec: spec["name"])
def test_persisted_member_sets_match_independent_frozen_v1(tmp_path, spec):
    genes = sorted({gene for vertex in spec["vertices"] for gene in vertex["genes"]})
    species = sorted({s for vertex in spec["vertices"] for s in vertex["species"]})
    gene_ids = {gene: i for i, gene in enumerate(genes)}
    index, isolates = [], []
    for nodes, members, _, _ in component_inputs(spec):
        component = nodes[0]["component_id"]
        directory = tmp_path / "hierarchy/components" / f"component={component:08d}"
        for node in nodes:
            child_count = sum(row["parent_id"] == node["cluster_id"] for row in nodes)
            node["child_count"] = child_count
            node["terminal_reason"] = None if child_count else "REFERENCE_TERMINAL"
        put(directory / "nodes.parquet", nodes)
        put(directory / "members.parquet", members)
        put(
            tmp_path / "evolution/components" / directory.name / "events.parquet",
            [
                dict(
                    component_id=component,
                    cluster_id=nodes[0]["cluster_id"],
                    network_event="FIXTURE",
                ),
            ],
        )
        index.extend(dict(component_id=component, protein_id=m["protein_id"]) for m in members)
        isolates.extend(
            dict(component_id=component, protein_id=m["protein_id"])
            for m in members
            if genes[m["protein_id"]] in spec.get("isolates", [])
        )
    put(
        tmp_path / "input/proteins.parquet",
        [
            dict(
                protein_id=i,
                species_id=species.index(g.split("|")[0]),
                original_id=g.split("|", 1)[1],
            )
            for g, i in gene_ids.items()
        ],
    )
    put(tmp_path / "input/species.parquet", [dict(species_id=i) for i in range(spec["n_species"])])
    put(tmp_path / "components/index.parquet", index)
    put(
        tmp_path / "components/singleton_terminal_families.parquet",
        isolates,
        pa.schema([("component_id", pa.int64()), ("protein_id", pa.int64())]),
    )
    if not spec["expected"]["duplicate_members"]:
        run_network_annotation_stage(tmp_path, 0.0, [])
    (tmp_path / "input/proteins.faa").write_text(
        "".join(f">OGP2P{i:012d}\nAAAA\n" for i in range(len(genes)))
    )
    config = OrthogroupConfig(species_overlap_count=spec.get("overlap_count", 0))
    if spec["expected"]["duplicate_members"]:
        with pytest.raises(HierarchyError):
            run_orthogroup_stage(tmp_path, config, [], retries=0)
        with pytest.raises(ExportError):
            run_export_stage(tmp_path, [], config=config)
        failures = list((tmp_path / "orthogroups/components").glob("*/og-failure.json"))
        assert failures
        assert sum(json.loads(path.read_text())["counts"]["duplicate"] for path in failures) == len(
            spec["expected"]["duplicate_members"]
        )
        return
    run_orthogroup_stage(tmp_path, config, [])
    actual, unassigned = [], []
    for directory in (tmp_path / "orthogroups/components").glob("component=*"):
        groups = pq.ParquetFile(directory / "groups.parquet").read().to_pylist()
        members = pq.ParquetFile(directory / "members.parquet").read().to_pylist()
        for group in groups:
            actual.append(
                (
                    group["processing_level"],
                    tuple(
                        sorted(
                            genes[m["protein_id"]]
                            for m in members
                            if m["local_group_id"] == group["local_group_id"]
                        )
                    ),
                )
            )
        unassigned.extend(
            genes[m["protein_id"]]
            for m in pq.ParquetFile(directory / "unassigned.parquet").read().to_pylist()
        )
    assert sorted(actual) == sorted(
        (g["level"], tuple(g["members"])) for g in spec["expected"]["groups"]
    )
    assert sorted(unassigned) == spec["expected"]["unassigned"]

    tables = assemble_orthogroup_export_tables(tmp_path, config)
    exported = []
    for group in tables.groups:
        exported.append(
            (group.processing_level, tuple(sorted(genes[m.protein_id] for m in group.members)))
        )
    assert sorted(exported) == sorted(
        (g["level"], tuple(g["members"])) for g in spec["expected"]["groups"]
    )
    canonical_expected = sorted(
        spec["expected"]["groups"],
        key=lambda group: tuple(
            sorted((species.index(g.split("|")[0]), g.split("|", 1)[1]) for g in group["members"])
        ),
    )
    assert [
        (g.family_id, tuple(sorted(genes[m.protein_id] for m in g.members))) for g in tables.groups
    ] == [
        (f"OG{rank:09d}", tuple(group["members"])) for rank, group in enumerate(canonical_expected)
    ]
    run_export_stage(tmp_path, [], config=config)
    assert (tmp_path / "results/unassigned.tsv").is_file()
