from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from ogprofiler.cli import main
from ogprofiler.hierarchy.loader import ComponentGraphLoader
from ogprofiler.similarity.io import write_retained_edges
from ogprofiler.similarity.models import RetainedEdge


def test_partitioned_component_to_production_hierarchy_and_resume(tmp_path: Path) -> None:
    proteomes = tmp_path / "proteomes"
    proteomes.mkdir()
    for species in range(3):
        records = "".join(
            f">g{species}_{index}\n{'A' * 20}\n" for index in range(5)
        )
        (proteomes / f"species_{species}.faa").write_text(records, encoding="utf-8")
    run = tmp_path / "run"
    assert main(["prepare", "--proteomes", str(proteomes), "--out", str(run)]) == 0

    edges: list[RetainedEdge] = []
    for group in range(3):
        vertices = range(group * 5, (group + 1) * 5)
        for left in vertices:
            for right in vertices:
                if left < right:
                    edges.append(
                        RetainedEdge(left, right, left // 5, right // 5, 10, 10, 10, 100, "TEST")
                    )
    edges.extend(
        [
            RetainedEdge(4, 5, 0, 1, 0.01, 0.01, 0.01, 100, "TEST"),
            RetainedEdge(9, 10, 1, 2, 0.01, 0.01, 0.01, 100, "TEST"),
        ]
    )
    write_retained_edges(run / "edges" / "retained_edges.parquet", edges)
    assert main(["components", "--run", str(run)]) == 0

    loaded = ComponentGraphLoader(run).load(0)
    assert loaded.local_to_global == tuple(range(15))
    assert loaded.global_to_local[14] == 14
    assert loaded.species_bitmap_by_protein[10] == 4
    assert isinstance(loaded.global_ids, np.memmap)
    assert isinstance(loaded.species_ids, np.memmap)
    assert loaded.global_ids.tolist() == list(range(15))

    command = [
        "hierarchy",
        "--run",
        str(run),
        "--component-id",
        "0",
        "--set",
        "hierarchy.stability_mode=fast",
        "--set",
        "hierarchy.gamma_min=0.1",
        "--set",
        "hierarchy.gamma_max=2.0",
        "--set",
        "hierarchy.max_child_fraction=0.8",
        "--set",
        "hierarchy.subtree_workers=2",
        "--set",
        "hierarchy.subtree_release_size=1",
    ]
    assert main(command) == 0
    output = run / "hierarchy" / "components" / "component=00000000"
    nodes_path = output / "nodes.parquet"
    first_mtime = nodes_path.stat().st_mtime_ns
    assert main(command) == 0
    assert nodes_path.stat().st_mtime_ns == first_mtime

    nodes = pq.read_table(nodes_path).to_pylist()
    assert nodes[0]["child_count"] == 3
    assert nodes[0]["n_species"] == 3
    assert pq.read_table(output / "members.parquet").num_rows == 15
    candidates = pq.read_table(output / "candidates.parquet")
    assert "stability" in candidates.column_names
    assert "tiny_fragment_fraction" in candidates.column_names
    manifest = json.loads((output / "hierarchy-manifest.json").read_text(encoding="utf-8"))
    assert manifest["algorithm_version"] == "hierarchical-leiden-v1"

    nodes_path.write_bytes(b"corrupt")
    assert main(command) == 0
    assert pq.read_table(nodes_path).num_rows == len(nodes)

    annotation_command = [
        "annotate-network",
        "--run",
        str(run),
        "--set",
        "evolution.network_overlap_threshold=0.5",
    ]
    assert main(annotation_command) == 0
    events_path = run / "evolution" / "components" / "component=00000000/events.parquet"
    first_event_mtime = events_path.stat().st_mtime_ns
    assert main(annotation_command) == 0
    assert events_path.stat().st_mtime_ns == first_event_mtime
    event_rows = pq.read_table(events_path).to_pylist()
    assert event_rows[0]["network_event"] == "POLYTOMY"
    assert event_rows[0]["legacy_event"] == "III-3"
    assert {row["network_event"] for row in event_rows[1:]} == {"SPECIES_SPECIFIC"}
    events_path.write_bytes(b"corrupt")
    assert main(annotation_command) == 0
    assert pq.read_table(events_path).num_rows == 4

    export_command = ["export", "--run", str(run), "--family-fasta", "OG000000000"]
    assert main(export_command) == 0
    results = run / "results"
    families_before = (results / "families.tsv").read_bytes()
    members_before = (results / "members.tsv").read_bytes()
    assert families_before.decode().splitlines()[0] == (
        "family_id\tcomponent_id\tcluster_id\tn_genes\tn_species\t"
        "terminal_reason\tnetwork_event"
    )
    assert len(families_before.decode().splitlines()) == 4
    assert len(members_before.decode().splitlines()) == 16
    assert (results / "hierarchy.tsv").is_file()
    assert (results / "events.tsv").is_file()
    assert (results / "fasta" / "OG000000000.faa").is_file()
    manifest = json.loads((results / "export-manifest.json").read_text(encoding="utf-8"))
    assert manifest["counts"] == {
        "events": 4,
        "families": 3,
        "fasta_files": 1,
        "hierarchy_nodes": 4,
        "members": 15,
    }

    first_family_mtime = (results / "families.tsv").stat().st_mtime_ns
    assert main(export_command) == 0
    assert (results / "families.tsv").stat().st_mtime_ns == first_family_mtime

    (results / "families.tsv").write_bytes(b"corrupt")
    assert main(export_command) == 0
    assert (results / "families.tsv").read_bytes() == families_before

    member_table = pq.read_table(output / "members.parquet")
    reordered = member_table.take(list(reversed(range(member_table.num_rows))))
    pq.write_table(reordered, output / "members.parquet")
    assert main(export_command) == 0
    assert (results / "families.tsv").read_bytes() == families_before
    assert (results / "members.tsv").read_bytes() == members_before

    assert main(["export", "graph", "--run", str(run), "--component", "0"]) == 0
    graphml = (results / "graphs" / "component=00000000.graphml").read_text(
        encoding="utf-8"
    )
    assert graphml.count("<node id=") == 15
    assert graphml.count("<edge id=") == len(edges)

    assert main(["orthologs", "--run", str(run)]) == 0
    ortholog_path = results / "ortholog_pairs.tsv.zst"
    assert not ortholog_path.exists()
    ortholog_command = [
        "orthologs",
        "--run",
        str(run),
        "--emit-pairwise-orthologs",
        "--set",
        "output.ortholog_pair_chunk_size=2",
    ]
    assert main(ortholog_command) == 0
    ortholog_mtime = ortholog_path.stat().st_mtime_ns
    ortholog_manifest = json.loads(
        (results / "ortholog-manifest.json").read_text(encoding="utf-8")
    )
    assert ortholog_manifest["counts"] == {
        "chunks": 38,
        "components": 1,
        "pairs": 75,
        "supporting_event_counts": {"POLYTOMY": 1, "SPECIATION_LIKE": 0},
        "supporting_nodes": 1,
    }
    assert main(ortholog_command) == 0
    assert ortholog_path.stat().st_mtime_ns == ortholog_mtime
    ortholog_path.write_bytes(b"corrupt")
    assert main(ortholog_command) == 0
    assert ortholog_path.read_bytes() != b"corrupt"

    fake_mafft = tmp_path / "fake-mafft"
    fake_mafft.write_text(
        "#!/bin/bash\n"
        "if [[ $1 == --version ]]; then echo 'fake-mafft 1' >&2; exit 0; fi\n"
        "/bin/cat \"${!#}\"\n",
        encoding="utf-8",
    )
    fake_tree = tmp_path / "fake-fasttree"
    fake_tree.write_text(
        "#!/bin/bash\n"
        "if [[ $# == 0 ]]; then echo 'fake-fasttree 1' >&2; exit 0; fi\n"
        "echo '(P000000000000:1,P000000000001:1,P000000000002:1,"
        "P000000000003:1,P000000000004:1);'\n",
        encoding="utf-8",
    )
    fake_mafft.chmod(0o755)
    fake_tree.chmod(0o755)
    phylogeny_command = [
        "annotate",
        "--run",
        str(run),
        "--phylogenetic-refinement",
        "--family",
        "OG000000000",
        "--set",
        f"phylogeny.alignment_executable={fake_mafft}",
        "--set",
        f"phylogeny.tree_executable={fake_tree}",
    ]
    assert main(phylogeny_command) == 0
    phylogeny_root = run / "evolution/phylogenetic"
    phylogeny_events = (phylogeny_root / "phylogenetic-events.tsv").read_text(
        encoding="utf-8"
    )
    phylogeny_manifest = json.loads(
        (phylogeny_root / "phylogenetic-manifest.json").read_text(encoding="utf-8")
    )
    assert phylogeny_manifest["counts"]["selected"] == 1
    assert "SPECIES_SPECIFIC\tDUPLICATION" in phylogeny_events
    rooted = phylogeny_root / "family=OG000000000/gene_tree.rooted.nwk"
    rooted_mtime = rooted.stat().st_mtime_ns
    assert main(phylogeny_command) == 0
    assert rooted.stat().st_mtime_ns == rooted_mtime
