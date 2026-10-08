"""Frozen-path audit: new contamination is not a same-leaf mixing claim."""

import json

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from benchmarks.og_extraction.refog_pollution_audit import audit
from ogprofiler.core.manifest import sha256_file


def fixture(root, merged):
    run = root / "new-hierarchy"
    hierarchy = run / "hierarchy/components/component=00000000"
    og = run / "orthogroups/components/component=00000000"
    for folder in (
        hierarchy,
        og,
        run / "input",
        run / "components",
        run / "results",
        root / "benchmark/RefOGs",
    ):
        folder.mkdir(parents=True)

    def write(folder, name, data):
        pq.write_table(pa.Table.from_pylist(data), folder / name)

    genes = ["A", "B", "X"]
    write(
        run / "input",
        "proteins.parquet",
        [dict(protein_id=i, original_id=g) for i, g in enumerate(genes)],
    )
    write(
        run / "components", "index.parquet", [dict(protein_id=i, component_id=0) for i in range(3)]
    )
    (root / "benchmark/RefOGs/RefOG001.txt").write_text("A\nB\n")
    leaves = [1, 2, 2] if merged else [1, 2, 3]
    nodes = [dict(cluster_id=0, parent_id=None, depth=0, child_count=len(set(leaves)), n_genes=3)]
    nodes += [
        dict(cluster_id=i, parent_id=0, depth=1, child_count=0, n_genes=leaves.count(i))
        for i in sorted(set(leaves))
    ]
    write(hierarchy, "nodes.parquet", nodes)
    write(
        hierarchy,
        "members.parquet",
        [dict(protein_id=i, terminal_cluster_id=leaf) for i, leaf in enumerate(leaves)],
    )
    write(hierarchy, "candidates.parquet", [dict(cluster_id=0, selected=True)])
    (hierarchy / "metrics.json").write_text("{}")
    write(
        og,
        "v1_events.parquet",
        [
            dict(
                cluster_id=n["cluster_id"],
                v1_event=("I" if merged else "III-3") if n["cluster_id"] == 0 else None,
            )
            for n in nodes
        ],
    )
    write(
        og,
        "selection_trace.parquet",
        [dict(cluster_id=0, status="SELECTED" if merged else "REJECTED")],
    )
    write(og, "groups.parquet", [dict(local_group_id=0, source_cluster_id=0, v1_event="I")])
    write(
        og,
        "members.parquet",
        [dict(protein_id=i, local_group_id=0 if merged else i) for i in range(3)],
    )
    (run / "results/members.tsv").write_text(
        "family_id\toriginal_id\n"
        + "".join(f"G{0 if merged else i}\t{g}\n" for i, g in enumerate(genes))
    )
    for folder, manifest in [(hierarchy, "hierarchy-manifest.json"), (og, "og-manifest.json")]:
        (folder / manifest).write_text(
            json.dumps(dict(output_checksums={p.name: sha256_file(p) for p in folder.iterdir()}))
        )
    if merged:
        d = root / "h5-metrics/new_bounded_hierarchy_v1_compatible"
        d.mkdir(parents=True)
        (d / "refog-diagnostics.tsv").write_text(
            "refog\tbest_group\tbest_group_extra_genes\tclassification\nRefOG001\tG0\t1\tOVERMERGE\n"
        )


def test_audit_locates_binary_i_rejoining_and_baseline_first_divergence(tmp_path):
    current, baseline = tmp_path / "current", tmp_path / "baseline"
    fixture(current, True)
    fixture(baseline, False)
    report = audit(current, baseline)
    (case,) = report["cases"]
    assert case["extra_genes"] == ["X"]
    assert case["newly_joined_pairs"] == 2
    assert case["baseline_equivalent_source_clades"][0]["event"] == "III-3"
    assert case["group"]["v1_event"] == "I"
    for pair in case["pairs"]:
        assert (
            pair["mechanism"]
            == {"A": "OG_REJOINS_DISTINCT_TERMINALS", "B": "SAME_TERMINAL"}[pair["left"]]
        )
        assert pair["first_partition_divergence"]["current_node"]["child_count"] == 2
        assert pair["first_partition_divergence"]["baseline_node"]["child_count"] == 3
    # Input mismatch must block positional protein-ID comparison across runs.
    (baseline / "new-hierarchy/components/index.parquet").write_bytes(b"changed")
    with pytest.raises(AssertionError):
        audit(current, baseline)
