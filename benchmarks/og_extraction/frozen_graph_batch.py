"""Frozen-component graph cuts first; OF evaluation only after predictions are saved."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import pyarrow.parquet as pq

from benchmarks.og_extraction.embleya_reference import evaluate, read_of
from benchmarks.og_extraction.embleya_representation import strata
from benchmarks.og_extraction.fixed_tree_graph import ALGORITHM, graph_cut
from ogprofiler.core.manifest import sha256_file


def run_batch(run, reference_dir, out, penalties):
    if (
        not penalties
        or len(set(penalties)) != len(penalties)
        or any(isinstance(p, bool) or not math.isfinite(p) or p < 0 for p in penalties)
    ):
        raise ValueError("Explicit unique finite nonnegative penalties required")
    out.mkdir(parents=True, exist_ok=False)
    frozen = {}

    def verify(path, expected=None):
        digest = sha256_file(path)
        if expected is not None and digest != expected:
            raise ValueError(f"Checksum mismatch: {path}")
        if path in frozen and frozen[path] != digest:
            raise ValueError(f"Input changed: {path}")
        frozen[path] = digest

    for relative in (
        "input/proteins.parquet",
        "input/species.parquet",
        "components/index.parquet",
        "run.yaml",
    ):
        verify(run / relative)
    for name in ("Orthogroups.tsv", "Orthogroups_UnassignedGenes.tsv"):
        verify(reference_dir / name)  # identity only, not read by scoring
    proteins = pq.read_table(run / "input/proteins.parquet").to_pylist()
    protein = {r["protein_id"]: r for r in proteins}
    rows = pq.read_table(run / "components/index.parquet").to_pylist()
    component = {r["protein_id"]: r["component_id"] for r in rows}
    if (
        len(protein) != len(proteins)
        or len(component) != len(rows)
        or set(protein) != set(component)
    ):
        raise ValueError("Protein/index universe mismatch")
    groups = defaultdict(set)
    for p, cid in component.items():
        groups[cid].add(p)
    predictions = [{} for _ in penalties]
    diagnostics = [[] for _ in penalties]
    protocol = dict(
        algorithm=ALGORITHM,
        penalties=penalties,
        automatic_selection=False,
        production_default_changed=False,
        source_run=str(run),
    )
    (out / "protocol.json").write_text(json.dumps(protocol, indent=2))
    for cid, ids in sorted(groups.items()):
        folder = run / "hierarchy/components" / f"component={cid:08d}"
        edge_dir = run / "components/edges" / f"component={cid:08d}"
        fragments = sorted(edge_dir.glob("*.parquet"))
        if len(ids) == 1 and not (folder / "nodes.parquet").exists():
            if fragments:
                raise ValueError("Unexpected singleton edge partition")
            nodes = [dict(cluster_id=0, parent_id=None, n_genes=1)]
            membership = [(next(iter(ids)), 0)]
            edges = []
        else:
            manifest_path = folder / "hierarchy-manifest.json"
            verify(manifest_path)
            manifest = json.loads(manifest_path.read_text())
            if (
                manifest["hierarchy_status"] != "RESOLVED"
                or not manifest["structural_validation_passed"]
            ):
                raise ValueError(f"Unverified hierarchy: {cid}")
            checks = manifest["input_checksums"]
            recorded_edges = {
                run / p for p in checks if p.startswith(f"components/edges/component={cid:08d}/")
            }
            if recorded_edges != set(fragments) or (len(ids) > 1 and not fragments):
                raise ValueError("Edge partition differs from hierarchy manifest")
            for relative, digest in checks.items():
                path = run / relative
                if Path(relative).is_absolute() or ".." in Path(relative).parts:
                    raise ValueError("Manifest path outside run")
                if path in frozen:
                    if frozen[path] != digest:
                        raise ValueError("Conflicting input checksum")
                else:
                    verify(path, digest)
            for name in ("nodes.parquet", "members.parquet"):
                verify(folder / name, manifest["output_checksums"][name])
            nodes = pq.ParquetFile(folder / "nodes.parquet").read().to_pylist()
            membership = [
                (r["protein_id"], r["terminal_cluster_id"])
                for r in pq.ParquetFile(folder / "members.parquet").read().to_pylist()
            ]
            edges = [
                (r["u"], r["v"], r["weight"])
                for path in fragments
                for r in pq.ParquetFile(path).read(columns=["u", "v", "weight"]).to_pylist()
            ]
        if {p for p, _ in membership} != ids:
            raise ValueError("Hierarchy/index component mismatch")
        for i, penalty in enumerate(penalties):
            cut = graph_cut(nodes, membership, edges, pair_penalty=penalty)
            for row in cut["members"]:
                p = row["protein_id"]
                if p in predictions[i]:
                    raise ValueError("Duplicate prediction")
                predictions[i][p] = f"C{cid}:N{row['cluster_id']}"
            diagnostics[i].append(
                dict(
                    component_id=cid,
                    score=cut["score"],
                    selected=cut["selected"],
                    **cut["diagnostics"],
                )
            )
    # Materialize ALL predictions before reference-based evaluation.
    for i, _penalty in enumerate(penalties):
        if set(predictions[i]) != set(protein):
            raise ValueError("Incomplete prediction")
        folder = out / f"candidate-{i:02d}"
        folder.mkdir()
        with (folder / "members.tsv").open("w") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["family_id", "protein_id", "species_id", "original_id"])
            for p in sorted(protein):
                writer.writerow(
                    [predictions[i][p], p, protein[p]["species_id"], protein[p]["original_id"]]
                )
        (folder / "cuts.json").write_text(json.dumps(diagnostics[i]))
    names = {
        r["species_name"]: r["species_id"]
        for r in pq.read_table(run / "input/species.parquet").to_pylist()
    }
    keys = {(r["species_id"], r["original_id"]): p for p, r in protein.items()}
    ref = {keys[k]: v for k, v in read_of(reference_dir / "Orthogroups.tsv", names, keys).items()}
    species = {p: r["species_id"] for p, r in protein.items()}
    reports = []
    for i, penalty in enumerate(penalties):
        folder = out / f"candidate-{i:02d}"
        result = evaluate(
            run, reference_dir, folder / "evaluation.json", members=folder / "members.tsv"
        )
        reports.append(
            dict(
                pair_penalty=penalty,
                candidate=folder.name,
                primary=result["primary_assigned_only"],
                secondary=result["secondary_unassigned_as_singletons"],
                strata=strata(ref, predictions[i], species),
                nontrivial_root_selected=sum(
                    d["root_selected"] for d in diagnostics[i] if d["proteins"] > 1
                ),
                nontrivial_all_leaves_selected=sum(
                    d["all_structural_leaves_selected"] for d in diagnostics[i] if d["proteins"] > 1
                ),
            )
        )
    for path, digest in frozen.items():
        if sha256_file(path) != digest:
            raise ValueError(f"Frozen input changed: {path}")
    summary = dict(
        completed=True,
        protocol=protocol,
        candidates=reports,
        input_hashes={str(p): h for p, h in frozen.items()},
        source_hashes={
            name: sha256_file(Path(__file__).with_name(name))
            for name in ("frozen_graph_batch.py", "fixed_tree_graph.py", "tree_cut.py")
        },
    )
    (out / "report.json").write_text(json.dumps(summary, indent=2, allow_nan=False))
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--reference-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--penalties", type=float, nargs="+", required=True)
    a = p.parse_args()
    run_batch(a.run, a.reference_dir, a.out, a.penalties)


if __name__ == "__main__":
    main()
