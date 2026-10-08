"""Frozen V1 hierarchy algorithm on selected components of an immutable mean SSN.

No V2 search or gates. V1's seed-free stochastic behavior is preserved, so this
is a recorded realization, not a claim of historical bitwise replay.
"""

from __future__ import annotations

import argparse
import ast
import json
import platform
import time
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from types import SimpleNamespace

import igraph
import leidenalg
import numpy as np

from benchmarks.og_extraction.reference_v1 import REFERENCE_PATH, REFERENCE_SHA256, run_reference
from benchmarks.og_extraction.refog_split_audit import load_component, rows
from ogprofiler.core.manifest import sha256_file
from ogprofiler.graph.partition import load_component_edge_table

TARGETS = ("RefOG012", "RefOG015", "RefOG019", "RefOG036")


def load_hierarchy_reference(find_partition):
    assert sha256_file(REFERENCE_PATH) == REFERENCE_SHA256
    tree = ast.parse(REFERENCE_PATH.read_text())
    methods = {"__init__", "SpecificBoard", "BipartiteGraphs", "RunCommunityDetection"}
    definitions = []
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == "HHN":
            definitions.append(
                ast.ClassDef(
                    name=node.name,
                    bases=node.bases,
                    keywords=[],
                    body=[
                        m for m in node.body if isinstance(m, ast.FunctionDef) and m.name in methods
                    ],
                    decorator_list=[],
                )
            )
        elif isinstance(node, ast.FunctionDef) and node.name in {
            "GetAttribution",
            "MappingGammaForCC",
        }:
            definitions.append(node)
    env = dict(
        __name__="frozen_v1_hierarchy",
        igraph=igraph,
        np=np,
        la=SimpleNamespace(
            find_partition=find_partition,
            ModularityVertexPartition=leidenalg.ModularityVertexPartition,
            RBConfigurationVertexPartition=leidenalg.RBConfigurationVertexPartition,
            RBERVertexPartition=leidenalg.RBERVertexPartition,
            CPMVertexPartition=leidenalg.CPMVertexPartition,
        ),
    )
    exec(
        compile(
            ast.fix_missing_locations(ast.Module(body=definitions, type_ignores=[])),
            str(REFERENCE_PATH),
            "exec",
        ),
        env,
    )
    return env


def original_tree(graph, ids, species, *, trace_path=None, find_partition=None):
    """BFS orchestrator around unmodified V1 RunCommunityDetection definitions."""
    graph = graph.copy()
    graph.vs["name"] = [f"{species[p]}|g{p}" for p in ids]
    graph.vs["index"] = graph.vs.indices
    index_of = {label: i for i, label in enumerate(graph.vs["name"])}
    records, calls, owner = [], [], {}
    delegate = find_partition or leidenalg.find_partition

    def observe(subgraph, *args, **kwargs):
        # Forward exactly the native arguments: no seed, robustness or gate insertion.
        assert "seed" not in kwargs
        started = time.perf_counter()
        partition = delegate(subgraph, *args, **kwargs)
        record = dict(
            call_index=len(calls),
            node=owner[frozenset(subgraph.vs["name"])],
            gamma=kwargs.get("resolution_parameter"),
            n_iterations=kwargs["n_iterations"],
            seed=None,
            membership=list(partition.membership),
            vertex_names=subgraph.vs["name"],
            child_count=len(partition.subgraphs()),
            quality=float(partition.q),
            seconds=time.perf_counter() - started,
        )
        calls.append(record)
        if trace_path:
            with trace_path.open("a") as stream:
                stream.write(json.dumps(record) + "\n")
        return partition

    env = load_hierarchy_reference(observe)
    obj = env["HHN"](graph)
    frontier, depth = [("0", list(range(len(ids))))], 0
    while frontier:
        owner = {
            frozenset(graph.vs[i]["name"] for i in indices): name for name, indices in frontier
        }
        connected = []
        following = obj.RunCommunityDetection(connected, frontier, "rber", "weight", 1.0)
        current = dict(frontier)
        for line in connected:
            name, child_names, gene_ids, n_genes, species_ids, n_species, quality = line.split(",")
            children = [] if child_names == "+" else child_names.split()
            proteins = [ids[index_of[g]] for g in gene_ids.split()]
            assert set(proteins) == {ids[i] for i in current[name]}
            terminal = (
                (
                    "SINGLETON"
                    if len(proteins) == 1
                    else "ONE_SPECIES"
                    if int(n_species) == 1
                    else "V1_SEARCH_NOT_EXACTLY_TWO"
                )
                if not children
                else None
            )
            records.append(
                dict(
                    name=name,
                    parent=name.rsplit("-", 1)[0] if "-" in name else None,
                    depth=depth,
                    members=sorted(proteins),
                    species=species_ids.split(),
                    children=children,
                    terminal_reason=terminal,
                    quality=quality,
                )
            )
        frontier, depth = following, depth + 1
    return records, calls


def paths_and_clades(nodes):
    by_name = {n["name"]: n for n in nodes}
    paths = {}
    for node in nodes:
        if node["children"]:
            child_members = [set(by_name[c]["members"]) for c in node["children"]]
            assert (
                sum(map(len, child_members))
                == len(set().union(*child_members))
                == len(node["members"])
            )
            assert set().union(*child_members) == set(node["members"])
            assert all(len(m) < len(node["members"]) for m in child_members)
            continue
        path, current = [], node["name"]
        while current is not None:
            path.append(current)
            current = by_name[current]["parent"]
        for protein in node["members"]:
            assert protein not in paths
            paths[protein] = list(reversed(path))
    return paths, {n["name"]: frozenset(n["members"]) for n in nodes}


def compare_component(run, cid, metadata, out):
    folder = run / f"hierarchy/components/component={cid:08d}"
    manifest = json.loads((folder / "hierarchy-manifest.json").read_text())
    checks = manifest["input_checksums"]
    for relative, digest in checks.items():
        assert sha256_file(run / relative) == digest
    table = load_component_edge_table(run / "components", cid)
    graph, ids = table.to_igraph()
    assert not graph.is_directed()
    species = {p: r["species_id"] for p, r in metadata.items()}
    originals = {p: r["original_id"] for p, r in metadata.items()}
    nodes, calls = original_tree(
        graph, ids, species, trace_path=out / f"component-{cid}-calls.jsonl"
    )
    paths, clades = paths_and_clades(nodes)
    assert set(paths) == set(ids)
    labels = {p: f"{species[p]}|g{p}" for p in ids}
    reverse = {g: p for p, g in labels.items()}
    spec = dict(
        n_species=len(set(species.values())),
        isolates=[],
        vertices=[
            dict(name=n["name"], genes=[labels[p] for p in n["members"]], species=n["species"])
            for n in nodes
        ],
        edges=[(n["parent"], n["name"]) for n in nodes if n["parent"] is not None],
    )
    og = run_reference(spec)
    groups = [{**g, "members": [reverse[p] for p in g["members"]]} for g in og["groups"]]
    occurrences = Counter(p for g in groups for p in g["members"])
    assert all(v == 1 for v in occurrences.values())
    v2 = load_component(run, cid, ids)
    v2clades = defaultdict(set)
    for p, path in v2["paths"].items():
        for node in path:
            v2clades[node].add(p)
    children = defaultdict(list)
    for node in v2["nodes"].values():
        if node["parent_id"] is not None:
            children[node["parent_id"]].append(node["cluster_id"])
    v1nodes = {n["name"]: n for n in nodes}
    divergence = {}
    for p in ids:
        for a, b in zip(paths[p], v2["paths"][p], strict=False):
            assert clades[a] == v2clades[b]
            left = {clades[c] for c in v1nodes[a]["children"]}
            right = {frozenset(v2clades[c]) for c in children[b]}
            if left != right:
                divergence[str(p)] = dict(
                    v1_node=a,
                    v1_depth=v1nodes[a]["depth"],
                    v1_event=og["raw_events"][a],
                    v1_children=len(left),
                    v1_members=sorted(clades[a]),
                    v2_node=b,
                    v2_depth=v2["nodes"][b]["depth"],
                    v2_event=v2["events"][b]["v1_event"],
                    v2_children=len(right),
                    v2_selection_kind=v2["nodes"][b].get("selection_kind"),
                )
                break
    report = dict(
        component_id=cid,
        proteins=len(ids),
        edges=graph.ecount(),
        nodes=nodes,
        leiden_calls=len(calls),
        raw_events=og["raw_events"],
        groups=groups,
        unassigned=[reverse[g] for g in og["unassigned"]],
        first_partition_divergence=divergence,
        input_checksums=checks,
        originals={str(p): originals[p] for p in ids},
    )
    for relative, digest in checks.items():
        assert sha256_file(run / relative) == digest
    (out / f"component-{cid}.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--origin", type=Path, required=True, help="Frozen H5 job root")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    run = args.origin / "new-hierarchy"
    edge_manifest = json.loads((run / "edges/edge-manifest.json").read_text())
    assert edge_manifest["parameters"]["symmetrization"] == "mean"
    meta = {p["protein_id"]: p for p in rows(run / "input/proteins.parquet")}
    original = {p["original_id"]: p["protein_id"] for p in meta.values()}
    components = {
        r["protein_id"]: r["component_id"] for r in rows(run / "components/index.parquet")
    }
    truth = {}
    uncertain = {}
    benchmark_checksums = {}
    for name in TARGETS:
        raw = set((args.origin / f"benchmark/RefOGs/{name}.txt").read_text().splitlines())
        low = args.origin / f"benchmark/RefOGs/low_certainty_assignments/{name}.txt"
        uncertain[name] = set(low.read_text().splitlines()) if low.exists() else set()
        truth[name] = raw - uncertain[name]
        benchmark_checksums[name] = dict(
            truth=sha256_file(args.origin / f"benchmark/RefOGs/{name}.txt"),
            low=sha256_file(low) if low.exists() else None,
        )
    wanted = {components[original[g]] for genes in truth.values() for g in genes}
    stats = {
        r["component_id"]: r["n_vertices"] for r in rows(run / "components/statistics.parquet")
    }
    reports = {
        cid: compare_component(run, cid, meta, args.out) for cid in sorted(wanted) if stats[cid] > 1
    }
    predicted = {
        p: (cid, i)
        for cid, r in reports.items()
        for i, g in enumerate(r["groups"])
        for p in g["members"]
    }
    for cid in wanted:
        if stats[cid] == 1:
            for p in original.values():
                if components[p] == cid:
                    predicted[p] = (cid, "ISOLATED_SSN_V1_LEVEL0")
    from benchmarks.og_extraction.refog_pollution_audit import export_members

    v2sets, v2groups = export_members(args.origin)
    v1paths = {cid: paths_and_clades(r["nodes"])[0] for cid, r in reports.items()}
    v2views = {cid: load_component(run, cid, list(v1paths[cid])) for cid in reports}
    from benchmarks.og_extraction.refog_pollution_audit import common_path

    results = []
    for name, genes in truth.items():
        retained_v1 = retained_v2 = ssn_lost = 0
        pair_records = []
        for a, b in combinations(sorted(genes), 2):
            p, q = original[a], original[b]
            retained_v1 += p in predicted and q in predicted and predicted[p] == predicted[q]
            retained_v2 += v2groups[a] == v2groups[b]
            ssn_lost += components[p] != components[q]
            record = dict(
                left=a,
                right=b,
                v1_retained=p in predicted and q in predicted and predicted[p] == predicted[q],
                v2_retained=v2groups[a] == v2groups[b],
            )
            if components[p] != components[q]:
                record["first_separation_stage"] = "SSN"
            else:
                cid = components[p]
                lca1 = common_path(v1paths[cid][p], v1paths[cid][q])[-1]
                lca2 = common_path(v2views[cid]["paths"][p], v2views[cid]["paths"][q])[-1]
                record.update(
                    component_id=cid,
                    v1_lca=lca1,
                    v1_lca_event=reports[cid]["raw_events"][lca1],
                    v2_lca=lca2,
                    v2_lca_event=v2views[cid]["events"][lca2]["v1_event"],
                    v1_same_terminal=v1paths[cid][p][-1] == v1paths[cid][q][-1],
                    v2_same_terminal=v2views[cid]["paths"][p][-1] == v2views[cid]["paths"][q][-1],
                )
            pair_records.append(record)
        (args.out / f"{name}-pairs.json").write_text(json.dumps(pair_records, indent=2) + "\n")
        target_ids = {original[g] for g in genes}
        low_ids = {original[g] for g in uncertain[name]}

        def best_diagnostic(sets, low_ids=low_ids, target_ids=target_ids):
            candidates = []
            for label, members in sets.items():
                confident = set(members) - low_ids
                overlap = len(confident & target_ids)
                if overlap:
                    candidates.append(
                        (
                            2 * overlap / (len(target_ids) + len(confident)),
                            str(label),
                            confident,
                            overlap,
                        )
                    )
            if not candidates:
                return None
            f1, label, members, overlap = max(candidates, key=lambda c: (c[0], c[1]))
            return dict(
                group=label,
                best_f1=f1,
                best_recall=overlap / len(target_ids),
                best_precision=overlap / len(members),
                extras=sorted(meta[p]["original_id"] for p in members - target_ids),
            )

        v1sets = {
            str((cid, i)): g["members"]
            for cid, r in reports.items()
            for i, g in enumerate(r["groups"])
        }
        for pid, label in predicted.items():
            if isinstance(label[1], str):
                v1sets[str(label)] = [pid]
        filtered_v2sets = {
            label: [original[g] for g in members]
            for label, members in v2sets.items()
            if members & genes
        }
        results.append(
            dict(
                refog=name,
                confident_genes=len(genes),
                truth_pairs=len(genes) * (len(genes) - 1) // 2,
                v1_retained_pairs=retained_v1,
                v2_retained_pairs=retained_v2,
                ssn_lost_pairs=ssn_lost,
                v1_best_group=best_diagnostic(v1sets),
                v2_best_group=best_diagnostic(filtered_v2sets),
                v1_missing=sorted(g for g in genes if original[g] not in predicted),
                components=sorted({components[original[g]] for g in genes}),
            )
        )
    report = dict(
        completed=True,
        scope="four targeted RefOG complete components; not full-dataset official benchmark",
        source_job="1410868",
        frozen_ssn=str(run),
        v1_source_sha256=REFERENCE_SHA256,
        audit_source_sha256=sha256_file(Path(__file__)),
        benchmark_checksums=benchmark_checksums,
        v1_seed=None,
        stochastic_realizations_per_component=1,
        original_v1_methods_unmodified=True,
        V2_parameters_unchanged=True,
        environment=dict(
            python=platform.python_version(),
            igraph=igraph.__version__,
            leidenalg=leidenalg.__version__,
            numpy=np.__version__,
        ),
        non_singleton_components=sorted(reports),
        singleton_components=sorted(wanted - set(reports)),
        refogs=results,
    )
    (args.out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
