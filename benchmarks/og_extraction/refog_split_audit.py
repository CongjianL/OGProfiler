"""Read-only RefOG first-separation and OG non-rejoining audit on frozen artifacts.

Physical separation in hierarchy is not automatically irreversible: OG selection
can join structural leaves. This report keeps those two questions separate.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from ogprofiler.core.manifest import sha256_file


def prepare_inputs(result_root: Path):
    """Build the small ID-only audit input from retrieved frozen metadata."""
    run = result_root / "new-hierarchy"
    reference = result_root / "benchmark/RefOGs"
    truth = {}
    for index in range(1, 71):
        name = f"RefOG{index:03d}"
        genes = set((reference / f"{name}.txt").read_text().splitlines())
        uncertain = reference / "low_certainty_assignments" / f"{name}.txt"
        truth[name] = genes - (
            set(uncertain.read_text().splitlines()) if uncertain.exists() else set()
        )
    target = set().union(*truth.values())
    protein = {
        p["original_id"]: p
        for p in rows(run / "input/proteins.parquet")
        if p["original_id"] in target
    }
    assert set(protein) == target
    ids = {p["protein_id"] for p in protein.values()}
    components = {
        p["protein_id"]: p["component_id"]
        for p in rows(run / "components/index.parquet")
        if p["protein_id"] in ids
    }
    groups = {}
    with (run / "results/members.tsv").open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["original_id"] in target:
                assert row["original_id"] not in groups
                groups[row["original_id"]] = row["family_id"]
    assert set(groups) == target
    selected = set()
    for genes in truth.values():
        counts = Counter(components[protein[g]["protein_id"]] for g in genes)
        selected.update(cid for cid, count in counts.items() if count > 1)
    data = dict(
        truth={k: sorted(v) for k, v in truth.items()},
        protein=protein,
        components=components,
        groups=groups,
        selected_components=sorted(selected),
    )
    (result_root / "refog-audit-input.json").write_text(json.dumps(data, indent=2) + "\n")
    return data


def rows(path, column=None, values=None):
    table = pq.ParquetFile(path).read()
    if column is not None:
        table = table.filter(pc.is_in(table[column], value_set=pa.array(values)))
    return table.to_pylist()


def load_component(root, component, proteins):
    hierarchy = root / "hierarchy/components" / f"component={component:08d}"
    og = root / "orthogroups/components" / f"component={component:08d}"
    for directory, names in (
        (hierarchy, ("nodes.parquet", "members.parquet", "candidates.parquet", "metrics.json")),
        (og, ("v1_events.parquet", "selection_trace.parquet", "groups.parquet", "members.parquet")),
    ):
        manifest_name = "hierarchy-manifest.json" if directory == hierarchy else "og-manifest.json"
        manifest = json.loads((directory / manifest_name).read_text())
        for name in names:
            if sha256_file(directory / name) != manifest["output_checksums"][name]:
                raise ValueError(f"Frozen artifact checksum mismatch: {directory / name}")
    nodes = {n["cluster_id"]: n for n in rows(hierarchy / "nodes.parquet")}
    membership = {
        m["protein_id"]: m["terminal_cluster_id"]
        for m in rows(hierarchy / "members.parquet", "protein_id", proteins)
    }
    paths = {}
    relevant = set()
    for protein, leaf in membership.items():
        path, current = [], leaf
        while current is not None:
            path.append(current)
            current = nodes[current]["parent_id"]
        paths[protein] = tuple(reversed(path))
        relevant.update(path)
    nodes = {key: nodes[key] for key in relevant}
    events = {
        n["cluster_id"]: n for n in rows(og / "v1_events.parquet", "cluster_id", sorted(relevant))
    }
    traces = defaultdict(list)
    for row in rows(og / "selection_trace.parquet", "cluster_id", sorted(relevant)):
        traces[row["cluster_id"]].append(row)
    candidate_table = pq.ParquetFile(hierarchy / "candidates.parquet").read()
    candidate_table = candidate_table.filter(candidate_table["selected"])
    candidate_table = candidate_table.filter(
        pc.is_in(candidate_table["cluster_id"], value_set=pa.array(sorted(relevant)))
    )
    selected = {c["cluster_id"]: c for c in candidate_table.to_pylist()}
    return dict(
        nodes=nodes,
        membership=membership,
        paths=paths,
        events=events,
        traces=traces,
        selected=selected,
    )


def audit(result_root: Path):
    if not (result_root / "refog-audit-input.json").exists():
        prepare_inputs(result_root)
    data = json.loads((result_root / "refog-audit-input.json").read_text())
    root = result_root / "new-hierarchy"
    protein = data["protein"]
    component = {int(k): v for k, v in data["components"].items()}
    wanted = defaultdict(list)
    for value in protein.values():
        wanted[component[value["protein_id"]]].append(value["protein_id"])
    loaded = {cid: load_component(root, cid, wanted[cid]) for cid in data["selected_components"]}
    summary = Counter()
    weighted = Counter()
    refogs, first_splits = [], {}
    examples = []
    for refog, genes in sorted(data["truth"].items()):
        n = len(genes)
        counts, reasons = Counter(), Counter()
        counts["genes"] = n
        counts["truth_pairs"] = n * (n - 1) // 2
        split_nodes = Counter()
        for left, right in combinations(genes, 2):
            p, q = protein[left]["protein_id"], protein[right]["protein_id"]
            cid, other = component[p], component[q]
            same_og = data["groups"][left] == data["groups"][right]
            counts["retained_pairs" if same_og else "lost_pairs"] += 1
            if cid != other:
                assert not same_og, "OG unexpectedly joins distinct SSN components"
                counts["ssn_lost_pairs"] += 1
                continue
            view = loaded[cid]
            a, b = view["paths"][p], view["paths"][q]
            common = []
            for x, y in zip(a, b, strict=False):
                if x != y:
                    break
                common.append(x)
            assert common
            lca = common[-1]
            if a[-1] == b[-1]:
                counts["same_terminal_pairs"] += 1
                assert same_og, "OG split a structural leaf; investigate stage behavior"
                continue
            counts["hierarchy_separated_pairs"] += 1
            if same_og:
                counts["hierarchy_separated_but_og_rejoined"] += 1
                continue
            counts["hierarchy_lost_pairs"] += 1
            split_nodes[(cid, lca)] += 1
            key = f"{cid}:{lca}"
            if key not in first_splits:
                node = view["nodes"][lca]
                candidate = view["selected"].get(lca, {})
                first_splits[key] = dict(
                    node,
                    v1_event=view["events"][lca]["v1_event"],
                    selected_stability=candidate.get("stability"),
                    lost_pairs=0,
                    refogs=set(),
                )
            first_splits[key]["lost_pairs"] += 1
            first_splits[key]["refogs"].add(refog)
            selectable = [
                c
                for c in reversed(common)
                if view["events"][c]["v1_event"] is None
                or (view["events"][c]["v1_event"] == "I" and view["nodes"][c]["n_species"] > 1)
            ]
            if not selectable:
                reason = "NO_COMMON_SELECTABLE_ANCESTOR"
                candidate = None
                statuses = []
            else:
                candidate = selectable[0]
                statuses = sorted({x["status"] for x in view["traces"][candidate]})
                reason = "COMMON_SELECTABLE_ANCESTOR_" + ("+".join(statuses) or "NO_TRACE")
            reasons[reason] += 1
            counts["lca_event_" + str(view["events"][lca]["v1_event"])] += 1
            if len([e for e in examples if e["refog"] == refog]) < 3:
                examples.append(
                    dict(
                        refog=refog,
                        left=left,
                        right=right,
                        component_id=cid,
                        lca=lca,
                        depth=view["nodes"][lca]["depth"],
                        child_count=view["nodes"][lca]["child_count"],
                        v1_event=view["events"][lca]["v1_event"],
                        reason=reason,
                        closest_selectable_ancestor=candidate,
                        trace_statuses=statuses,
                        left_og=data["groups"][left],
                        right_og=data["groups"][right],
                    )
                )
        assert counts["lost_pairs"] == counts["ssn_lost_pairs"] + counts["hierarchy_lost_pairs"]
        assert (
            counts["retained_pairs"]
            == counts["same_terminal_pairs"] + counts["hierarchy_separated_but_og_rejoined"]
        )
        assert counts["truth_pairs"] == counts["lost_pairs"] + counts["retained_pairs"]
        summary.update({k: v for k, v in counts.items() if k != "genes"})
        summary.update(reasons)
        for key, value in counts.items():
            if (
                key.endswith("pairs")
                or key.startswith("lca_event_")
                or key == "hierarchy_separated_but_og_rejoined"
            ):
                weighted[key] += value / (n - 1)
        refogs.append(
            dict(
                refog=refog,
                **counts,
                ssn_components=len({component[protein[g]["protein_id"]] for g in genes}),
                predicted_fragments=len({data["groups"][g] for g in genes}),
                pair_recall=counts["retained_pairs"] / counts["truth_pairs"],
                non_rejoining_reasons=dict(reasons),
                first_split_nodes=[
                    dict(component_id=cid, cluster_id=cluster, lost_pairs=count)
                    for (cid, cluster), count in split_nodes.most_common()
                ],
            )
        )
    first_splits = [{**v, "refogs": sorted(v["refogs"])} for v in first_splits.values()]
    report = dict(
        read_only=True,
        source_job=(result_root.parent / "job_id.txt").read_text().strip()
        if (result_root.parent / "job_id.txt").is_file()
        else None,
        refogs=len(refogs),
        confident_genes=len(protein),
        selected_components=len(loaded),
        counts=dict(summary),
        official_even_weighted_pairs=dict(weighted),
        reconstructed_official_recall=weighted["retained_pairs"] / weighted["truth_pairs"],
        ssn_connected_pair_recall_ceiling=1 - weighted["ssn_lost_pairs"] / weighted["truth_pairs"],
        hierarchy_terminal_pair_recall=weighted["same_terminal_pairs"] / weighted["truth_pairs"],
        input_sha256=sha256_file(result_root / "refog-audit-input.json"),
        audit_source_sha256=sha256_file(Path(__file__)),
        interpretation=(
            "First physical separation and final non-rejoining are separate; "
            "no counterfactual causal proof."
        ),
    )
    expected = json.loads((result_root / "h5-report-compact.json").read_text())["report"][
        "results"
    ][3]["official_recall"]
    assert abs(report["reconstructed_official_recall"] - expected) < 1e-12
    return report, refogs, sorted(first_splits, key=lambda x: -x["lost_pairs"]), examples


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    report, refs, splits, examples = audit(args.results)
    for name, value in [
        ("summary", report),
        ("refogs", refs),
        ("first-split-nodes", splits),
        ("examples", examples),
    ]:
        (args.out / f"{name}.json").write_text(json.dumps(value, indent=2) + "\n")
    with (args.out / "refogs.tsv").open("w") as stream:
        keys = [
            "refog",
            "genes",
            "ssn_components",
            "predicted_fragments",
            "pair_recall",
            "truth_pairs",
            "ssn_lost_pairs",
            "hierarchy_lost_pairs",
            "same_terminal_pairs",
            "hierarchy_separated_but_og_rejoined",
        ]
        writer = csv.DictWriter(
            stream, fieldnames=keys, delimiter="\t", extrasaction="ignore", restval=0
        )
        writer.writeheader()
        writer.writerows(refs)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
