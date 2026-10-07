"""Bacteria OG concordance against OF3 assigned groups, not official QFO accuracy."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from datetime import datetime
from pathlib import Path

import pyarrow.parquet as pq

from benchmarks.og_extraction.embleya_reference import compare
from ogprofiler.core.manifest import sha256_file


def rows(path):
    with path.open() as h:
        return list(csv.DictReader(h, delimiter="\t"))


def partition(path, id_column, group_column, expected):
    result = {}
    for r in rows(path):
        p = r[id_column]
        if p in result or p not in expected:
            raise ValueError(f"Duplicate/unknown protein in {path}: {p}")
        result[p] = r[group_column]
    if set(result) != expected:
        raise ValueError(f"Incomplete coverage in {path}")
    return result


def of_assigned(path, expected, species):
    result = {}
    with path.open() as h:
        reader = csv.reader(h, delimiter="\t")
        header = next(reader)
        if set(header[1:]) != set(species) or len(header[1:]) != len(species):
            raise ValueError("OF3 species columns differ from frozen inputs")
        groups = set()
        for r in reader:
            if len(r) != len(header) or r[0] in groups:
                raise ValueError("Invalid/duplicate OF3 group")
            groups.add(r[0])
            for cell in r[1:]:
                for p in filter(None, (p.strip() for p in cell.split(","))):
                    if p not in expected or p in result:
                        raise ValueError("Unknown/duplicate OF3 assignment")
                    result[p] = r[0]
    if not result:
        raise ValueError("Empty assigned reference")
    return result


def evaluate(v2, of3, v1, out):
    roots = dict(V2=v2, OF3=of3, V1=v1)
    inputs, artifacts = {}, {}
    for name, root in roots.items():
        completion = json.loads((root / "qfo-methods-completion.json").read_text())
        if (
            not completion["methods_completed"]
            or completion["mode"] != "full"
            or completion["method"] != name.lower()
        ):
            raise ValueError(f"Incomplete/wrong method: {name}")
        inputs[name] = json.loads((root / "campaign/campaign-input.json").read_text())
        info = inputs[name]
        if info["collection"] != "bacteria" or info["mode"] != "full":
            raise ValueError("Only full bacterial collection is supported")
        actual = {p.name: sha256_file(p) for p in (root / "campaign/input").glob("*.fasta")}
        if actual != info["input_sha256"]:
            raise ValueError("Frozen method inputs changed")
    if any(inputs[n]["input_sha256"] != inputs["V2"]["input_sha256"] for n in roots):
        raise ValueError("Methods used different frozen inputs")
    proteins = pq.ParquetFile(v2 / "campaign/v2/input/proteins.parquet").read().to_pylist()
    expected = {r["original_id"] for r in proteins}
    if len(expected) != len(proteins) or len(expected) != inputs["V2"]["n_proteins"]:
        raise ValueError("V2 metadata identities differ")
    mapping = {int(r["protein_id"]): r["original_id"] for r in proteins}
    index = pq.ParquetFile(v2 / "campaign/v2/components/index.parquet").read().to_pylist()
    components = {mapping[int(r["protein_id"])]: int(r["component_id"]) for r in index}
    if len(index) != len(expected) or set(components) != expected:
        raise ValueError("Component index coverage differs")
    predictions = {}
    for name, p, idcol, groupcol in [
        ("V2", v2 / "campaign/v2/results/members.tsv", "original_id", "family_id"),
        ("V1", v1 / "campaign/v1-groups.tsv", "protein_id", "group_id"),
        ("OF3", of3 / "campaign/of3-groups.tsv", "protein_id", "group_id"),
    ]:
        predictions[name] = partition(p, idcol, groupcol, expected)
        artifacts[str(p)] = sha256_file(p)
    mp = rows(v1 / "campaign/v1-id-map.tsv")
    if (
        len(mp) != len(expected)
        or len({r["adapted_id"] for r in mp}) != len(mp)
        or {r["original_id"] for r in mp} != expected
    ):
        raise ValueError("V1 mapping is not bijective")
    primary = Path((of3 / "campaign/of3-assigned-primary.txt").read_text().strip())
    ref = of_assigned(primary, expected, [Path(p).stem for p in inputs["V2"]["input_sha256"]])
    artifacts[str(primary)] = sha256_file(primary)
    unassigned = expected - set(ref)
    metrics, coverage, resources = {}, {}, {}
    for name, pred in predictions.items():
        sizes = Counter(pred.values())
        assigned_groups = {pred[p] for p in ref}
        coverage[name] = dict(
            proteins=len(pred),
            groups=len(sizes),
            singleton_groups=sum(n == 1 for n in sizes.values()),
            largest_group=max(sizes.values()),
            of3_unassigned_proteins=len(unassigned),
            of3_unassigned_attached_to_assigned=sum(pred[p] in assigned_groups for p in unassigned),
        )
        if name != "OF3":
            result = compare(ref, pred, components)
            bp, br = result["bcubed_precision"], result["bcubed_recall"]
            result["bcubed_f1"] = 2 * bp * br / (bp + br) if bp + br else None
            result["component_diagnostic_scope"] = (
                "V2 retained graph components, including for V1 predictions"
            )
            metrics[name] = result
        key = name.lower()
        p = roots[name] / f"campaign/{key}-timing/metadata.tsv"
        resource = {r["key"]: r["value"] for r in rows(p)}
        if resource["exit_code"] != "0":
            raise ValueError("Method timing reports failure")
        resources[name] = resource
        artifacts[str(p)] = sha256_file(p)
    starts = [
        datetime.fromisoformat(r["start_time_utc"].replace("Z", "+00:00"))
        for r in resources.values()
    ]
    ends = [
        datetime.fromisoformat(r["end_time_utc"].replace("Z", "+00:00")) for r in resources.values()
    ]
    report = dict(
        scope=(
            "OG concordance against OF3 assigned partition; "
            "not biological truth or official QFO accuracy"
        ),
        dataset_sha256=inputs["V2"]["dataset_sha256"],
        n_species=inputs["V2"]["n_species"],
        n_proteins=len(expected),
        of3_assigned=len(ref),
        of3_unassigned=len(unassigned),
        v1_id_restoration_verified=True,
        coverage=coverage,
        metrics=metrics,
        resources=resources,
        campaign_method_span_seconds=(max(ends) - min(starts)).total_seconds(),
        source_runs={k: str(v) for k, v in roots.items()},
        artifact_sha256=artifacts,
        official_qfo_scores_produced=False,
    )
    out.mkdir(parents=True, exist_ok=False)
    (out / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    fields = [
        "method",
        "pair_precision",
        "pair_recall",
        "pair_f1",
        "bcubed_precision",
        "bcubed_recall",
        "bcubed_f1",
        "split_reference_groups",
        "merged_predicted_groups",
        "lost_pairs",
        "discordant_pairs",
    ]
    with (out / "metrics.tsv").open("w") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for name, m in metrics.items():
            w.writerow(dict(method=name, **{k: m[k] for k in fields[1:]}))
    print(
        json.dumps(
            {k: v for k, v in report.items() if k not in {"artifact_sha256", "metrics"}}, indent=2
        )
    )
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("v2", "of3", "v1", "out"):
        p.add_argument("--" + key, type=Path, required=True)
    a = p.parse_args()
    evaluate(a.v2, a.of3, a.v1, a.out)


if __name__ == "__main__":
    main()
