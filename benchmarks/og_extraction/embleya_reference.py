"""Partition concordance with frozen OrthoFinder OGs, not biological ground truth.

Linear-sized contingency tables; never enumerate all protein pairs. OF unassigned
proteins are excluded from the primary score and treated as distinct singletons
only in the explicitly named secondary score.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import pyarrow.parquet as pq

from ogprofiler.core.manifest import sha256_file


def choose2(n):
    return n * (n - 1) // 2


def compare(reference, prediction, components):
    universe = set(reference)
    extra = set(prediction) - universe
    # Missing assignments are distinct singleton predictions, never one fake OG.
    pred = {g: ("assigned", prediction[g]) if g in prediction else ("missing", g) for g in universe}
    cells = Counter((reference[g], pred[g]) for g in universe)
    nr, np = Counter(reference.values()), Counter(pred.values())
    tp = sum(choose2(n) for n in cells.values())
    rp, pp = sum(map(choose2, nr.values())), sum(map(choose2, np.values()))
    precision = tp / pp if pp else None
    recall = tp / rp if rp else None
    fragments, merged = defaultdict(dict), defaultdict(dict)
    for (r, p), n in cells.items():
        fragments[r][p] = n
        merged[p][r] = n
    details = []
    for r, sizes in fragments.items():
        retained = sum(map(choose2, sizes.values()))
        best = max(2 * n / (nr[r] + np[p]) for p, n in sizes.items())
        details.append(
            dict(
                reference_og=r,
                size=nr[r],
                fragments=len(sizes),
                lost_pairs=choose2(nr[r]) - retained,
                best_f1=best,
            )
        )
    contamination = []
    for p, sizes in merged.items():
        contamination.append(
            dict(
                predicted_og=str(p[1]),
                size=np[p],
                reference_ogs=len(sizes),
                discordant_pairs=choose2(np[p]) - sum(map(choose2, sizes.values())),
            )
        )
    cc = Counter((reference[g], components[g]) for g in universe)
    supported = sum(map(choose2, cc.values()))
    # Evaluate whether predictions unexpectedly cross component boundaries too.
    pc = Counter((pred[g], components[g]) for g in universe)
    cross_component_predictions = sum(map(choose2, np.values())) - sum(map(choose2, pc.values()))
    return dict(
        proteins=len(universe),
        reference_groups=len(nr),
        predicted_groups=len(np),
        missing_predictions=len(universe - set(prediction)),
        excluded_predictions=len(extra),
        true_positive_pairs=tp,
        reference_pairs=rp,
        predicted_pairs=pp,
        discordant_pairs=pp - tp,
        lost_pairs=rp - tp,
        pair_precision=precision,
        pair_recall=recall,
        pair_f1=2 * tp / (pp + rp) if pp + rp else None,
        bcubed_precision=sum(n * n / np[p] for (r, p), n in cells.items()) / len(universe),
        bcubed_recall=sum(n * n / nr[r] for (r, p), n in cells.items()) / len(universe),
        macro_best_f1=sum(d["best_f1"] for d in details) / len(details),
        exact_reference_groups=sum(
            len(fragments[r]) == 1 and np[next(iter(fragments[r]))] == n for r, n in nr.items()
        ),
        split_reference_groups=sum(len(v) > 1 for v in fragments.values()),
        merged_predicted_groups=sum(len(v) > 1 for v in merged.values()),
        component_pair_recall_ceiling=supported / rp if rp else None,
        reference_pairs_across_components=rp - supported,
        cross_component_predicted_pairs=cross_component_predictions,
        top_splits=sorted(details, key=lambda d: d["lost_pairs"], reverse=True)[:30],
        top_merges=sorted(contamination, key=lambda d: d["discordant_pairs"], reverse=True)[:30],
    )


def read_of(path, species, keys):
    groups = {}
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        if set(reader.fieldnames[1:]) != set(species):
            raise ValueError("OrthoFinder and V2 species differ")
        for row in reader:
            for name in reader.fieldnames[1:]:
                for gene in filter(None, (g.strip() for g in row[name].split(","))):
                    key = (species[name], gene)
                    if key not in keys or key in groups:
                        raise ValueError(f"Unknown or duplicate OF protein: {key}")
                    groups[key] = row[reader.fieldnames[0]]
    return groups


def evaluate(run, reference_dir, out, *, members=None):
    species = {
        r["species_name"]: r["species_id"]
        for r in pq.read_table(run / "input/species.parquet").to_pylist()
    }
    proteins = pq.read_table(run / "input/proteins.parquet").to_pylist()
    keys = {(r["species_id"], r["original_id"]): r["protein_id"] for r in proteins}
    if len(keys) != len(proteins):
        raise ValueError("Duplicate species/protein identity")
    index = {
        r["protein_id"]: r["component_id"]
        for r in pq.read_table(run / "components/index.parquet").to_pylist()
    }
    components = {g: index[p] for g, p in keys.items()}
    assigned = read_of(reference_dir / "Orthogroups.tsv", species, keys)
    unassigned = read_of(reference_dir / "Orthogroups_UnassignedGenes.tsv", species, keys)
    if set(assigned) & set(unassigned) or set(assigned) | set(unassigned) != set(keys):
        raise ValueError("Reference does not cover the identical protein universe")
    prediction = {}
    members = members or run / "results/members.tsv"
    with members.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            key = (int(r["species_id"]), r["original_id"])
            if key not in keys or key in prediction or keys[key] != int(r["protein_id"]):
                raise ValueError(f"Unknown, duplicate or inconsistent prediction: {key}")
            prediction[key] = r["family_id"]
    all_reference = {g: ("og", r) for g, r in assigned.items()}
    all_reference.update({g: ("unassigned", g) for g in unassigned})
    paths = [
        run / "input/proteins.parquet",
        run / "input/species.parquet",
        run / "components/index.parquet",
        run / "run.yaml",
        members,
        reference_dir / "Orthogroups.tsv",
        reference_dir / "Orthogroups_UnassignedGenes.tsv",
    ]
    report = dict(
        reference_kind="OrthoFinder partition, not biological truth",
        run=str(run),
        species=len(species),
        input_proteins=len(keys),
        of_unassigned=len(unassigned),
        prediction_assigned=len(prediction),
        hashes={str(p): sha256_file(p) for p in paths},
        primary_assigned_only=compare(assigned, prediction, components),
        secondary_unassigned_as_singletons=compare(all_reference, prediction, components),
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--members", type=Path)
    args = parser.parse_args()
    evaluate(args.run, args.reference_dir, args.out, members=args.members)


if __name__ == "__main__":
    main()
