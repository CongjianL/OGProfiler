"""P6 scoring: supplied official scorer plus explicitly separate diagnostics."""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
import re
import runpy
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from ogprofiler.core.manifest import sha256_file

OFFICIAL_SCORER_SHA256 = "81eb1e660c17819549b07eea8a54b4fb42a89180cafeb4569d92195d282f5e6f"


def read_groups(path, *, exchange=False):
    result = defaultdict(list)
    with path.open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            result[row["family_id" if exchange else "group_id"]].append(
                row["original_id" if exchange else "protein_id"]
            )
    return dict(result)


def validate_mapping(groups, expected):
    counts = Counter(p for genes in groups.values() for p in genes)
    if not groups or any(not members for members in groups.values()):
        raise ValueError("Empty predictions")
    if any(n != 1 for n in counts.values()):
        raise ValueError("Duplicate assignment in evaluation mapping")
    if not set(counts) <= expected:
        raise ValueError("Unknown original protein IDs in evaluation mapping")
    # Missing predictions remain missing, never fill singleton groups to improve coverage.
    return set(counts)


def refog_diagnostics(groups, truth, uncertain):
    """Raw and confident coverage, best-group diagnostics; not official metrics."""
    assigned = set().union(*map(set, groups.values()))
    rows = []
    for name, raw in truth.items():
        ref = raw - uncertain.get(name, set())
        candidates = []
        for group, genes in groups.items():
            pred = set(genes) - uncertain.get(name, set())
            overlap = len(ref & pred)
            if overlap:
                candidates.append((2 * overlap / (len(ref) + len(pred)), group, pred, overlap))
        if candidates:
            f1, group, pred, overlap = max(candidates, key=lambda x: (x[0], x[1]))
        else:
            f1, group, pred, overlap = 0.0, "", set(), 0
        rows.append(
            dict(
                refog=name,
                confident_genes=len(ref),
                best_group=group,
                best_precision=overlap / len(pred) if pred else 0,
                best_recall=overlap / len(ref) if ref else 0,
                best_f1=f1,
                fragments=len(candidates),
                split_excess=max(0, len(candidates) - 1),
                best_group_extra_genes=len(pred - ref),
                best_extra_examples=";".join(sorted(pred - ref)[:20]),
                fragment_group_examples=";".join(sorted(c[1] for c in candidates)[:20]),
                contaminated_fragments=sum(bool(c[2] - ref) for c in candidates),
                missing_confident=len(ref - assigned),
                covered_raw=len(raw & assigned),
                raw_genes=len(raw),
                classification="MERGE_AND_SPLIT"
                if any(c[2] - ref for c in candidates) and len(candidates) > 1
                else "OVERMERGE"
                if any(c[2] - ref for c in candidates)
                else "SPLIT"
                if len(candidates) > 1
                else "MISSING"
                if ref - assigned
                else "EXACT",
            )
        )
    return rows


def score(groups, benchmark, out, label):
    out.mkdir(parents=True, exist_ok=False)
    scorer = benchmark / "benchmark.py"
    if sha256_file(scorer) != OFFICIAL_SCORER_SHA256:
        raise ValueError("Official scorer differs from preflight source")
    namespace = runpy.run_path(str(scorer))
    expected = namespace["get_expected_genes"]()
    assigned = validate_mapping(groups, expected)
    prediction = out / "prediction.txt"
    prediction.write_text(
        "".join(f"{g}: {' '.join(sorted(m))}\n" for g, m in sorted(groups.items()))
    )
    parsed = namespace["read_orthogroups_smart"](str(prediction))
    # Pin actual default official reader behavior; reject dropped/split IDs and groups.
    if Counter(frozenset(m) for m in parsed) != Counter(frozenset(m) for m in groups.values()):
        raise ValueError("Official smart reader differs from original-ID prediction groups")
    with (
        (out / "official.stdout").open("w") as stdout,
        (out / "official.stderr").open("w") as stderr,
    ):
        subprocess.run(
            [sys.executable, str(scorer), str(prediction)], stdout=stdout, stderr=stderr, check=True
        )
    official_output = (out / "official.stdout").read_text()
    if "ERROR:" in official_output or "% F-score" not in official_output:
        raise ValueError("Official scorer did not produce valid metrics")
    ref_root = benchmark / "RefOGs"
    refs = namespace["read_refogs"](str(ref_root) + "/")
    low = namespace["read_uncertain_refogs"](str(ref_root / "low_certainty_assignments") + "/")
    capture = io.StringIO()
    with contextlib.redirect_stdout(capture):
        f1, precision, recall = namespace["calculate_benchmarks_pairwise"](refs, low, parsed)
    rounded = dict(
        (name, float(value))
        for value, name in re.findall(r"([0-9.]+)% (F-score|Precision|Recall)", official_output)
    )
    for name, value in [("F-score", f1), ("Precision", precision), ("Recall", recall)]:
        if name not in rounded or abs(rounded[name] - value) > 0.051:
            raise ValueError("Official CLI and actual function metrics differ")
    (out / "official-function.stdout").write_text(capture.getvalue())
    truth = {f"RefOG{i:03d}": genes for i, genes in enumerate(refs, 1)}
    uncertain = {f"RefOG{i:03d}": genes for i, genes in enumerate(low, 1)}
    diagnostics = refog_diagnostics(groups, truth, uncertain)
    with (out / "refog-diagnostics.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(diagnostics[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(diagnostics)
    universe = set().union(*refs)
    result = dict(
        label=label,
        official_precision=precision / 100,
        official_recall=recall / 100,
        official_f1=f1 / 100,
        input_proteins=len(expected),
        assigned_proteins=len(assigned),
        coverage=len(assigned) / len(expected),
        refog_raw_coverage=len(assigned & universe) / len(universe),
        groups=len(groups),
        mapping_validated=True,
        official_reader_exact=True,
        official_source_sha256=sha256_file(scorer),
        prediction_sha256=sha256_file(prediction),
        diagnostic_classifications=dict(Counter(r["classification"] for r in diagnostics)),
        diagnostics_scope="best-group confident RefOG metrics; not official pairwise score",
    )
    (out / "metrics.json").write_text(json.dumps(result, indent=2))
    return result


def compare(benchmark, frozen, repaired, v1_groups, v1_manifest, out):
    out.mkdir(parents=True, exist_ok=False)
    results = []
    for label, run in [("historical_frozen", frozen), ("repaired_mean", repaired)]:
        for strategy, members in [
            ("terminal", run / "results/terminal-families/members.tsv"),
            ("v1_compatible", run / "results/members.tsv"),
        ]:
            results.append(
                score(
                    read_groups(members, exchange=True),
                    benchmark,
                    out / f"{label}_{strategy}",
                    f"{label}_{strategy}",
                )
            )
    manifest = dict(
        line.rstrip("\n").split("\t", 1) for line in v1_manifest.read_text().splitlines()
    )
    if (
        manifest["dataset_digest"]
        != "380c85d9548c607f5df6daf656d74b21597fe215ef78d4f8ff296776bb1fdd07"
    ):
        raise ValueError("Full historical V1 dataset digest differs")
    result = score(
        read_groups(v1_groups), benchmark, out / "historical_full_v1", "historical_full_v1"
    )
    result["version_context"] = manifest
    result["groups_source_sha256"] = sha256_file(v1_groups)
    result["comparison_scope"] = (
        "independent full historical first-version pipeline; "
        "different search/SSN/hierarchy/version, not extraction-only control"
    )
    results.append(result)
    deltas = []
    for i, arm in [(0, "historical_frozen"), (2, "repaired_mean")]:
        terminal, og = results[i : i + 2]
        deltas.append(
            dict(
                arm=arm,
                **{
                    key: og[key] - terminal[key]
                    for key in (
                        "official_precision",
                        "official_recall",
                        "official_f1",
                        "coverage",
                        "refog_raw_coverage",
                    )
                },
            )
        )
        by_strategy = []
        for strategy in ("terminal", "v1_compatible"):
            with (out / f"{arm}_{strategy}" / "refog-diagnostics.tsv").open() as stream:
                by_strategy.append({r["refog"]: r for r in csv.DictReader(stream, delimiter="\t")})
        paired = []
        for refog in sorted(by_strategy[0]):
            before, after = (table[refog] for table in by_strategy)
            row = dict(
                refog=refog,
                terminal_class=before["classification"],
                v1_compatible_class=after["classification"],
                attribution="same_hierarchy_selection_effect",
                terminal_best_group=before["best_group"],
                v1_compatible_best_group=after["best_group"],
            )
            row.update(
                {
                    f"delta_{k}": float(after[k]) - float(before[k])
                    for k in (
                        "best_precision",
                        "best_recall",
                        "best_f1",
                        "split_excess",
                        "best_group_extra_genes",
                        "missing_confident",
                    )
                }
            )
            paired.append(row)
        with (out / f"{arm}-paired-refog-deltas.tsv").open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(paired[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(paired)
    report = dict(
        evaluation_completed=True,
        results=results,
        within_hierarchy_deltas=deltas,
        success_criterion="report metrics and signed deltas; no predetermined accuracy improvement",
        attribution=(
            "within-arm effects isolate selection; across arms change upstream SSN/hierarchy"
        ),
    )
    (out / "report.json").write_text(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for name in ("benchmark", "frozen", "repaired", "v1-groups", "v1-manifest", "out"):
        p.add_argument("--" + name, type=Path, required=True)
    args = p.parse_args()
    print(
        json.dumps(
            compare(
                args.benchmark,
                args.frozen,
                args.repaired,
                args.v1_groups,
                args.v1_manifest,
                args.out,
            ),
            indent=2,
        )
    )
