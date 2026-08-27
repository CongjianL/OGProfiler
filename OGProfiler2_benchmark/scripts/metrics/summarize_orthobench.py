#!/usr/bin/env python3
"""Create complete RefOG summary, split diagnostics, and two-panel QC figure."""
from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


def read_groups(path: Path) -> dict[str, set[str]]:
    result: dict[str, set[str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result.setdefault(row["group_id"], set()).add(row["protein_id"])
    return result


def read_metrics(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def median(values: list[float]) -> float:
    return float(statistics.median(values))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--groups", required=True, type=Path)
    parser.add_argument("--refogs", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--diagnostics", required=True, type=Path)
    parser.add_argument("--figure-prefix", required=True, type=Path)
    parser.add_argument("--vi", required=True, type=float)
    args = parser.parse_args()

    rows = read_metrics(args.metrics)
    if len(rows) != 70:
        raise ValueError(f"Expected 70 RefOG rows, observed {len(rows)}")
    groups = read_groups(args.groups)
    truth = {
        path.stem: {line.strip() for line in path.read_text().splitlines() if line.strip()}
        for path in sorted(args.refogs.glob("RefOG*.txt"))
    }

    diagnostics: list[dict[str, object]] = []
    for row in rows:
        refog = row["refog"]
        reference = truth[refog]
        fragments = sorted(
            (len(reference & members) for members in groups.values() if reference & members),
            reverse=True,
        )
        best = row["best_predicted_group"]
        diagnostics.append(
            {
                "RefOG": refog,
                "true_size": len(reference),
                "n_intersecting_predicted_groups": len(fragments),
                "best_group_size": len(groups[best]) if best else 0,
                "best_group_intersection": int(row["intersection"]),
                "F1": float(row["F1"]),
                "split_count": int(row["split_count"]),
                "largest_fragment_fraction": fragments[0] / len(reference) if fragments else 0.0,
                "second_largest_fragment_fraction": fragments[1] / len(reference) if len(fragments) > 1 else 0.0,
            }
        )

    f1 = [float(row["F1"]) for row in rows]
    splits = [float(row["split_count"]) for row in rows]
    contamination = [float(row["contamination"]) for row in rows]
    missing = [float(row["missing_fraction"]) for row in rows]
    exact = sum(row["exact"].lower() == "true" for row in rows)
    summary = [
        ("macro_best_group_F1", statistics.fmean(f1)),
        ("median_refog_F1", median(f1)),
        ("exact_recovery_n", exact),
        ("exact_recovery_fraction", exact / len(rows)),
        ("mean_split", statistics.fmean(splits)),
        ("median_split", median(splits)),
        ("fraction_refogs_split", sum(value > 0 for value in splits) / len(rows)),
        ("mean_contamination", statistics.fmean(contamination)),
        ("median_contamination", median(contamination)),
        ("mean_missing", statistics.fmean(missing)),
        ("median_missing", median(missing)),
        ("VI_refog_universe", args.vi),
        ("ARI", "DEFERRED"),
    ]
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    with args.summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows(summary)
    with args.diagnostics.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(diagnostics[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(diagnostics)

    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.size": 7,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "savefig.dpi": 300,
    })
    ordered = sorted(diagnostics, key=lambda row: (float(row["F1"]), str(row["RefOG"])))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), constrained_layout=True)
    axes[0].plot(range(len(ordered)), [float(row["F1"]) for row in ordered], "o-", color="#2563A6", markersize=2.5, linewidth=0.8)
    axes[0].set(title="Best-group recovery varies across 70 RefOGs", xlabel="RefOG (ordered by F1)", ylabel="Best-group F1")
    ticks = list(range(0, len(ordered), 7))
    axes[0].set_xticks(ticks, [str(ordered[index]["RefOG"]) for index in ticks], rotation=45, ha="right")
    axes[0].set_ylim(-0.03, 1.03)
    axes[0].text(0.02, 0.04, "higher = better", transform=axes[0].transAxes, color="#2563A6")
    axes[1].scatter([int(row["true_size"]) for row in diagnostics], [int(row["split_count"]) for row in diagnostics], s=16, color="#D97706", alpha=0.85, edgecolors="white", linewidths=0.3)
    axes[1].set(title="Fragmentation increases for some RefOGs", xlabel="True RefOG size", ylabel="Split count")
    axes[1].margins(0.06)
    args.figure_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.figure_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(args.figure_prefix.with_suffix(".pdf"))
    plt.close(fig)


if __name__ == "__main__":
    main()
