#!/usr/bin/env python3
"""Build frozen Stage 2A comparison tables, paired statistics, and report."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy import stats


METHODS = ["OGProfiler2", "OrthoFinder3", "FastOMA", "SonicParanoid2", "Proteinortho6"]
VERSIONS = {
    "OGProfiler2": "2.0.0a1",
    "OrthoFinder3": "3.1.5",
    "FastOMA": "0.5.1",
    "SonicParanoid2": "2.0.9",
    "Proteinortho6": "6.3.6",
}


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def kv(path: Path) -> dict[str, str]:
    return {r["metric"] if "metric" in r else r["key"]: r["value"] for r in rows(path)}


def write(path: Path, records: list[dict[str, object]], fields: list[str] | None = None) -> None:
    if not records:
        raise ValueError(f"no rows for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def f(values: list[dict[str, str]], field: str) -> np.ndarray:
    return np.asarray([float(r[field]) for r in values], dtype=float)


def wall_seconds(value: str) -> float:
    parts = [float(x) for x in value.split(":")]
    if len(parts) == 2:
        return parts[0] * 60 + parts[1]
    if len(parts) == 3:
        return parts[0] * 3600 + parts[1] * 60 + parts[2]
    raise ValueError(f"unexpected wall clock: {value}")


def holm(pvalues: list[float]) -> list[float]:
    order = np.argsort(pvalues)
    adjusted = np.zeros(len(pvalues), dtype=float)
    running = 0.0
    m = len(pvalues)
    for rank, index in enumerate(order):
        running = max(running, (m - rank) * pvalues[index])
        adjusted[index] = min(1.0, running)
    return adjusted.tolist()


def rank_biserial(diff: np.ndarray) -> float:
    nonzero = diff[diff != 0]
    if not len(nonzero):
        return 0.0
    ranks = stats.rankdata(np.abs(nonzero), method="average")
    positive = float(ranks[nonzero > 0].sum())
    negative = float(ranks[nonzero < 0].sum())
    return (positive - negative) / (positive + negative)


def percentile_ci(values: np.ndarray) -> tuple[float, float]:
    low, high = np.quantile(values, [0.025, 0.975])
    return float(low), float(high)


def md_table(records: list[dict[str, object]], fields: list[str]) -> str:
    out = ["| " + " | ".join(fields) + " |", "|" + "|".join(["---"] * len(fields)) + "|"]
    for record in records:
        out.append("| " + " | ".join(str(record[x]) for x in fields) + " |")
    return "\n".join(out)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--metrics-out", type=Path, required=True)
    parser.add_argument("--figure-data", type=Path, required=True)
    parser.add_argument("--bootstrap-seed", type=int, default=20260901)
    parser.add_argument("--bootstrap-n", type=int, default=10000)
    args = parser.parse_args()

    method_rows: dict[str, list[dict[str, str]]] = {}
    official_summary: list[dict[str, object]] = []
    extended_summary: list[dict[str, object]] = []
    error_summary: list[dict[str, object]] = []
    resource_summary: list[dict[str, object]] = []

    for method in METHODS:
        source = args.input / method
        official = kv(source / "official.tsv")
        refs = rows(source / "refog_metrics.tsv")
        if len(refs) != 70:
            raise ValueError(f"{method}: expected 70 RefOG rows, observed {len(refs)}")
        if [r["refog"] for r in refs] != sorted(r["refog"] for r in refs):
            raise ValueError(f"{method}: RefOG rows are not sorted")
        method_rows[method] = refs
        official_summary.append({
            "method": method, "version": VERSIONS[method],
            "precision": official["precision"], "recall": official["recall"],
            "F_score": official["f_score"], "official_exact": official["orthogroups_exactly_correct"],
            "status": "PASS",
        })
        split = f(refs, "split_count")
        contamination = f(refs, "contamination")
        missing = f(refs, "missing_fraction")
        f1 = f(refs, "F1")
        exact = np.asarray([r["exact"].lower() == "true" for r in refs], dtype=bool)
        extended_summary.append({
            "method": method,
            "macro_best_group_F1": float(np.mean(f1)),
            "median_refog_F1": float(np.median(f1)),
            "strict_exact_n": int(exact.sum()),
            "strict_exact_fraction": float(exact.mean()),
            "fraction_refogs_split": float(np.mean(split > 0)),
            "mean_split": float(np.mean(split)), "median_split": float(np.median(split)),
            "mean_contamination": float(np.mean(contamination)), "median_contamination": float(np.median(contamination)),
            "mean_missing": float(np.mean(missing)), "median_missing": float(np.median(missing)),
            "VI": float(kv(source / "extended_summary_raw.tsv")["variation_of_information_refog_universe"]),
        })
        fp = kv(source / "pairwise_fp_concentration.tsv")
        family = kv(source / "family_size_summary.tsv")
        error_summary.append({
            "method": method,
            "top1_FP_fraction": fp["top_1_family_FP_fraction"],
            "top5_FP_fraction": fp["top_5_family_FP_fraction"],
            "top10_FP_fraction": fp["top_10_family_FP_fraction"],
            "largest_predicted_family": family["largest_family_size"],
            "n_families_gt1000": family["n_size_gt1000"],
            "singleton_fraction": family["fraction_singletons"],
        })
        metadata = kv(source / "metadata.tsv")
        wall = wall_seconds(metadata["wall_clock"])
        cpu = float(metadata["user_cpu_seconds"]) + float(metadata["system_cpu_seconds"])
        resource_summary.append({
            "method": method, "wall_seconds": wall, "total_cpu_seconds": cpu,
            "effective_mean_cores": cpu / wall,
            "peak_rss_gib": float(metadata["peak_rss_kb"]) / 1024**2,
            "disk_gib": float(metadata["run_disk_bytes"]) / 1024**3,
            "status": metadata["status"],
        })

    args.metrics_out.mkdir(parents=True, exist_ok=True)
    write(args.metrics_out / "B1_official_summary.tsv", official_summary)
    write(args.metrics_out / "B1_extended_summary.tsv", extended_summary)
    write(args.metrics_out / "B1_error_structure_summary.tsv", error_summary)
    write(args.metrics_out / "B1_resource_summary.tsv", resource_summary)

    matrix_specs = {
        "B1_refog_F1_matrix.tsv": "F1",
        "B1_refog_split_matrix.tsv": "split_count",
        "B1_refog_contamination_matrix.tsv": "contamination",
        "B1_refog_missing_matrix.tsv": "missing_fraction",
    }
    refogs = [r["refog"] for r in method_rows[METHODS[0]]]
    for filename, field in matrix_specs.items():
        matrix = []
        for i, refog in enumerate(refogs):
            if any(method_rows[m][i]["refog"] != refog for m in METHODS):
                raise ValueError(f"RefOG alignment mismatch at {refog}")
            matrix.append({"RefOG": refog, **{m: method_rows[m][i][field] for m in METHODS}})
        write(args.metrics_out / filename, matrix)

    arrays = {m: f(method_rows[m], "F1") for m in METHODS}
    friedman = stats.friedmanchisquare(*(arrays[m] for m in METHODS))
    stats_rows: list[dict[str, object]] = [{
        "test": "Friedman", "comparison": "all_methods", "n_refogs": 70,
        "statistic": float(friedman.statistic), "raw_p": float(friedman.pvalue),
        "Holm_adjusted_p": "NA", "paired_median_difference": "NA", "rank_biserial_effect_size": "NA",
        "discordant_OGProfiler2_only": "NA", "discordant_competitor_only": "NA",
    }]
    wilcoxon_temp = []
    mcnemar_temp = []
    og_exact = np.asarray([r["exact"].lower() == "true" for r in method_rows["OGProfiler2"]])
    for competitor in METHODS[1:]:
        diff = arrays["OGProfiler2"] - arrays[competitor]
        test = stats.wilcoxon(arrays["OGProfiler2"], arrays[competitor], alternative="two-sided", zero_method="wilcox")
        wilcoxon_temp.append((competitor, float(test.statistic), float(test.pvalue), float(np.median(diff)), rank_biserial(diff)))
        comp_exact = np.asarray([r["exact"].lower() == "true" for r in method_rows[competitor]])
        b = int(np.sum(og_exact & ~comp_exact)); c = int(np.sum(~og_exact & comp_exact))
        pvalue = float(stats.binomtest(b, b + c, 0.5, alternative="two-sided").pvalue) if b + c else 1.0
        mcnemar_temp.append((competitor, b, c, pvalue))
    w_adj = holm([x[2] for x in wilcoxon_temp])
    m_adj = holm([x[3] for x in mcnemar_temp])
    for value, adjusted in zip(wilcoxon_temp, w_adj):
        competitor, statistic, raw_p, median_diff, effect = value
        stats_rows.append({
            "test": "Wilcoxon_signed_rank", "comparison": f"OGProfiler2_vs_{competitor}", "n_refogs": 70,
            "statistic": statistic, "raw_p": raw_p, "Holm_adjusted_p": adjusted,
            "paired_median_difference": median_diff, "rank_biserial_effect_size": effect,
            "discordant_OGProfiler2_only": "NA", "discordant_competitor_only": "NA",
        })
    for value, adjusted in zip(mcnemar_temp, m_adj):
        competitor, b, c, raw_p = value
        stats_rows.append({
            "test": "McNemar_exact", "comparison": f"OGProfiler2_vs_{competitor}", "n_refogs": 70,
            "statistic": min(b, c), "raw_p": raw_p, "Holm_adjusted_p": adjusted,
            "paired_median_difference": "NA", "rank_biserial_effect_size": "NA",
            "discordant_OGProfiler2_only": b, "discordant_competitor_only": c,
        })
    write(args.metrics_out / "B1_paired_statistics.tsv", stats_rows)

    rng = np.random.default_rng(args.bootstrap_seed)
    indices = rng.integers(0, 70, size=(args.bootstrap_n, 70))
    bootstrap_rows = []
    for method in METHODS:
        exact = np.asarray([r["exact"].lower() == "true" for r in method_rows[method]], dtype=float)
        for metric, source, function in (
            ("macro_best_group_F1", arrays[method], lambda x: np.mean(x, axis=1)),
            ("median_refog_F1", arrays[method], lambda x: np.median(x, axis=1)),
            ("strict_exact_fraction", exact, lambda x: np.mean(x, axis=1)),
        ):
            estimates = function(source[indices])
            low, high = percentile_ci(estimates)
            point = float(function(source[None, :])[0])
            bootstrap_rows.append({
                "method": method, "metric": metric, "estimate": point,
                "ci95_low": low, "ci95_high": high,
                "n_refogs": 70, "n_resamples": args.bootstrap_n, "seed": args.bootstrap_seed,
            })
    write(args.metrics_out / "B1_bootstrap_CI.tsv", bootstrap_rows)

    args.figure_data.mkdir(parents=True, exist_ok=True)
    write(args.figure_data / "official_precision_recall.tsv", official_summary)
    long_f1 = [{"RefOG": r["refog"], "method": m, "F1": r["F1"]} for m in METHODS for r in method_rows[m]]
    write(args.figure_data / "refog_F1_long.tsv", long_f1)
    split_contam = [{
        "method": r["method"], "median_split": r["median_split"], "mean_split": r["mean_split"],
        "median_contamination": r["median_contamination"], "mean_contamination": r["mean_contamination"],
    } for r in extended_summary]
    write(args.figure_data / "split_contamination.tsv", split_contam)
    fp_long = []
    for record in error_summary:
        for rank in (1, 5, 10):
            fp_long.append({"method": record["method"], "top_n": rank, "cumulative_FP_fraction": record[f"top{rank}_FP_fraction"]})
    write(args.figure_data / "fp_concentration_long.tsv", fp_long)

    ext_by = {r["method"]: r for r in extended_summary}
    err_by = {r["method"]: r for r in error_summary}
    report = [
        "# B1 comparative report", "",
        "Status: **COMPLETE**. Dataset digest: `380c85d9548c607f5df6daf656d74b21597fe215ef78d4f8ff296776bb1fdd07`.", "",
        "## A. Official Open Orthobench", "", md_table(official_summary, ["method", "version", "precision", "recall", "F_score", "official_exact", "status"]), "",
        "## B. Extended RefOG metrics", "", md_table(extended_summary, ["method", "macro_best_group_F1", "median_refog_F1", "strict_exact_n", "median_split", "median_contamination", "median_missing", "VI"]), "",
        "Official exact counts follow `benchmark.py`; strict exact counts require exact set identity in the common 70-RefOG implementation and are therefore reported separately.", "",
        "## C. Pairwise false-positive burden", "", md_table(error_summary, ["method", "top1_FP_fraction", "top5_FP_fraction", "top10_FP_fraction", "largest_predicted_family", "n_families_gt1000", "singleton_fraction"]), "",
        "## D. Paired statistics", "",
        f"The Friedman test across five paired methods used 70 RefOGs (statistic={float(friedman.statistic):.6g}, p={float(friedman.pvalue):.6g}).",
        "Pairwise Wilcoxon signed-rank and exact McNemar tests use Holm correction across the four pre-specified OGProfiler2 comparisons; bootstrap intervals use 10,000 RefOG-level resamples with seed 20260901.", "",
        md_table(stats_rows[1:], ["test", "comparison", "raw_p", "Holm_adjusted_p", "paired_median_difference", "rank_biserial_effect_size", "discordant_OGProfiler2_only", "discordant_competitor_only"]), "",
        "Positive paired median differences and rank-biserial effects favor OGProfiler2. See `B1_bootstrap_CI.tsv` for all percentile intervals.", "",
        "## E. Resource summary", "", md_table(resource_summary, ["method", "wall_seconds", "total_cpu_seconds", "effective_mean_cores", "peak_rss_gib", "disk_gib", "status"]), "",
        "Runtime and resources are descriptive only: each method has one formal run, so no performance inference is made.", "",
        "## F. Error-profile comparison", "",
        f"OrthoFinder3 has the highest official F-score ({official_summary[1]['F_score']}%) and macro per-RefOG F1 ({float(ext_by['OrthoFinder3']['macro_best_group_F1']):.3f}).",
        f"OGProfiler2 combines median split {float(ext_by['OGProfiler2']['median_split']):.3g}, median contamination {float(ext_by['OGProfiler2']['median_contamination']):.3g}, and median missing fraction {float(ext_by['OGProfiler2']['median_missing']):.3g}. Its error burden is heavy-tail concentrated: the largest family contains {err_by['OGProfiler2']['largest_predicted_family']} proteins and the top family accounts for {float(err_by['OGProfiler2']['top1_FP_fraction']):.1%} of pairwise false positives.",
        "- Versus OrthoFinder3, the OGProfiler2 deficit combines more splitting and missing assignments with a much larger heavy-tail family, despite slightly lower mean contamination.",
        "- Versus FastOMA, OGProfiler2 has less splitting and missingness and higher macro F1, but substantially greater contamination concentration; the methods fail in different directions.",
        "- Versus SonicParanoid2, splitting is similar, while OGProfiler2 has more mean contamination and missingness; both show concentrated pairwise false positives.",
        "- Versus Proteinortho6, OGProfiler2 has less splitting and missingness and higher macro F1, whereas Proteinortho6 trades very high precision for low recall and diffuse false positives.",
        "Relative differences are reported as benchmark observations only; no competitor-driven parameter changes are proposed.", "",
        "## G. Benchmark judgment", "",
        "Stage 2A measures accuracy and descriptive resources under frozen inputs, default/recommended competitor modes, and one formal run per method. It does not support runtime significance claims or algorithm retuning.",
    ]
    (args.metrics_out / "B1_COMPARATIVE_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    (args.figure_data / "FIGURE_CONTRACT.md").write_text(
        "Core conclusion: OrthoFinder3 leads overall accuracy, while method-specific gaps separate splitting/missing errors from heavy-tail false-positive concentration.\n"
        "Evidence: A official precision-recall; B 70 paired RefOG F1 values; C median split versus contamination; D top-1/5/10 cumulative FP burden.\n"
        "Archetype: quantitative grid with panel B as the distributional hero panel.\n"
        "Export: 183 mm wide PDF and 600 dpi PNG; editable vector text; all plotted values supplied as TSV.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
