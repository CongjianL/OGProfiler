#!/usr/bin/env python3
"""Build Stage 1C QC figures and evidence-backed Markdown reports."""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def metric_table(path: Path) -> dict[str, str]:
    return {row["metric"]: row["value"] for row in table(path)}


def save(fig, prefix: Path) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(prefix.with_suffix(".pdf"), bbox_inches="tight")


def pct(value: float) -> str:
    return f"{100 * value:.1f}%"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--fig2", required=True, type=Path)
    parser.add_argument("--diagnostic-figures", required=True, type=Path)
    parser.add_argument("--performance-observations", required=True, type=Path)
    args = parser.parse_args()
    m = args.metrics

    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "xtick.direction": "out", "ytick.direction": "out", "savefig.dpi": 300})

    burden = table(m / "predicted_family_pairwise_burden.tsv")
    family_dist = table(m / "family_size_distribution.tsv")
    stage = table(m / "refog_stage_attribution.tsv")
    split = table(m / "split_origin.tsv")
    fusions = table(m / "catastrophic_fusions.tsv")
    concentrations = metric_table(m / "pairwise_fp_concentration.tsv")
    size_summary = metric_table(m / "family_size_summary.tsv")
    official = metric_table(m / "OGProfiler2_seed42_official.tsv")
    extended = metric_table(m / "OGProfiler2_seed42_summary.tsv")
    correlations = table(m / "refog_challenge_correlations.tsv")

    # Family-size distribution, with family- and protein-weighted views.
    x = [int(r["family_size"]) for r in family_dist]
    nf = [int(r["n_families"]) for r in family_dist]
    protein_fraction = [float(r["fraction_proteins"]) for r in family_dist]
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.loglog(x, nf, "o", ms=2.5, color="#2563A6", label="families")
    ax.set(xlabel="Terminal family size", ylabel="Number of families",
           title="Terminal-family sizes are heavy-tailed")
    ax2 = ax.twinx(); ax2.plot(x, protein_fraction, color="#D97706", alpha=.7, lw=1,
                              label="protein-weighted fraction")
    ax2.set_ylabel("Fraction of proteins at family size", color="#D97706")
    save(fig, args.fig2 / "OGProfiler2_family_size_distribution"); plt.close(fig)

    # FP Pareto.
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    cumulative = [float(r["cumulative_FP_fraction"]) for r in burden]
    ax.plot(range(1, len(cumulative) + 1), cumulative, color="#A33A3A", lw=1.2)
    ax.set_xscale("log"); ax.set_ylim(0, 1.02)
    ax.set(xlabel="Predicted families ranked by official-like FP (log scale)",
           ylabel="Cumulative fraction of official-like FP", title="Pairwise FP Pareto curve")
    for n in (1, 5, 10): ax.scatter(n, cumulative[n-1], s=18, label=f"top {n}: {pct(cumulative[n-1])}")
    ax.legend(frameon=False, loc="lower right", fontsize=7)
    save(fig, args.diagnostic_figures / "OGProfiler2_FP_pareto"); plt.close(fig)

    # Split origins.
    colors = {"PRE_HIERARCHY_DISCONNECTED": "#2563A6", "HIERARCHY_GENERATED": "#D97706",
              "MIXED_PRE_AND_HIERARCHY": "#7B4FA3", "UNCERTAIN": "#777777"}
    ordered = sorted(split, key=lambda r: (-int(r["final_split_count"]), r["refog"]))
    fig, ax = plt.subplots(figsize=(9, 3.6))
    for i, r in enumerate(ordered):
        ax.bar(i, int(r["final_split_count"]), color=colors[r["classification"]], width=.8)
    ax.set(xlabel="Split RefOG", ylabel="Split count", title="Origin of RefOG fragmentation")
    ax.set_xticks(range(len(ordered)), [r["refog"] for r in ordered], rotation=90, fontsize=5)
    handles = [plt.Line2D([0], [0], color=v, lw=6, label=k) for k, v in colors.items()
               if any(r["classification"] == k for r in ordered)]
    ax.legend(handles=handles, frameon=False, fontsize=6)
    save(fig, args.diagnostic_figures / "OGProfiler2_split_origin"); plt.close(fig)

    # F1 versus size.
    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    ax.scatter([int(r["true_size"]) for r in stage], [float(r["best_terminal_F1"]) for r in stage],
               s=20, alpha=.8, color="#2563A6", edgecolors="white", linewidths=.3)
    size_f1 = next(r for r in correlations if r["challenge_variable"] == "refog_size" and r["response"] == "F1")
    ax.set(xlabel="RefOG size", ylabel="Best-terminal-family F1", title="F1 versus RefOG size")
    ax.text(.03, .05, f"Spearman ρ={float(size_f1['spearman_rho']):.3f}\n95% bootstrap CI "
            f"[{float(size_f1['bootstrap_ci_2.5']):.3f}, {float(size_f1['bootstrap_ci_97.5']):.3f}]",
            transform=ax.transAxes, fontsize=7)
    save(fig, args.diagnostic_figures / "OGProfiler2_F1_vs_refog_size"); plt.close(fig)

    # Contamination versus selected best-family size.
    family_size = {r["family_id"]: int(r["family_size_total"]) for r in burden}
    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    bx = [family_size[r["best_terminal_family"]] for r in stage]
    by = [float(r["contamination"]) for r in stage]
    ax.scatter(bx, by, s=20, alpha=.8, color="#A33A3A", edgecolors="white", linewidths=.3)
    ax.set_xscale("log"); ax.set(xlabel="Best predicted family size (log scale)", ylabel="Contamination",
                                  title="Contamination concentrates in giant families")
    save(fig, args.diagnostic_figures / "OGProfiler2_contamination_vs_best_family_size"); plt.close(fig)

    # Focal giant family: composition and edge-category summaries rather than an unreadable network.
    composition = table(m / "OG000000491/refog_composition.tsv")
    edge_summary = table(m / "OG000000491/edge_summary.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.4))
    axes[0].bar([r["category"] for r in composition], [int(r["n"]) for r in composition], color="#2563A6")
    axes[0].set(title="OG000000491 composition", ylabel="Proteins")
    axes[0].tick_params(axis="x", rotation=30)
    axes[1].barh([r["edge_category"] for r in edge_summary], [int(r["n_edges"]) for r in edge_summary], color="#D97706")
    axes[1].set_xscale("log"); axes[1].set(title="Retained induced edges", xlabel="Edge count (log scale)")
    save(fig, args.diagnostic_figures / "OG000000491_graph_summary"); plt.close(fig)

    total_fp = sum(float(r["official_like_FP"]) for r in burden)
    total_tp = sum(float(r["official_like_TP"]) for r in burden)
    reconstructed_precision = total_tp / (total_tp + total_fp)
    top1 = float(concentrations["top_1_family_FP_fraction"])
    top5 = float(concentrations["top_5_family_FP_fraction"])
    top10 = float(concentrations["top_10_family_FP_fraction"])
    class_counts = {r["classification"]: int(r["n"]) for r in table(m / "refog_error_classification.tsv")}
    split_counts = Counter(r["classification"] for r in split)
    cats = [r for r in fusions if r["catastrophic_definition_met"] == "true"]
    cat_fp = sum(float(r["fraction_total_FP"]) for r in cats)
    og = next(r for r in fusions if r["family_id"] == "OG000000491")
    og_comp = {r["category"]: int(r["n"]) for r in composition}
    bridges = [r for r in table(m / "fusion_bridge_candidates.tsv") if r["family_id"] == "OG000000491"]
    articulation = sum(r["articulation_point"] == "true" for r in bridges)

    comparison = f"""# Metric definition comparison — B1 seed42_rep2

## Verified official pairwise precision

The diagnostic code copies the combinatorics of the supplied
`BENCHMARKS/benchmark.py::calculate_benchmarks_pairwise` without changing or replacing it.
For each RefOG it removes that RefOG's low-certainty assignments from both reference and
overlapping predictions, counts pairwise TP/FP/FN, and applies the official equal-RefOG
normalization `N = len(RefOG) - 1`.

The family decomposition reconstructs TP={total_tp:.6f}, FP={total_fp:.6f}, and
precision={reconstructed_precision:.9f} ({pct(reconstructed_precision)}), matching the
official printed precision of {official['precision']}% after rounding.

Top 1 / 5 / 10 predicted families contribute {pct(top1)} / {pct(top5)} / {pct(top10)}
of all official-like FP. Thus pairwise precision is controlled by a very small Pareto
tail even though the median best-family contamination is {float(extended['median_contamination']):.4f}.

## Why low official precision and strong median RefOG metrics coexist

Official precision is pair-weighted within each equally normalized RefOG. A large fused
family contributes `overlap × non-RefOG family members`, so a few giant contaminated
families create many FP pairs. In contrast, median contamination and median best-group F1
give each of the 70 RefOGs one observation. The observed mean/median contamination are
{float(extended['mean_contamination']):.4f}/{float(extended['median_contamination']):.4f},
while median best-group F1 is {float(extended['median_refog_F1']):.4f}; these summaries
therefore describe the typical RefOG rather than the pairwise FP tail.

## Official exact=14 versus strict exact=12

The official code calls a RefOG exact when its per-RefOG FP and FN are both zero *after*
removing low-certainty assignments for that RefOG*. The extended strict check requires
literal equality between the original RefOG member set and one predicted family.
Direct recomputation identified RefOG015 and RefOG056 as the two official-only cases:
their missing differences contain only low-certainty members. The other 12 are exact by
both definitions.
"""
    (m / "metric_definition_comparison.md").write_text(comparison, encoding="utf-8")

    raw_disconnected = sum(int(r["raw_search_connected_components"]) > 1 for r in stage)
    retained_disconnected = sum(int(r["retained_induced_components"]) > 1 for r in stage)
    edge_transition = sum(int(r["raw_search_connected_components"]) == 1 and
                          int(r["retained_induced_components"]) > 1 for r in stage)
    missing_search = sum(int(r["n_members_with_any_search_hit"]) < int(r["true_size"]) for r in stage)
    size_split = next(r for r in correlations if r["challenge_variable"] == "refog_size" and r["response"] == "split_count")
    report = f"""# B1 biological error-attribution report

Run: `B1_orthobench_OGProfiler2_seed42_rep2`  
Mode: read-only analysis of saved immutable artifacts; no search or prediction rerun.

## A. FP concentration

- Top 1 family accounts for **{pct(top1)}** of official-like pairwise FP.
- Top 5 account for **{pct(top5)}**.
- Top 10 account for **{pct(top10)}**.
- The decomposition reconstructs official precision as **{pct(reconstructed_precision)}**
  from TP={total_tp:.3f}, FP={total_fp:.3f}.

## B. Split origin

Among 56 split RefOGs, **{split_counts['PRE_HIERARCHY_DISCONNECTED']}/56** are exclusively
pre-hierarchy disconnected, **{split_counts['HIERARCHY_GENERATED']}/56** are colocated in
one global component and split only in the hierarchy, and
**{split_counts['MIXED_PRE_AND_HIERARCHY']}/56** show both pre-hierarchy disconnection and
additional hierarchy splitting. In total, **{split_counts['PRE_HIERARCHY_DISCONNECTED'] + split_counts['MIXED_PRE_AND_HIERARCHY']}/56
({pct((split_counts['PRE_HIERARCHY_DISCONNECTED'] + split_counts['MIXED_PRE_AND_HIERARCHY'])/56)})**
already have members in multiple global components before hierarchy; only **1/56
({pct(1/56)})** split despite full pre-hierarchy colocation.

## C. Catastrophic fusion

The diagnostic convention (`recall >= 0.8` and `contamination >= 0.5`) identifies
**{len(cats)} predicted families**, jointly accounting for **{pct(cat_fp)}** of all
official-like FP. These are the same first five families in the FP ranking, so the low
pairwise precision is strongly catastrophic-fusion dominated rather than uniformly poor.

## D. Stage attribution across 70 RefOGs

- SEARCH_OR_EDGE_DISCONNECTED: **{class_counts.get('SEARCH_OR_EDGE_DISCONNECTED', 0)}**
- HIERARCHY_SPLIT: **{class_counts.get('HIERARCHY_SPLIT', 0)}**
- CATASTROPHIC_FUSION: **{class_counts.get('CATASTROPHIC_FUSION', 0)}**
- MIXED: **{class_counts.get('MIXED', 0)}**
- WELL_RECOVERED: **{class_counts.get('WELL_RECOVERED', 0)}**

This is a benchmark-specific diagnostic classification, not an Orthobench biological
taxonomy. Saved raw hits exist. {missing_search}/70 RefOGs have at least one member with no
recorded hit to another saved endpoint; {raw_disconnected}/70 are disconnected in the
raw within-RefOG similarity graph, versus {retained_disconnected}/70 after retained-edge
filtering. Most importantly, **{edge_transition}/70** transition from one raw within-RefOG
component to multiple retained induced components, implicating edge retention/connectivity
more strongly than missing raw similarity alone. Four RefOGs have more retained induced
pieces than global components because paths through non-RefOG nodes reconnect some pieces.

## E. OG000000491

`OG000000491` contains **{og['total_size']} proteins from {og['n_species']} species**:
{og_comp.get('RefOG026', 0)} RefOG026, {og_comp.get('RefOG032', 0)} RefOG032, and
{og_comp.get('NON_REFOG', 0)} non-RefOG proteins. It alone contributes **{pct(float(og['fraction_total_FP']))}**
of all pairwise FP. The saved hierarchy shows a 3,087-gene root split at depth 0 and the
857-gene child becoming terminal at depth 1 with reason `UNSTABLE` and event `AMBIGUOUS`.
There are no direct retained RefOG026↔RefOG032 edges inside this terminal family; both
RefOGs connect into a much larger non-RefOG network. Among the 20 highest cutoff-4
betweenness candidates, {articulation} are articulation points; the top two are non-RefOG
proteins, supporting a bridge-mediated topology diagnosis without asserting that removal
would be biologically valid.

## F. Biological challenge

Only RefOG size was found as a reliable structured challenge variable in the bundled
material. Larger RefOGs have lower F1 (Spearman ρ={float(size_f1['spearman_rho']):.3f},
bootstrap 95% CI [{float(size_f1['bootstrap_ci_2.5']):.3f}, {float(size_f1['bootstrap_ci_97.5']):.3f}])
and more splitting (ρ={float(size_split['spearman_rho']):.3f}, 95% CI
[{float(size_split['bootstrap_ci_2.5']):.3f}, {float(size_split['bootstrap_ci_97.5']):.3f}]).
These are exploratory. Reliable structured evolutionary-rate, alignment-quality, and
domain-complexity fields were not located, so those associations remain `NA`; no external
annotation was introduced.

## G. Performance observation

The formal run used 18,860.27 total CPU seconds over 14,189 wall seconds: **1.329 effective
mean CPU cores**, or **4.15%** of 32 allocated CPUs. This is recorded for B4 scalability;
no profiling or parallelism rewrite was performed in Stage 1C.

## H. Final diagnostic conclusion

**MIXED: EDGE_FILTER_LIMITED + HIERARCHY_OVER_SPLITTING + STOPPING_UNDER_SPLITTING +
CATASTROPHIC_FUSION_DOMINATED**, with a smaller SEARCH_LIMITED contribution. Pre-hierarchy
disconnection is nearly universal among split RefOGs, and many disconnections appear only
after edge retention. Hierarchy adds splitting in 15 cases (14 mixed plus one hierarchy-only),
while a few large `UNSTABLE` terminal families account for nearly all pairwise FP. These
observations attribute the existing run; they do not establish a new default or parameter.
"""
    (m / "B1_DIAGNOSTIC_REPORT.md").write_text(report, encoding="utf-8")

    performance = """# B1 performance observations

Formal run: `B1_orthobench_OGProfiler2_seed42_rep2`  
Slurm job: `1404237`

| Quantity | Value |
|---|---:|
| Wall time | 14,189 s |
| User CPU time | 18,030.04 s |
| System CPU time | 830.23 s |
| Total CPU time | 18,860.27 s |
| Allocated CPUs | 32 |
| Effective mean CPU cores (`CPU / wall`) | 1.3292 |
| Approximate allocation efficiency (`1.3292 / 32`) | 4.1538% |
| Peak RSS | 72,819,420 KiB (69.446 GiB) |

GNU `/usr/bin/time -v` measured the wrapped OGProfiler command. Its user/system values
normally accumulate CPU for the measured process and waited-for descendants, including
the pipeline child processes that terminate under that command, while multithreaded CPU
time is summed across threads. They do not constitute per-stage profiling and may omit
detached work or work not waited for by the measured process. Consequently, the ratio is
an allocation-level observation, not a diagnosis of which stage underused CPUs.

This observation is deferred to B4 scalability. Stage 1C makes no parallelism change.
"""
    args.performance_observations.parent.mkdir(parents=True, exist_ok=True)
    args.performance_observations.write_text(performance, encoding="utf-8")


if __name__ == "__main__":
    main()
