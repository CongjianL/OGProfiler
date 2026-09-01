# B1 Stage 2A completion record

`STAGE2A_COMPLETE=YES`

## Frozen scope

- Dataset digest: `380c85d9548c607f5df6daf656d74b21597fe215ef78d4f8ff296776bb1fdd07`
- OGProfiler version: `2.0.0a1`
- OGProfiler algorithm commit: `a7838963de4882d867203e60ce641b5689bde45e`
- Frozen OGProfiler result: `B1_orthobench_OGProfiler2_seed42_rep2`
- Analysis source commit: `4fb6ecff9abc1c55b9b5d75d2a4ff3ff009efbb3`
- Bootstrap: 10,000 RefOG-level resamples, seed `20260901`

No OGProfiler algorithm source, parameter, or seed was changed in Stage 2A.

## Formal execution ledger

| Method | Formal inference job | Postprocess job | Final status |
|---|---:|---:|---|
| OGProfiler2 | `1404237` | included | PASS (frozen baseline) |
| OrthoFinder3 | `1404266` | `1404275` | PASS |
| FastOMA | `1404543` | `1404570` | PASS |
| SonicParanoid2 | `1404443` | included | PASS |
| Proteinortho6 | `1404268` | `1404272` | PASS |

Failed attempts and their diagnoses remain preserved in `B1_STAGE2A_RETRIES.md` and are excluded from formal accuracy and resource tables.

## Frozen outputs

- `05_metrics/orthobench/B1_official_summary.tsv`
- `05_metrics/orthobench/B1_extended_summary.tsv`
- `05_metrics/orthobench/B1_error_structure_summary.tsv`
- `05_metrics/orthobench/B1_refog_F1_matrix.tsv`
- `05_metrics/orthobench/B1_refog_split_matrix.tsv`
- `05_metrics/orthobench/B1_refog_contamination_matrix.tsv`
- `05_metrics/orthobench/B1_refog_missing_matrix.tsv`
- `05_metrics/orthobench/B1_paired_statistics.tsv`
- `05_metrics/orthobench/B1_bootstrap_CI.tsv`
- `05_metrics/orthobench/B1_resource_summary.tsv`
- `05_metrics/orthobench/B1_COMPARATIVE_REPORT.md`
- `07_figures/Fig2_candidate/Fig2_candidate.pdf`
- `07_figures/Fig2_candidate/Fig2_candidate.png`
- `07_figures/Fig2_candidate/data/*.tsv`

The Fig. 2 candidate is preliminary. B1 runtime and resources are descriptive because each method has one formal run.

## Benchmark judgment

OrthoFinder3 has the highest official F-score and macro best-group RefOG F1. OGProfiler2 exceeds FastOMA and Proteinortho6 in macro RefOG F1 but trails OrthoFinder3; its largest weakness is the concentration of pairwise false positives in a small number of very large families, combined with more splitting and missingness than OrthoFinder3. These observations are benchmark findings only and did not trigger algorithm retuning.
