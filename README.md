# OGProfiler 2

OGProfiler 2 builds deterministic, hierarchical protein-family assignments from
sequence-similarity networks. It combines a directional homology search,
legacy-compatible normalized bit scores, sparse connected components,
hierarchical Leiden subdivision, stable terminal-family IDs, and optional
phylogenetic evidence. The implementation lives in `src/ogprofiler`; the frozen
V1 reference remains under `legacy/` for regression comparison.

## Scientific model

The default pipeline is:

```text
proteomes → all-vs-all search → normalized retained edges → connected components
          → hierarchical Leiden → network-event annotation → stable result tables
```

A terminal family is a hierarchy node at which subdivision stops because of an
explicit rule such as one species, minimum size, depth, absent edges, instability,
or exhausted resolution search. The network hierarchy describes nested
community structure in a similarity graph. It is not a gene tree.

## Installation

Python 3.10 or newer is required. DIAMOND is the default search backend.

### Conda / micromamba

```bash
micromamba create -f environment.yml
micromamba activate ogprofiler
python -m pip install -e .
```

### Python development environment

Install DIAMOND separately, then run:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[dev]'
```

Confirm the installation with `ogprofiler --version`.

## Quick start

Run or resume the complete default workflow:

```bash
ogprofiler run --proteomes proteomes/ --out run/ \
  --set search.threads=8 \
  --set runtime.workers=4
```

Preview the exact stage commands without executing them:

```bash
ogprofiler run --proteomes proteomes/ --out run/ --dry-run
```

Resume an existing workspace from a stage:

```bash
ogprofiler run --out run/ --from-stage hierarchy
```

Monitor and inspect it:

```bash
ogprofiler status --run run/
ogprofiler status --run run/ --json
ogprofiler inspect --run run/
ogprofiler inspect --run run/ --component 0 --json
```

Every pipeline stage also remains directly callable:

```bash
ogprofiler prepare --proteomes proteomes/ --out run/
ogprofiler search --run run/ --backend diamond
ogprofiler edges --run run/ --method lrb
ogprofiler components --run run/
ogprofiler hierarchy-all --run run/ --set runtime.workers=4
ogprofiler annotate-network --run run/
ogprofiler export --run run/
```

DIAMOND is the production default. MMseqs2 and NCBI BLAST+ are interchangeable
compatibility backends and emit the same directional Parquet schema:

```bash
ogprofiler search --run run/ --backend mmseqs
ogprofiler search --run run/ --backend blastp
```

## Configuration and reproducibility

Configuration is YAML-based. `prepare` writes the fully resolved `run.yaml`;
`--set section.key=value` overrides a setting. The random seed defaults to 42.
Each stage records input checksums, parameters, algorithm versions, commands,
and output checksums. `ogprofiler run` additionally writes `provenance.json`
with package, Python, platform, configuration, dataset-manifest, and external
executable versions.

For formal runs, preserve the complete workspace, use a clean tagged source
revision, pin `environment.yml`, and retain every manifest. See
[`docs/reproducibility.md`](docs/reproducibility.md).

## Outputs

The stable exchange tables under `run/results/` are:

| File | Meaning |
|---|---|
| `families.tsv` | One row per terminal family with stable dataset-scoped OG ID |
| `members.tsv` | Protein-to-family membership with species and original IDs |
| `hierarchy.tsv` | Parent/child hierarchy nodes, resolution, quality, and stop reason |
| `events.tsv` | Network-derived event labels and confidence |

`export-manifest.json` checksums these outputs. Per-family FASTA, component
GraphML, pairwise ortholog candidates, and phylogenetic refinement are explicit
opt-ins because they may be large or require external tools:

```bash
ogprofiler export --run run/ --family-fasta OG000000123
ogprofiler export graph --run run/ --component 0
ogprofiler orthologs --run run/ --emit-pairwise-orthologs
ogprofiler annotate --run run/ --phylogenetic-refinement --family OG000000123
```

## Scientific caveats

- Network hierarchy is not gene genealogy and does not establish duplication,
  loss, or speciation history by itself.
- `network_event` is a species-overlap heuristic. It must remain separate from
  `phylo_event`, which is derived from alignment, tree inference, rooting, and
  reconciliation.
- Similarity thresholds, proteome completeness, fusion/domain architecture,
  divergence, lineage expansion, and search sensitivity affect family recovery.
- Pairwise ortholog output can grow quadratically and is disabled by default.
- A successful process or Slurm job establishes execution success, not
  biological validity. Interpret results against controls and benchmark truth.

## Development and release verification

```bash
ruff check .
mypy src/ogprofiler
pytest -q
python benchmarks/generate_phase0_datasets.py --check
```

The phased design, benchmark evidence, and release checklist are documented in
`docs/`. Current package version: **2.0.0a1**.
