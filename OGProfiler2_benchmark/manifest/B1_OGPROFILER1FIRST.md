# B1 OGProfiler v1 first-version comparator

## Frozen version identity

- Comparator label: `OGProfiler1First`
- Definition: first historical version in this repository
- Commit: `729675fee0dabc02f63ae5c2cda7571748751deb`
- Commit date: 2022-04-02
- Commit subject: `Initial commit add script`
- Historical release tag: none present
- Exported `OGProfiler.py` SHA-256:
  `65ec43d269b410956fadc4f8215a3e0a0772eac8bcc1060ab011b090823cca41`

The exported historical source is byte-identical to `git show
729675fee0dabc02f63ae5c2cda7571748751deb:OGProfiler.py`. It is not patched.

## Input compatibility adapter

The frozen Open Orthobench input remains the scientific input and is audited
against digest
`380c85d9548c607f5df6daf656d74b21597fe215ef78d4f8ff296776bb1fdd07`.
The first historical implementation requires pipe-delimited identifiers while
writing its ancillary pairwise output. A run-local, bijective adapter changes
each tool-facing ID to `OGPV1_NN|ORIGINAL_ID`; `id_map.tsv` records the mapping,
and the converter restores original IDs before partition validation and scoring.
Sequences and proteome membership are unchanged. The frozen input tree itself
is never modified.

## Formal configuration

- Search: DIAMOND 2.0.14, `--more-sensitive`, E-value `1e-5`
- DIAMOND build compatibility: `--ignore-warnings` for `makedb` only
- Distance rule: `lrb`
- Weight: `NBS`
- Community method: `rber`
- Gamma: `1.0`
- Species overlap: `0`
- Search/network threads: `32`/`32`
- Slurm request: 32 CPU, 128 GB, 72 hours
- v1 random-seed control: absent from the historical CLI

All scientific parameters are the historical defaults. Only thread counts are
raised from 8 to the benchmark-wide 32-CPU allocation.

## Preflight history

| Job | State | Diagnosis | Action |
|---|---|---|---|
| 1404586 | FAILED | Raw Open Orthobench IDs lack the pipe delimiter assumed only by v1's final ancillary pairwise writer; the terminal family output had already been produced. | Added a bijective run-local ID adapter and reverse mapping without altering historical source. |
| 1404587 | CANCELLED | A 10-record-per-species smoke sample left two proteomes with only sub-20-aa sequences; DIAMOND refused to create those two databases and v1 waited for failed worker results. | Restricted smoke sampling to the first 10 records of at least 20 aa per species; formal input remains unfiltered. |
| 1404590 | COMPLETED | The 12-species compatibility smoke, terminal-family conversion, reverse ID mapping, and standardized partition validation all passed in 65 seconds. | Formal submission gate passed. |
| 1404591 | TIMEOUT | DIAMOND 2.0.14 rejected the complete human and mouse proteomes as nucleotide-like during `makedb`. Twenty-four directed searches were absent; two matrix workers then exited before reporting to the historical queue, leaving the main process waiting until the 72-hour limit. | Added a run-local DIAMOND launcher that supplies `--ignore-warnings` to `makedb` only. Historical OGProfiler source, sequences, search sensitivity, E-value, and downstream parameters remain unchanged. Retry starts from a fresh immutable run directory. |
| 1407221 | FAILED | The retry prepended the launcher by replacing `PATH` inside `micromamba run`, which also removed the v1 environment's Python path; system Python then reported a missing `igraph` import. | Invoke the v1 environment's Python by absolute path and prepend both the launcher and the complete v1 environment `bin` directory. A dependency-import gate now runs before input preparation. |

The formal run is submitted only after a complete smoke run, standardized
partition conversion, and coverage validation succeed.
