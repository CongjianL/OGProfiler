# Similarity search backends

OGProfiler exposes one `SearchBackend` protocol with three implementations:

- **DIAMOND**: production default;
- **MMseqs2**: high-throughput compatibility path using `createdb` and
  `easy-search`;
- **NCBI BLAST+**: conservative compatibility/validation path using
  `makeblastdb` and `blastp`.

All backends consume the prepared global FASTA and produce the same directional
hit table. Select a backend with:

```bash
ogprofiler search --run run/ --backend diamond
ogprofiler search --run run/ --backend mmseqs
ogprofiler search --run run/ --backend blastp
```

An explicit `--backend` selects the conventional executable of the same name.
Use `--set search.executable=/path/to/tool` for a custom installation. For a
custom BLAST path, `makeblastdb` is resolved beside the configured `blastp`.

## Common parameters

```yaml
search:
  backend: diamond
  executable: diamond
  evalue: 0.001
  threads: 8
  sensitivity: more-sensitive
  max_target_seqs: 0
  max_hsps: 0
  mmseqs_sensitivity: sensitive
```

`max_target_seqs` defaults to `0`, meaning the target-count flag is **omitted**
so each backend uses its own default, matching OrthoFinder3/V1 (which do not
pass the flag). A positive value maps to DIAMOND `--max-target-seqs`, MMseqs2
`--max-seqs`, or BLAST+ `-max_target_seqs`.

`max_hsps` defaults to `0` (do not pass the flag, matching OrthoFinder3/V1). A
positive value maps to DIAMOND/BLAST+ `--max-hsps`; MMseqs2 has no direct
equivalent.

`sensitivity` is the DIAMOND sensitivity name and maps directly to DIAMOND
flags (`more-sensitive` by default). MMseqs2 uses `mmseqs_sensitivity`
(`sensitive` → `-s 5.7` by default, matching OrthoFinder3's untuned MMseqs2
call) and ignores the shared `sensitivity`. BLAST+ has no direct equivalent and
ignores these convenience settings.

Every subprocess uses an argument vector with `check=True`; shell parsing and
redirection are absent. `max_hsps` controls whether an HSP-limit flag is
passed; by default (`0`) it is omitted, matching OrthoFinder3/V1. MMseqs2
creates a unique temporary directory beside the raw hit output and removes it
on both success and command failure.

## Standard directional hit table

Each backend requests these eight raw fields in this order:

```text
query ID, target ID, identity, aligned length,
query length, target length, e-value, bit score
```

The shared bounded-memory parser writes `run/search/hits.parquet`:

```text
query_id          int64
target_id         int64
query_species     int32
target_species    int32
bitscore          float64
identity          float32
query_coverage    float32
target_coverage   float32
evalue            float64
```

Coverage is `100 * aligned_length / sequence_length`. Coverage filtering,
normalization, reciprocal-hit logic, and edge symmetrization remain downstream
edge-stage responsibilities.

## Provenance and resume

`run/search/search-manifest.json` records backend/version, database and search
commands, resolved parameters, input checksums, every database artifact, raw
hits, normalized Parquet hits, row count, and schema. Verified resume requires
matching algorithm version, backend version, parameters, prepared inputs, and
current normalized hit checksum. Backend or parameter changes invalidate reuse.
