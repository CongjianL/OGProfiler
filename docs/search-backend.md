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
  sensitivity: sensitive
  max_target_seqs: 0
```

`max_target_seqs: 0` means no OGProfiler target-count cap. DIAMOND receives an
explicit zero; MMseqs2 and BLAST+ omit their positive-limit option. A positive
value maps to DIAMOND `--max-target-seqs`, MMseqs2 `--max-seqs`, or BLAST+
`-max_target_seqs`.

DIAMOND sensitivity names map directly to DIAMOND flags. MMseqs2 maps the same
ordered names to `-s` values 2.0, 4.0, 5.7, 6.5, 7.5, and 8.5. BLAST+ has no
direct equivalent and ignores this cross-backend convenience setting; it is
not presented as the performance path.

Every subprocess uses an argument vector with `check=True`; shell parsing and
redirection are absent. One HSP/alignment is requested per query-target pair.
MMseqs2 creates a unique temporary directory beside the raw hit output and
removes it on both success and command failure.

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
