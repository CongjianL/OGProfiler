# Similarity search backend

Phase 3 introduces an implementation-independent `SearchBackend` protocol and
the production-path `DiamondBackend`. Search consumes a prepared run workspace:

```bash
ogprofiler prepare --proteomes proteomes/ --out run/
ogprofiler search --run run/ --backend diamond
```

## DIAMOND behavior

The defaults are resolved into `run.yaml`:

```yaml
search:
  backend: diamond
  executable: diamond
  evalue: 0.001
  threads: 8
  sensitivity: sensitive
  max_target_seqs: 0
```

`max_target_seqs: 0` is the OGProfiler all-vs-all default and is passed
explicitly to DIAMOND; it requests no per-query target-count cap. A positive
value opts into truncation and becomes part of the search manifest and resume
identity. Supported sensitivity values are `fast`, `mid-sensitive`,
`sensitive`, `more-sensitive`, `very-sensitive`, and `ultra-sensitive`.

DIAMOND is invoked with argument vectors through
`subprocess.run(..., check=True, capture_output=True, text=True)`. Database and
hit output paths are explicit arguments; shell parsing and redirection are not
used. One HSP is retained per query-target pair (`--max-hsps 1`).

## Directional hit table

`run/search/hits.parquet` contains one directional row per reported hit:

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
normalization, reciprocal-hit logic, and edge symmetrization belong to Phase 4
and are not performed by the search layer. Raw rows are converted to Parquet in
bounded batches so the parser does not retain the complete hit table in memory.

## Provenance and resume

`run/search/search-manifest.json` records:

- backend and reported version;
- top-level, database-build, and search commands;
- resolved parameters and their hash;
- prepared FASTA/metadata checksums;
- database, raw-hit, and Parquet output checksums;
- hit count, schema, and directional semantics.

Search is reused only when the algorithm version, backend and its currently
reported version, parameters, input checksums, and current `hits.parquet`
checksum match the manifest. Missing, changed, or corrupt output triggers a
new database build and search.

MMseqs and BLAST remain protocol slots for later Phase 3 compatibility work;
the Diamond path is the 2.0 MVP search backend.
