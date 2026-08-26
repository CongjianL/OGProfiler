from __future__ import annotations

import subprocess
from collections.abc import Sequence
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from ogprofiler.exceptions import SearchError
from ogprofiler.search.base import SearchBackend, SearchParameters
from ogprofiler.search.blast import BLAST_FIELDS, BlastBackend
from ogprofiler.search.diamond import DIAMOND_FIELDS, DiamondBackend
from ogprofiler.search.hits import HIT_COLUMNS, parse_diamond_hits
from ogprofiler.search.mmseqs import MMSEQS_FIELDS, MmseqsBackend


def test_diamond_backend_builds_explicit_argument_vectors(tmp_path: Path) -> None:
    commands: list[tuple[str, ...]] = []

    def runner(command: Sequence[str]) -> subprocess.CompletedProcess[str]:
        captured = tuple(command)
        commands.append(captured)
        if "version" in captured:
            return subprocess.CompletedProcess(captured, 0, "diamond version 2.test\n", "")
        if "makedb" in captured:
            database = Path(captured[captured.index("--db") + 1]).with_suffix(".dmnd")
            database.write_text("db", encoding="utf-8")
        if "blastp" in captured:
            output = Path(captured[captured.index("--out") + 1])
            output.write_text("", encoding="utf-8")
        return subprocess.CompletedProcess(captured, 0, "", "")

    backend = DiamondBackend("diamond-test", runner)
    fasta = tmp_path / "proteins.faa"
    fasta.write_text(">OGP2P000000000000\nAAAA\n", encoding="utf-8")
    database = tmp_path / "database" / "proteins"
    output = tmp_path / "hits.tsv"
    parameters = SearchParameters(4, 1e-5, "very-sensitive", 0)

    assert backend.version() == "diamond version 2.test"
    backend.build_database(fasta, database)
    search_command = backend.search(fasta, database, output, parameters)

    assert search_command[:2] == ("diamond-test", "blastp")
    assert search_command[search_command.index("--outfmt") + 2 :][: len(DIAMOND_FIELDS)] == (
        DIAMOND_FIELDS
    )
    assert search_command[search_command.index("--max-target-seqs") + 1] == "0"
    assert "--very-sensitive" in search_command
    assert all("<" not in item and ">" not in item for command in commands for item in command)


def test_diamond_failure_is_search_error() -> None:
    def runner(command: Sequence[str]) -> subprocess.CompletedProcess[str]:
        raise subprocess.CalledProcessError(1, command, stderr="synthetic failure")

    with pytest.raises(SearchError, match="synthetic failure"):
        DiamondBackend(runner=runner).version()


def test_parse_diamond_hits_writes_directional_standard_schema(tmp_path: Path) -> None:
    proteins = pa.table(
        {
            "protein_id": pa.array([0, 1], type=pa.int64()),
            "species_id": pa.array([3, 7], type=pa.int32()),
        }
    )
    protein_path = tmp_path / "proteins.parquet"
    pq.write_table(proteins, protein_path)
    raw_path = tmp_path / "hits.tsv"
    raw_path.write_text(
        "OGP2P000000000001\tOGP2P000000000000\t75\t8\t8\t10\t1e-20\t42\n"
        "OGP2P000000000000\tOGP2P000000000001\t50\t4\t10\t8\t1e-5\t20\n",
        encoding="utf-8",
    )
    output = tmp_path / "hits.parquet"

    assert parse_diamond_hits(raw_path, protein_path, output, batch_size=1) == 2
    table = pq.read_table(output)
    assert tuple(table.column_names) == HIT_COLUMNS
    rows = table.to_pylist()
    assert rows[0]["query_id"] == 1
    assert rows[0]["query_species"] == 7
    assert rows[0]["target_species"] == 3
    assert rows[0]["query_coverage"] == pytest.approx(100.0)
    assert rows[0]["target_coverage"] == pytest.approx(80.0)
    assert rows[1]["query_id"] == 0
    assert rows[1]["target_id"] == 1


def test_diamond_satisfies_search_backend_protocol() -> None:
    assert isinstance(DiamondBackend(), SearchBackend)


def test_mmseqs_backend_uses_easy_search_and_cleans_temporary_directory(
    tmp_path: Path,
) -> None:
    commands: list[tuple[str, ...]] = []

    def runner(command: Sequence[str]) -> subprocess.CompletedProcess[str]:
        captured = tuple(command)
        commands.append(captured)
        if captured[1] == "version":
            return subprocess.CompletedProcess(captured, 0, "15.6f452\n", "")
        if captured[1] == "createdb":
            Path(captured[3]).write_text("db", encoding="utf-8")
            Path(captured[3] + ".dbtype").write_text("0", encoding="utf-8")
        if captured[1] == "easy-search":
            Path(captured[4]).write_text("", encoding="utf-8")
            Path(captured[5]).mkdir(parents=True)
        return subprocess.CompletedProcess(captured, 0, "", "")

    backend = MmseqsBackend("mmseqs-test", runner)
    fasta = tmp_path / "proteins.faa"
    fasta.write_text(">OGP2P000000000000\nAAAA\n", encoding="utf-8")
    database = tmp_path / "database/proteins"
    output = tmp_path / "hits.tsv"
    parameters = SearchParameters(4, 1e-5, "very-sensitive", 25)
    assert backend.version() == "15.6f452"
    backend.build_database(fasta, database)
    command = backend.search(fasta, database, output, parameters)
    assert command[1] == "easy-search"
    assert command[command.index("--format-output") + 1] == MMSEQS_FIELDS
    assert command[command.index("--max-seqs") + 1] == "25"
    assert not Path(command[5]).exists()
    assert isinstance(backend, SearchBackend)


def test_blast_backend_builds_database_and_omits_zero_target_limit(tmp_path: Path) -> None:
    commands: list[tuple[str, ...]] = []

    def runner(command: Sequence[str]) -> subprocess.CompletedProcess[str]:
        captured = tuple(command)
        commands.append(captured)
        if "-version" in captured:
            return subprocess.CompletedProcess(captured, 0, "blastp: 2.16.0+\n", "")
        if captured[0] == "makeblastdb-test":
            prefix = Path(captured[captured.index("-out") + 1])
            prefix.with_suffix(".pin").write_text("db", encoding="utf-8")
        if captured[0] == "blastp-test" and "-query" in captured:
            Path(captured[captured.index("-out") + 1]).write_text("", encoding="utf-8")
        return subprocess.CompletedProcess(captured, 0, "", "")

    backend = BlastBackend("blastp-test", "makeblastdb-test", runner)
    fasta = tmp_path / "proteins.faa"
    fasta.write_text(">OGP2P000000000000\nAAAA\n", encoding="utf-8")
    database = tmp_path / "database/proteins"
    output = tmp_path / "hits.tsv"
    parameters = SearchParameters(2, 1e-3, "sensitive", 0)
    assert backend.version() == "blastp: 2.16.0+"
    backend.build_database(fasta, database)
    command = backend.search(fasta, database, output, parameters)
    assert command[command.index("-outfmt") + 1] == BLAST_FIELDS
    assert "-max_target_seqs" not in command
    assert isinstance(backend, SearchBackend)
