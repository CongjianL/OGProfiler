"""NCBI BLAST+ compatibility backend."""

from __future__ import annotations

import subprocess
from collections.abc import Callable, Sequence
from pathlib import Path

from ogprofiler.exceptions import SearchError
from ogprofiler.search.base import SearchParameters

CommandRunner = Callable[[Sequence[str]], subprocess.CompletedProcess[str]]
BLAST_FIELDS = "6 qseqid sseqid pident length qlen slen evalue bitscore"


def _default_runner(command: Sequence[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(list(command), check=True, capture_output=True, text=True)


class BlastBackend:
    """Execute makeblastdb and blastp as an explicit compatibility path."""

    name = "blastp"

    def __init__(
        self,
        executable: str = "blastp",
        makeblastdb_executable: str = "makeblastdb",
        runner: CommandRunner = _default_runner,
    ) -> None:
        self.executable = executable
        self.makeblastdb_executable = makeblastdb_executable
        self._runner = runner

    def _run(self, command: Sequence[str]) -> subprocess.CompletedProcess[str]:
        try:
            return self._runner(command)
        except FileNotFoundError as error:
            raise SearchError(f"BLAST+ executable was not found: {command[0]}") from error
        except subprocess.CalledProcessError as error:
            detail = (error.stderr or error.stdout or str(error)).strip()
            raise SearchError(f"BLAST+ command failed: {detail}") from error

    def version(self) -> str:
        completed = self._run((self.executable, "-version"))
        output = (completed.stdout or completed.stderr).strip()
        if not output:
            raise SearchError("BLAST+ version command produced no output")
        return output.splitlines()[0]

    def build_database(self, fasta_path: Path, database_path: Path) -> tuple[str, ...]:
        database_path.parent.mkdir(parents=True, exist_ok=True)
        command = (
            self.makeblastdb_executable,
            "-in",
            str(fasta_path),
            "-dbtype",
            "prot",
            "-out",
            str(database_path),
        )
        self._run(command)
        if not self.database_files(database_path):
            raise SearchError(f"BLAST+ database was not created: {database_path}")
        return command

    def database_files(self, database_path: Path) -> tuple[Path, ...]:
        return tuple(
            path
            for path in sorted(database_path.parent.glob(f"{database_path.name}.*"))
            if path.is_file()
        )

    def search(
        self,
        query_path: Path,
        database_path: Path,
        output_path: Path,
        parameters: SearchParameters,
    ) -> tuple[str, ...]:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        command = (
            self.executable,
            "-query",
            str(query_path),
            "-db",
            str(database_path),
            "-out",
            str(output_path),
            "-outfmt",
            BLAST_FIELDS,
            "-evalue",
            str(parameters.evalue),
            "-num_threads",
            str(parameters.threads),
            *(
                ("-max_target_seqs", str(parameters.max_target_seqs))
                if parameters.max_target_seqs
                else ()
            ),
            *(
                ("-max_hsps", str(parameters.max_hsps))
                if parameters.max_hsps
                else ()
            ),
        )
        self._run(command)
        if not output_path.is_file():
            raise SearchError(f"BLAST+ output was not created: {output_path}")
        return command
