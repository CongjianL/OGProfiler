"""MMseqs2 easy-search backend with bounded temporary storage."""

from __future__ import annotations

import shutil
import subprocess
import uuid
from collections.abc import Callable, Sequence
from pathlib import Path

from ogprofiler.exceptions import SearchError
from ogprofiler.search.base import SearchParameters

CommandRunner = Callable[[Sequence[str]], subprocess.CompletedProcess[str]]

MMSEQS_FIELDS = "query,target,pident,alnlen,qlen,tlen,evalue,bits"
_SENSITIVITY = {
    "fast": "2.0",
    "mid-sensitive": "4.0",
    "sensitive": "5.7",
    "more-sensitive": "6.5",
    "very-sensitive": "7.5",
    "ultra-sensitive": "8.5",
}


def _default_runner(command: Sequence[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(list(command), check=True, capture_output=True, text=True)


class MmseqsBackend:
    """Build an MMseqs target database and execute all-vs-all easy-search."""

    name = "mmseqs"

    def __init__(
        self, executable: str = "mmseqs", runner: CommandRunner = _default_runner
    ) -> None:
        self.executable = executable
        self._runner = runner

    def _run(self, command: Sequence[str]) -> subprocess.CompletedProcess[str]:
        try:
            return self._runner(command)
        except FileNotFoundError as error:
            raise SearchError(f"MMseqs2 executable was not found: {self.executable}") from error
        except subprocess.CalledProcessError as error:
            detail = (error.stderr or error.stdout or str(error)).strip()
            raise SearchError(f"MMseqs2 command failed: {detail}") from error

    def version(self) -> str:
        completed = self._run((self.executable, "version"))
        output = (completed.stdout or completed.stderr).strip()
        if not output:
            raise SearchError("MMseqs2 version command produced no output")
        return output

    def build_database(self, fasta_path: Path, database_path: Path) -> tuple[str, ...]:
        database_path.parent.mkdir(parents=True, exist_ok=True)
        command = (self.executable, "createdb", str(fasta_path), str(database_path))
        self._run(command)
        if not self.database_files(database_path):
            raise SearchError(f"MMseqs2 database was not created: {database_path}")
        return command

    def database_files(self, database_path: Path) -> tuple[Path, ...]:
        return tuple(
            path
            for path in sorted(database_path.parent.glob(f"{database_path.name}*"))
            if path.is_file()
        )

    def search(
        self,
        query_path: Path,
        database_path: Path,
        output_path: Path,
        parameters: SearchParameters,
    ) -> tuple[str, ...]:
        try:
            sensitivity = _SENSITIVITY[parameters.sensitivity]
        except KeyError as error:
            raise SearchError(f"Unknown MMseqs2 sensitivity: {parameters.sensitivity}") from error
        output_path.parent.mkdir(parents=True, exist_ok=True)
        temporary = output_path.parent / f".mmseqs-tmp-{uuid.uuid4().hex}"
        command: tuple[str, ...] = (
            self.executable,
            "easy-search",
            str(query_path),
            str(database_path),
            str(output_path),
            str(temporary),
            "--format-output",
            MMSEQS_FIELDS,
            "-e",
            str(parameters.evalue),
            "--threads",
            str(parameters.threads),
            "-s",
            sensitivity,
            *(
                ("--max-seqs", str(parameters.max_target_seqs))
                if parameters.max_target_seqs
                else ()
            ),
        )
        try:
            self._run(command)
        finally:
            shutil.rmtree(temporary, ignore_errors=True)
        if not output_path.is_file():
            raise SearchError(f"MMseqs2 output was not created: {output_path}")
        return command
