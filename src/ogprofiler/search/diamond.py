"""DIAMOND all-vs-all search backend."""

from __future__ import annotations

import subprocess
from collections.abc import Callable, Sequence
from pathlib import Path

from ogprofiler.exceptions import SearchError
from ogprofiler.search.base import SearchParameters

CommandRunner = Callable[[Sequence[str]], subprocess.CompletedProcess[str]]

_SENSITIVITY_FLAGS = {
    "fast": "--fast",
    "mid-sensitive": "--mid-sensitive",
    "sensitive": "--sensitive",
    "more-sensitive": "--more-sensitive",
    "very-sensitive": "--very-sensitive",
    "ultra-sensitive": "--ultra-sensitive",
}

DIAMOND_FIELDS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "slen",
    "evalue",
    "bitscore",
)


def _default_runner(command: Sequence[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(command),
        check=True,
        capture_output=True,
        text=True,
    )


class DiamondBackend:
    """Execute DIAMOND with argument vectors and explicit output paths."""

    name = "diamond"

    def __init__(
        self,
        executable: str = "diamond",
        runner: CommandRunner = _default_runner,
    ) -> None:
        self.executable = executable
        self._runner = runner

    def _run(self, command: Sequence[str]) -> subprocess.CompletedProcess[str]:
        try:
            return self._runner(command)
        except FileNotFoundError as error:
            raise SearchError(f"DIAMOND executable was not found: {self.executable}") from error
        except subprocess.CalledProcessError as error:
            detail = (error.stderr or error.stdout or str(error)).strip()
            raise SearchError(f"DIAMOND command failed: {detail}") from error

    def version(self) -> str:
        completed = self._run((self.executable, "version"))
        output = (completed.stdout or completed.stderr).strip()
        if not output:
            raise SearchError("DIAMOND version command produced no output")
        return output

    def build_database(self, fasta_path: Path, database_path: Path) -> tuple[str, ...]:
        database_path.parent.mkdir(parents=True, exist_ok=True)
        command = (
            self.executable,
            "makedb",
            "--in",
            str(fasta_path),
            "--db",
            str(database_path),
        )
        self._run(command)
        expected = database_path.with_suffix(".dmnd")
        if not expected.is_file():
            raise SearchError(f"DIAMOND database was not created: {expected}")
        return command

    def database_files(self, database_path: Path) -> tuple[Path, ...]:
        return (database_path.with_suffix(".dmnd"),)

    def search(
        self,
        query_path: Path,
        database_path: Path,
        output_path: Path,
        parameters: SearchParameters,
    ) -> tuple[str, ...]:
        try:
            sensitivity_flag = _SENSITIVITY_FLAGS[parameters.sensitivity]
        except KeyError as error:
            allowed = ", ".join(sorted(_SENSITIVITY_FLAGS))
            raise SearchError(f"Unknown DIAMOND sensitivity; expected: {allowed}") from error
        output_path.parent.mkdir(parents=True, exist_ok=True)
        command = (
            self.executable,
            "blastp",
            "--query",
            str(query_path),
            "--db",
            str(database_path),
            "--out",
            str(output_path),
            "--outfmt",
            "6",
            *DIAMOND_FIELDS,
            "--evalue",
            str(parameters.evalue),
            "--threads",
            str(parameters.threads),
            sensitivity_flag,
            *(
                ("--max-target-seqs", str(parameters.max_target_seqs))
                if parameters.max_target_seqs
                else ()
            ),
            *(
                ("--max-hsps", str(parameters.max_hsps))
                if parameters.max_hsps
                else ()
            ),
        )
        self._run(command)
        if not output_path.is_file():
            raise SearchError(f"DIAMOND output was not created: {output_path}")
        return command
