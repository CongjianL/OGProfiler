"""External alignment and tree backend protocols and adapters."""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

from ogprofiler.exceptions import PhylogenyError


@dataclass(frozen=True, slots=True)
class BackendRun:
    command: tuple[str, ...]
    version: str
    stdout: str
    stderr: str


class AlignmentBackend(Protocol):
    @property
    def name(self) -> str: ...

    def version(self) -> str: ...

    def align(self, input_fasta: Path, output_fasta: Path, threads: int) -> BackendRun: ...


class TreeBackend(Protocol):
    @property
    def name(self) -> str: ...

    def version(self) -> str: ...

    def infer(self, alignment_fasta: Path, output_newick: Path) -> BackendRun: ...


def _resolve(executable: str) -> str:
    resolved = shutil.which(executable)
    if resolved is None:
        raise PhylogenyError(f"Executable not found on PATH: {executable}")
    return resolved


def _run(command: list[str], *, allow_nonzero: bool = False) -> subprocess.CompletedProcess[str]:
    try:
        result = subprocess.run(command, text=True, capture_output=True, check=False)
    except OSError as error:
        raise PhylogenyError(f"Failed to execute {command[0]}: {error}") from error
    if result.returncode != 0 and not allow_nonzero:
        detail = result.stderr.strip() or result.stdout.strip() or f"exit {result.returncode}"
        raise PhylogenyError(f"External command failed: {' '.join(command)}: {detail}")
    return result


@dataclass(frozen=True, slots=True)
class MafftBackend:
    executable: str = "mafft"
    name: str = "mafft"

    def version(self) -> str:
        executable = _resolve(self.executable)
        result = _run([executable, "--version"], allow_nonzero=True)
        value = (result.stdout or result.stderr).strip().splitlines()
        return value[0] if value else "unknown"

    def align(self, input_fasta: Path, output_fasta: Path, threads: int) -> BackendRun:
        executable = _resolve(self.executable)
        command = [executable, "--thread", str(threads), "--auto", str(input_fasta)]
        result = _run(command)
        if not result.stdout.lstrip().startswith(">"):
            raise PhylogenyError("MAFFT produced no FASTA alignment")
        output_fasta.write_text(result.stdout, encoding="utf-8", newline="\n")
        return BackendRun(tuple(command), self.version(), result.stdout, result.stderr)


@dataclass(frozen=True, slots=True)
class FastTreeBackend:
    executable: str = "FastTree"
    name: str = "fasttree"

    def version(self) -> str:
        executable = _resolve(self.executable)
        result = _run([executable], allow_nonzero=True)
        lines = (result.stderr or result.stdout).strip().splitlines()
        return lines[0] if lines else "unknown"

    def infer(self, alignment_fasta: Path, output_newick: Path) -> BackendRun:
        executable = _resolve(self.executable)
        command = [executable, "-quiet", str(alignment_fasta)]
        result = _run(command)
        newick = result.stdout.strip()
        if not newick.endswith(";"):
            raise PhylogenyError("FastTree produced invalid Newick output")
        output_newick.write_text(newick + "\n", encoding="utf-8", newline="\n")
        return BackendRun(tuple(command), self.version(), result.stdout, result.stderr)
