"""Search backend protocol and implementation-independent parameters."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, runtime_checkable


@dataclass(frozen=True, slots=True)
class SearchParameters:
    threads: int
    evalue: float
    sensitivity: str
    max_target_seqs: int

    def to_dict(self) -> dict[str, int | float | str]:
        return {
            "threads": self.threads,
            "evalue": self.evalue,
            "sensitivity": self.sensitivity,
            "max_target_seqs": self.max_target_seqs,
        }


@dataclass(frozen=True, slots=True)
class SearchStageResult:
    hits_path: Path
    manifest_path: Path
    reused: bool


@runtime_checkable
class SearchBackend(Protocol):
    """Minimal replaceable interface for homolog-search tools."""

    @property
    def name(self) -> str: ...

    def version(self) -> str: ...

    def build_database(self, fasta_path: Path, database_path: Path) -> tuple[str, ...]: ...

    def search(
        self,
        query_path: Path,
        database_path: Path,
        output_path: Path,
        parameters: SearchParameters,
    ) -> tuple[str, ...]: ...
