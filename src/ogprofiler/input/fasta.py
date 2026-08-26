"""Strict streaming FASTA parsing and residue normalization."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from ogprofiler.exceptions import InputError

VALID_AMINO_ACIDS = frozenset("ABCDEFGHIKLMNPQRSTVWXYZJUO*")


@dataclass(frozen=True, slots=True)
class FastaRecord:
    identifier: str
    sequence: str


def _normalize_sequence(
    sequence: str,
    *,
    path: Path,
    identifier: str,
    illegal_character_policy: str,
) -> str:
    normalized = "".join(sequence.split()).upper()
    if not normalized:
        raise InputError(f"Empty sequence for {identifier!r} in {path}")
    illegal = sorted(set(normalized) - VALID_AMINO_ACIDS)
    if illegal and illegal_character_policy == "error":
        symbols = "".join(illegal)
        raise InputError(f"Illegal sequence characters {symbols!r} for {identifier!r} in {path}")
    if illegal:
        illegal_set = set(illegal)
        normalized = "".join("X" if residue in illegal_set else residue for residue in normalized)
    return normalized


def parse_fasta(path: Path, illegal_character_policy: str = "error") -> Iterator[FastaRecord]:
    identifier: str | None = None
    sequence_parts: list[str] = []
    seen: set[str] = set()
    found_record = False

    try:
        handle = path.open(encoding="utf-8")
    except OSError as error:
        raise InputError(f"Failed to open FASTA {path}: {error}") from error

    with handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    yield FastaRecord(
                        identifier,
                        _normalize_sequence(
                            "".join(sequence_parts),
                            path=path,
                            identifier=identifier,
                            illegal_character_policy=illegal_character_policy,
                        ),
                    )
                header = line[1:].strip()
                identifier = header.split(maxsplit=1)[0] if header else ""
                if not identifier:
                    raise InputError(f"Empty FASTA identifier in {path} at line {line_number}")
                if identifier in seen:
                    raise InputError(f"Duplicate protein ID {identifier!r} in {path}")
                seen.add(identifier)
                sequence_parts = []
                found_record = True
            elif identifier is None:
                raise InputError(
                    f"Sequence data before first FASTA header in {path} at line {line_number}"
                )
            else:
                sequence_parts.append(line)

    if identifier is not None:
        yield FastaRecord(
            identifier,
            _normalize_sequence(
                "".join(sequence_parts),
                path=path,
                identifier=identifier,
                illegal_character_policy=illegal_character_policy,
            ),
        )
    if not found_record:
        raise InputError(f"Empty proteome: {path}")
