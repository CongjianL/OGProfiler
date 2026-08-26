"""Deterministic ID assignment helpers."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path


def sorted_input_paths(paths: Iterable[Path]) -> list[Path]:
    """Return paths in a platform-independent lexical order."""

    return sorted(paths, key=lambda path: path.name.encode("utf-8"))


def sorted_original_ids(original_ids: Iterable[str]) -> list[str]:
    """Return FASTA identifiers in UTF-8 byte order."""

    return sorted(original_ids, key=lambda value: value.encode("utf-8"))
