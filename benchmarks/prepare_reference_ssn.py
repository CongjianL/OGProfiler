#!/usr/bin/env python3
"""Prepare a reference SSN for frozen-V1 hierarchy-only comparison."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import igraph as ig


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ssn", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    destination = args.out / "ssn.gml"
    shutil.copyfile(args.ssn, destination)
    graph = ig.Graph.Read_GML(str(destination))
    if "name" not in graph.vs.attributes():
        raise ValueError("Reference SSN requires a unique vertex name attribute")
    names = [str(value) for value in graph.vs["name"]]
    if len(set(names)) != len(names):
        raise ValueError("Reference SSN vertex names must be unique")
    (args.out / "SequenceIDs.txt").write_text(
        "".join(f"{name}\t{name}\n" for name in names), encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
