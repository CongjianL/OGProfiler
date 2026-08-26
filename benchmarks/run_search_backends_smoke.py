#!/usr/bin/env python3
"""Run the same small prepared dataset through every search backend."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pyarrow.parquet as pq

from ogprofiler.cli import main


def run_backend(proteomes: Path, output_root: Path, backend: str) -> dict[str, object]:
    run_root = output_root / backend
    if main(["prepare", "--proteomes", str(proteomes), "--out", str(run_root)]) != 0:
        raise RuntimeError(f"prepare failed for {backend}")
    if main(["search", "--run", str(run_root), "--backend", backend]) != 0:
        raise RuntimeError(f"search failed for {backend}")
    manifest = json.loads(
        (run_root / "search/search-manifest.json").read_text(encoding="utf-8")
    )
    table = pq.read_table(run_root / "search/hits.parquet")
    return {
        "backend": backend,
        "version": manifest["backend_version"],
        "hit_count": table.num_rows,
        "schema": table.column_names,
    }


def main_script() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--proteomes", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit(f"Output already exists: {args.out}")
    args.out.mkdir(parents=True)
    results = [
        run_backend(args.proteomes, args.out, backend)
        for backend in ("diamond", "mmseqs", "blastp")
    ]
    schemas = {tuple(result["schema"]) for result in results}
    if len(schemas) != 1:
        raise RuntimeError("Search backends emitted different normalized schemas")
    (args.out / "summary.json").write_text(
        json.dumps({"results": results}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"results": results}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main_script())
