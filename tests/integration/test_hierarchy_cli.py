from __future__ import annotations

import json
from pathlib import Path

import igraph as ig
import pyarrow.parquet as pq

from ogprofiler.cli import main


def test_legacy_ssn_to_hierarchy_cli(tmp_path: Path) -> None:
    ssn = tmp_path / "ssn.gml"
    graph = ig.Graph(n=7, edges=[(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5)])
    graph.vs["name"] = [f"legacy_{index}" for index in range(7)]
    graph.es["NBS"] = [1.0] * graph.ecount()
    graph.write_gml(str(ssn))
    run = tmp_path / "run"
    exit_code = main(
        [
            "prototype-hierarchy",
            "--ssn",
            str(ssn),
            "--out",
            str(run),
            "--set",
            "hierarchy.gamma_max=1.0",
        ]
    )
    assert exit_code == 0
    manifest = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["component_count"] == 3
    assert manifest["algorithm_version"] == "hierarchy-prototype-v1"
    assert pq.read_table(run / "hierarchy/components/00000000/members.parquet").num_rows == 3
    assert pq.read_table(run / "hierarchy/components/00000002/members.parquet").num_rows == 1
