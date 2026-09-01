from __future__ import annotations
import subprocess,sys
from pathlib import Path

SCRIPT=Path(__file__).parents[1]/"scripts/preprocess/prepare_ogprofiler_v1_input.py"

def test_bijective_v1_adapter_and_limit(tmp_path):
 src=tmp_path/"source";src.mkdir()
 (src/"B.fa").write_text(">b1 desc\nAA\n>b2\nBB\n")
 (src/"A.fa").write_text(">a1\nCC\n>a2\nDD\n")
 out=tmp_path/"adapted";mapping=tmp_path/"map.tsv"
 subprocess.run([sys.executable,SCRIPT,"--input",src,"--out",out,"--map",mapping,"--limit-per-file","1"],check=True)
 assert (out/"A.fa").read_text()==">OGPV1_01|a1\nCC\n"
 assert (out/"B.fa").read_text()==">OGPV1_02|b1\nAA\n"
 assert mapping.read_text().splitlines()==["adapted_id\toriginal_id","OGPV1_01|a1\ta1","OGPV1_02|b1\tb1"]

def test_adapter_can_select_longer_smoke_records(tmp_path):
 src=tmp_path/"source";src.mkdir();(src/"A.fa").write_text(">short\nAA\n>long\nAAAAA\n")
 out=tmp_path/"adapted";mapping=tmp_path/"map.tsv"
 subprocess.run([sys.executable,SCRIPT,"--input",src,"--out",out,"--map",mapping,"--limit-per-file","1","--min-length","5"],check=True)
 assert (out/"A.fa").read_text()==">OGPV1_01|long\nAAAAA\n"
 assert "short" not in mapping.read_text() and "long" in mapping.read_text()
