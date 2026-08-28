from __future__ import annotations
import subprocess,sys
from pathlib import Path

SCRIPT=Path(__file__).parents[1]/"scripts/converters/prepare_proteinortho_input.py"

def test_removes_only_stop_markers_from_copy(tmp_path):
 src=tmp_path/"src";src.mkdir();f=src/"a.fa";f.write_text(">p1\nAA*BC\n>p2\nM**K\n")
 out=tmp_path/"out";report=tmp_path/"report.tsv"
 subprocess.run([sys.executable,SCRIPT,"--input",src,"--outdir",out,"--report",report],check=True)
 assert f.read_text()==">p1\nAA*BC\n>p2\nM**K\n"
 assert (out/"a.fa").read_text()==">p1\nAABC\n>p2\nMK\n"
 assert "\t3\t*" in report.read_text()

def test_rejects_unexpected_symbols(tmp_path):
 src=tmp_path/"src";src.mkdir();(src/"a.fa").write_text(">p1\nAA?BC\n")
 r=subprocess.run([sys.executable,SCRIPT,"--input",src,"--outdir",tmp_path/"out","--report",tmp_path/"r.tsv"])
 assert r.returncode!=0
