from __future__ import annotations
import subprocess,sys
from pathlib import Path

SCRIPT=Path(__file__).parents[1]/"scripts/converters/convert_competitor_groups.py"
def fasta(path):path.write_text(">p1\nAA\n>p2\nAA\n>p3\nAA\n")
def test_orthofinder_and_unassigned(tmp_path):
 f=tmp_path/"in.fa";fasta(f);raw=tmp_path/"Orthogroups.tsv";raw.write_text("Orthogroup\tA\tB\nOG1\tp1, p2\t\n")
 out=tmp_path/"groups.tsv";subprocess.run([sys.executable,SCRIPT,"--tool","orthofinder","--input",raw,"--fasta",f,"--out",out],check=True)
 s=out.read_text();assert "OG1\tp1" in s and "ORTHOFINDER_UNASSIGNED_SINGLETON_p3\tp3" in s
 assert b"\r" not in out.read_bytes()
def test_proteinortho(tmp_path):
 f=tmp_path/"in.fa";fasta(f);raw=tmp_path/"x.tsv";raw.write_text("# Species\tGenes\tAlg.-Conn.\tA\tB\n2\t3\t1\tp1,p2\tp3\n")
 out=tmp_path/"groups.tsv";subprocess.run([sys.executable,SCRIPT,"--tool","proteinortho","--input",raw,"--fasta",f,"--out",out],check=True)
 assert out.read_text().count("\n")==4

def test_sonicparanoid_interleaved_score_columns(tmp_path):
 f=tmp_path/"in.fa";fasta(f);raw=tmp_path/"ortholog_groups.tsv"
 raw.write_text(
  "group_id\tgroup_size\tsp_in_grp\tseed_ortholog_cnt\tA\tavg_score_sp1\tB\tavg_score_sp2\tconflict\n"
  "1\t2\t2\t2\tp1:0.8,p2\t0.8\tp3\t1\tno\n"
 )
 out=tmp_path/"groups.tsv";subprocess.run([sys.executable,SCRIPT,"--tool","sonicparanoid","--input",raw,"--fasta",f,"--out",out],check=True)
 s=out.read_text();assert "SONICPARANOID_1\tp1" in s and "SONICPARANOID_1\tp3" in s and "0.8" not in s

def test_fastoma_roothog_header_is_not_a_group(tmp_path):
 f=tmp_path/"in.fa";fasta(f);raw=tmp_path/"RootHOGs.tsv"
 raw.write_text("RootHOG\tMembers\nHOG:0001\tp1,p2\n")
 out=tmp_path/"groups.tsv";subprocess.run([sys.executable,SCRIPT,"--tool","fastoma","--input",raw,"--fasta",f,"--out",out],check=True)
 assert "HOG:0001\tp1" in out.read_text()
