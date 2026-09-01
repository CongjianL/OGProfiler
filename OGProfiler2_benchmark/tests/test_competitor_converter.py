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
 s=out.read_text();assert "SONICPARANOID_00000001\tp1" in s and "SONICPARANOID_00000001\tp3" in s and "0.8" not in s

def test_sonicparanoid_209_compact_columns_and_nonunique_group_ids(tmp_path):
 f=tmp_path/"in.fa";fasta(f);raw=tmp_path/"ortholog_groups.tsv"
 raw.write_text(
  "group_id\tgroup_size\tsp_in_grp\tseed_ortholog_cnt\tA.fa\tB.fa\n"
  "2\t2\t2\t2\tp1\tp2\n"
  "2\t2\t2\t2\tp3\t\n"
 )
 out=tmp_path/"groups.tsv";subprocess.run([sys.executable,SCRIPT,"--tool","sonicparanoid","--input",raw,"--fasta",f,"--out",out],check=True)
 s=out.read_text()
 assert "SONICPARANOID_00000001\tp1" in s and "SONICPARANOID_00000001\tp2" in s
 assert "SONICPARANOID_00000002\tp3" in s

def test_fastoma_roothog_header_is_not_a_group(tmp_path):
 f=tmp_path/"in.fa";fasta(f);raw=tmp_path/"RootHOGs.tsv"
 raw.write_text("RootHOG\tMembers\nHOG:0001\tp1,p2\n")
 out=tmp_path/"groups.tsv";subprocess.run([sys.executable,SCRIPT,"--tool","fastoma","--input",raw,"--fasta",f,"--out",out],check=True)
 assert "HOG:0001\tp1" in out.read_text()

def test_fastoma_header_only_primary_becomes_singletons(tmp_path):
 f=tmp_path/"in.fa";fasta(f);raw=tmp_path/"RootHOGs.tsv"
 raw.write_text("RootHOG\tProtein\tOMAmerRootHOG\n")
 out=tmp_path/"groups.tsv";subprocess.run([sys.executable,SCRIPT,"--tool","fastoma","--input",raw,"--fasta",f,"--out",out],check=True)
 rows=out.read_text().splitlines()
 assert len(rows)==4
 assert all(f"FASTOMA_UNASSIGNED_SINGLETON_p{i}\tp{i}" in rows for i in range(1,4))

def test_fastoma_051_assignment_rows_ignore_omamer_annotation(tmp_path):
 f=tmp_path/"in.fa";fasta(f);raw=tmp_path/"RootHOGs.tsv"
 raw.write_text(
  "RootHOG\tProtein\tOMAmerRootHOG\n"
  "HOG:0000001\tp1\tHOG:F0000010\n"
  "HOG:0000001\tp2\tHOG:F0000010\n"
 )
 out=tmp_path/"groups.tsv";subprocess.run([sys.executable,SCRIPT,"--tool","fastoma","--input",raw,"--fasta",f,"--out",out],check=True)
 rows=out.read_text().splitlines()
 assert "HOG:0000001\tp1" in rows and "HOG:0000001\tp2" in rows
 assert all("HOG:F0000010" not in row for row in rows)
 assert "FASTOMA_UNASSIGNED_SINGLETON_p3\tp3" in rows
