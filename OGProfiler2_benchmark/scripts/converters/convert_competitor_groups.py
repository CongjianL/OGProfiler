#!/usr/bin/env python3
"""Convert frozen primary competitor outputs into the common terminal partition."""
from __future__ import annotations
import argparse,csv,re,sys
from collections import defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent))
from validate_groups import fasta_ids,read_groups,validate

def tokens(cell:str)->list[str]:
 return [x for x in re.split(r"[,;\s]+",cell.strip()) if x and x not in {"*","-"}]

def parse_wide(path:Path,skip_columns:int=1)->dict[str,list[str]]:
 out={}
 with path.open(newline="",encoding="utf-8") as h:
  for i,row in enumerate(csv.reader(h,delimiter="\t")):
   if not row:continue
   if i==0 and (row[0].lower().startswith("orthogroup") or row[0].startswith("#")):continue
   gid=row[0].strip();members=[]
   for cell in row[skip_columns:]:members.extend(tokens(cell))
   if gid and members:out[gid]=members
 return out

def parse_proteinortho(path:Path)->dict[str,list[str]]:
 out={}
 with path.open(newline="",encoding="utf-8") as h:
  reader=csv.reader(h,delimiter="\t");n=0
  for row in reader:
   if not row or row[0].startswith("#"):continue
   n+=1;members=[]
   for cell in row[3:]:members.extend(tokens(cell))
   if members:out[f"PROTEINORTHO_{n:08d}"]=members
 return out

def main()->int:
 p=argparse.ArgumentParser();p.add_argument("--tool",choices=["orthofinder","fastoma","sonicparanoid","proteinortho"],required=True);p.add_argument("--input",type=Path,required=True);p.add_argument("--fasta",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
 if a.tool=="proteinortho":groups=parse_proteinortho(a.input)
 else:groups=parse_wide(a.input)
 expected=fasta_ids(a.fasta);seen={x for members in groups.values() for x in members};unknown=seen-expected
 if unknown:raise ValueError(f"unknown input proteins ({len(unknown)}): {', '.join(sorted(unknown)[:20])}")
 prefix={"orthofinder":"ORTHOFINDER","fastoma":"FASTOMA","sonicparanoid":"SONICPARANOID","proteinortho":"PROTEINORTHO"}[a.tool]
 for pid in sorted(expected-seen):groups[f"{prefix}_UNASSIGNED_SINGLETON_{pid}"]=[pid]
 a.out.parent.mkdir(parents=True,exist_ok=True)
 with a.out.open("w",newline="",encoding="utf-8") as h:
  w=csv.writer(h,delimiter="\t");w.writerow(["group_id","protein_id"])
  for gid in sorted(groups):
   for pid in sorted(groups[gid]):w.writerow([gid,pid])
 errors=validate(read_groups(a.out),expected)
 if errors:raise ValueError("\n".join(errors))
 return 0
if __name__=="__main__":raise SystemExit(main())
