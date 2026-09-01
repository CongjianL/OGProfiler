#!/usr/bin/env python3
"""Create a bijective FASTA-ID adapter required by historical OGProfiler v1."""
from __future__ import annotations
import argparse,csv
from pathlib import Path

EXTENSIONS={".fa",".faa",".fas",".fasta"}

def fasta_records(path:Path):
 header=None;sequence=[]
 with path.open(encoding="utf-8") as h:
  for line_number,line in enumerate(h,1):
   if line.startswith(">"):
    if header is not None:yield header,sequence
    header=line[1:].strip().split()[0]
    if not header:raise ValueError(f"empty FASTA ID: {path}:{line_number}")
    sequence=[]
   elif header is not None:sequence.append(line.rstrip("\r\n"))
 if header is not None:yield header,sequence

def adapt(input_dir:Path,out_dir:Path,map_path:Path,limit:int|None=None,min_length:int=0)->None:
 files=sorted(p for p in input_dir.iterdir() if p.is_file() and p.suffix.lower() in EXTENSIONS)
 if not files:raise ValueError(f"no FASTA files in {input_dir}")
 out_dir.mkdir(parents=True,exist_ok=False)
 rows=[];original_seen=set();adapted_seen=set()
 for file_index,src in enumerate(files,1):
  dst=out_dir/src.name;written=0
  with dst.open("w",encoding="utf-8",newline="\n") as oh:
   for original,sequence in fasta_records(src):
     if sum(map(len,sequence))<min_length:continue
     if limit is not None and written>=limit:break
     if "|" in original:raise ValueError(f"input ID already contains '|': {original}")
     if original in original_seen:raise ValueError(f"duplicate original FASTA ID: {original}")
     adapted=f"OGPV1_{file_index:02d}|{original}"
     if adapted in adapted_seen:raise ValueError(f"duplicate adapted FASTA ID: {adapted}")
     original_seen.add(original);adapted_seen.add(adapted);rows.append((adapted,original));written+=1
     oh.write(f">{adapted}\n")
     for line in sequence:oh.write(line+"\n")
  if not written:raise ValueError(f"no FASTA records in {src}")
 map_path.parent.mkdir(parents=True,exist_ok=True)
 with map_path.open("w",encoding="utf-8",newline="") as h:
  w=csv.writer(h,delimiter="\t",lineterminator="\n");w.writerow(["adapted_id","original_id"]);w.writerows(rows)

def main()->int:
 p=argparse.ArgumentParser();p.add_argument("--input",type=Path,required=True);p.add_argument("--out",type=Path,required=True);p.add_argument("--map",type=Path,required=True);p.add_argument("--limit-per-file",type=int);p.add_argument("--min-length",type=int,default=0);a=p.parse_args()
 if a.limit_per_file is not None and a.limit_per_file<1:raise ValueError("--limit-per-file must be positive")
 if a.min_length<0:raise ValueError("--min-length must be non-negative")
 adapt(a.input,a.out,a.map,a.limit_per_file,a.min_length);return 0
if __name__=="__main__":raise SystemExit(main())
