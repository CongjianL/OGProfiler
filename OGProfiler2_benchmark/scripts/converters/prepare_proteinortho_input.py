#!/usr/bin/env python3
"""Create a Proteinortho-compatible copy without changing source FASTA files."""
from __future__ import annotations
import argparse,csv,hashlib,re
from pathlib import Path

ALLOWED=re.compile(r"[^XOUBZACDEFGHIKLMNPQRSTVWYxoubzacdefghiklmnpqrstvwy]")

def digest(path:Path)->str:
 h=hashlib.sha256();h.update(path.read_bytes());return h.hexdigest()

def prepare(src:Path,dst:Path)->tuple[int,str]:
 removed=[];out=[]
 for lineno,line in enumerate(src.read_text().splitlines(),1):
  if line.startswith(">"):
   out.append(line);continue
  bad=ALLOWED.findall(line);removed.extend(bad)
  clean=ALLOWED.sub("",line)
  if clean:out.append(clean)
 dst.write_text("\n".join(out)+"\n")
 unexpected=sorted(set(removed)-{"*"})
 if unexpected:raise ValueError(f"{src}: unexpected removed symbols: {unexpected}")
 return len(removed),"".join(sorted(set(removed)))

def main()->int:
 p=argparse.ArgumentParser();p.add_argument("--input",type=Path,required=True);p.add_argument("--outdir",type=Path,required=True);p.add_argument("--report",type=Path,required=True);a=p.parse_args()
 a.outdir.mkdir(parents=True,exist_ok=False);rows=[]
 for src in sorted(a.input.glob("*.fa")):
  dst=a.outdir/src.name;n,symbols=prepare(src,dst)
  rows.append((src.name,digest(src),digest(dst),n,symbols or "NONE"))
 if not rows:raise ValueError("no input FASTA files")
 a.report.parent.mkdir(parents=True,exist_ok=True)
 with a.report.open("w",newline="") as h:
  w=csv.writer(h,delimiter="\t");w.writerow(["file","source_sha256","prepared_sha256","removed_n","removed_symbols"]);w.writerows(rows)
 return 0
if __name__=="__main__":raise SystemExit(main())
