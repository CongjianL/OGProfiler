#!/usr/bin/env python3
"""Parse the official benchmark.py printed precision/recall/F1 values without rewriting them."""
from __future__ import annotations
import argparse,re,csv
from pathlib import Path
PAT=re.compile(r'(?i)\b(precision|recall|f1(?:[- ]score)?)\b[^0-9]*([0-9]+(?:\.[0-9]+)?)')
def main():
 p=argparse.ArgumentParser();p.add_argument('stdout',type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args(); rows=[]
 for line in a.stdout.read_text(errors='replace').splitlines():
  for metric,value in PAT.findall(line): rows.append({'metric':metric.lower().replace(' ','_'),'value':value,'source_line':line})
 if not rows: raise SystemExit('No official scalar metrics parsed; retain raw stdout and update parser against observed format.')
 with a.out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
if __name__=='__main__':main()
