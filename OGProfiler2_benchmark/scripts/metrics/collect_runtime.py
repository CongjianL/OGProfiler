#!/usr/bin/env python3
"""Collect run_timed metadata into the stable scalability schema."""
from __future__ import annotations
import argparse,csv,re
from pathlib import Path
FIELDS='dataset tool version replicate seed threads n_species n_proteins n_edges wall_seconds cpu_seconds peak_rss_bytes disk_bytes exit_code status'.split()
def meta(path):
 d={}
 with path.open() as f:
  for r in csv.DictReader(f,delimiter='\t'): d[r['key']]=r['value']
 return d
def seconds(value):
 if not value:return ''
 p=[float(x) for x in value.split(':')]
 return p[0]*3600+p[1]*60+p[2] if len(p)==3 else p[0]*60+p[1] if len(p)==2 else p[0]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('runs',type=Path);ap.add_argument('--out',required=True,type=Path);args=ap.parse_args(); rows=[]
 for p in sorted(args.runs.rglob('metadata.tsv')):
  d=meta(p); version='UNKNOWN'
  command=(p.parent/'command.txt').read_text() if (p.parent/'command.txt').exists() else ''
  m=re.search(r'ogprofiler[^\n]*',command); version='2.0.0a1' if m else version
  rows.append(dict(dataset=d.get('dataset',''),tool=d.get('tool',''),version=version,replicate=d.get('replicate',''),seed=d.get('seed',''),threads=d.get('threads',''),n_species='',n_proteins='',n_edges='',wall_seconds=seconds(d.get('wall_clock','')),cpu_seconds=sum(float(d.get(x,'0') or 0) for x in ('user_cpu_seconds','system_cpu_seconds')),peak_rss_bytes=str(int(d.get('peak_rss_kb','0') or 0)*1024),disk_bytes=d.get('disk_bytes',''),exit_code=d.get('exit_code',''),status=d.get('status','')))
 args.out.parent.mkdir(parents=True,exist_ok=True)
 with args.out.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t');w.writeheader();w.writerows(rows)
if __name__=='__main__':main()
