#!/usr/bin/env python3
"""RefOG best-group metrics; never replaces the official benchmark.py score."""
from __future__ import annotations
import argparse,csv,math
from pathlib import Path

def groups(path):
 d={}
 with path.open() as f:
  for r in csv.reader(f,delimiter='\t'):
   if not r or r[0]=='group_id':continue
   d.setdefault(r[0],set()).add(r[1])
 return d
def refogs(path): return {p.stem:set(x.strip() for x in p.read_text().splitlines() if x.strip()) for p in sorted(path.glob('RefOG*.txt'))}
def vi(pred,true,universe):
 n=len(universe); p=[set(v)&universe for v in pred.values()]; t=[set(v)&universe for v in true.values()]; val=0
 for a in p:
  for b in t:
   k=len(a&b)
   if k: val-=2*k/n*math.log((k*k)/(len(a)*len(b)))
 return val
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--refogs',type=Path,required=True);ap.add_argument('--groups',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--summary',type=Path,required=True);a=ap.parse_args();truth=refogs(a.refogs);pred=groups(a.groups); rows=[]
 for name,r in truth.items():
  candidates=[(2*len(r&p)/(len(r)+len(p)),gid,p) for gid,p in pred.items() if r&p] or [(0,'',set())]
  f,g,p=max(candidates);i=len(r&p); rows.append(dict(refog=name,n_true=len(r),best_predicted_group=g,intersection=i,precision=i/len(p) if p else 0,recall=i/len(r),F1=f,exact=str(p==r).lower(),split_count=max(0,sum(bool(r&q) for q in pred.values())-1),contamination=(len(p-r)/len(p)) if p else 0,missing_fraction=len(r-p)/len(r)))
 a.out.parent.mkdir(parents=True,exist_ok=True)
 with a.out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
 macro=sum(float(x['F1']) for x in rows)/len(rows); universe=set().union(*truth.values());a.summary.write_text('metric\tvalue\nmacro_refog_f1\t%s\nexact_recovery_count\t%s\nvariation_of_information_refog_universe\t%s\n' %(macro,sum(x['exact']=='true' for x in rows),vi(pred,truth,universe)))
if __name__=='__main__':main()
