#!/usr/bin/env python3
"""Convert OGProfiler 2 results exchange tables to strict benchmark TSV tables."""
from __future__ import annotations
import argparse,csv,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent)); from validate_groups import read_groups, validate, fasta_ids

def rows(path: Path):
 with path.open(newline='',encoding='utf-8') as f: yield from csv.DictReader(f,delimiter='\t')
def main() -> int:
 ap=argparse.ArgumentParser(); ap.add_argument('--run',required=True,type=Path); ap.add_argument('--groups-out',required=True,type=Path); ap.add_argument('--hierarchy-dir',type=Path); ap.add_argument('--fasta',type=Path); args=ap.parse_args()
 members=args.run/'results'/'members.tsv'; hierarchy=args.run/'results'/'hierarchy.tsv'
 if not members.is_file(): raise FileNotFoundError(members)
 out={}
 for r in rows(members):
  fid=r.get('family_id'); pid=r.get('original_id')
  if not fid or not pid: raise ValueError('members.tsv requires family_id and original_id')
  out.setdefault(fid,[]).append(pid)
 args.groups_out.parent.mkdir(parents=True,exist_ok=True)
 with args.groups_out.open('w',newline='',encoding='utf-8') as f:
  w=csv.writer(f,delimiter='\t'); w.writerow(['group_id','protein_id'])
  for fid in sorted(out):
   for pid in sorted(out[fid]): w.writerow([fid,pid])
 expected=fasta_ids(args.fasta) if args.fasta else None
 errors=validate(read_groups(args.groups_out),expected)
 if errors: raise ValueError('\n'.join(errors))
 if args.hierarchy_dir:
  args.hierarchy_dir.mkdir(parents=True,exist_ok=True)
  nodes=args.hierarchy_dir/'hierarchy_nodes.tsv'; mem=args.hierarchy_dir/'hierarchy_members.tsv'
  with nodes.open('w',newline='',encoding='utf-8') as f:
   w=csv.writer(f,delimiter='\t'); w.writerow(['node_id','parent_id','level','n_members'])
   for r in rows(hierarchy): w.writerow([f"{r['component_id']}:{r['cluster_id']}", '' if not r['parent_id'] else f"{r['component_id']}:{r['parent_id']}",r['depth'],r['n_genes']])
  with mem.open('w',newline='',encoding='utf-8') as f:
   w=csv.writer(f,delimiter='\t'); w.writerow(['node_id','protein_id'])
   for r in rows(members): w.writerow([r['family_id'],r['original_id']])
 return 0
if __name__=='__main__': raise SystemExit(main())
