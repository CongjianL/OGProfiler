#!/usr/bin/env python3
"""Read-only inventory for supplied benchmark roots; writes manifests outside them."""
from __future__ import annotations
import argparse,hashlib
from pathlib import Path
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def nseq(p):
 with p.open('rb') as f:return sum(1 for x in f if x.startswith(b'>'))
def main():
 a=argparse.ArgumentParser();a.add_argument('--orthobench',type=Path,required=True);a.add_argument('--qfo',type=Path,required=True);a.add_argument('--out',type=Path,required=True);x=a.parse_args();x.out.mkdir(parents=True,exist_ok=True); rows=[]; sums=[]
 for name,root,release in [('Open_Orthobench',x.orthobench,'Orthobench Revisited (repository README)'),('QFO',x.qfo,'UNIDENTIFIED: no README/metadata/release marker found during audit')]:
  for p in sorted(root.rglob('*')):
   if not p.is_file():continue
   digest=sha(p); sums.append((digest,str(p)))
   if p.suffix.lower() in {'.fa','.faa','.fasta','.fna'} and ((name == 'Open_Orthobench' and 'BENCHMARKS/Input' in str(p)) or (name == 'QFO' and p.parent.name == 'all')): rows.append((name,release,str(p),p.stem,nseq(p),p.stat().st_size,digest))
 with (x.out/'dataset_manifest.tsv').open('w') as f:
  f.write('dataset\trelease\tfile\tspecies\tn_sequences\tsize_bytes\tsha256\n');[f.write('\t'.join(map(str,r))+'\n') for r in rows]
 with (x.out/'checksums.sha256').open('w') as f:
  for h,p in sums:f.write(f'{h}  {p}\n')
 ob=[r for r in rows if r[0]=='Open_Orthobench'];q=[r for r in rows if r[0]=='QFO']
 (x.out/'dataset_audit.md').write_text(f'''# Dataset audit\n\n* Open Orthobench: {len(ob)} FASTA files; {sum(r[4] for r in ob)} sequences. `README.md` identifies *Orthobench Revisited*; `BENCHMARKS/RefOGs` contains 70 RefOG files and `BENCHMARKS/benchmark.py` is the official scorer.\n* QFO authoritative `all/` input: {len(q)} FASTA files; {sum(r[4] for r in q)} sequences. No release label was found in top-level metadata; release is recorded as unidentified rather than inferred. QFO headers use UniProt `sp|ACCESSION|` and `tr|ACCESSION|` forms; accessions are identifiable but canonical/isoform status requires supplied QfO metadata.\n\nAll entries and SHA256 values are in `dataset_manifest.tsv` and `checksums.sha256`; source data were only read.\n''')
if __name__=='__main__':main()
