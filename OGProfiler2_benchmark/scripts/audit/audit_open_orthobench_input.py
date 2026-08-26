#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib
from pathlib import Path

def sha(path: Path) -> str:
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1<<20),b''): h.update(block)
 return h.hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--refogs',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 files=sorted(a.input.glob('*.fa'))
 n=sum(sum(1 for line in x.open('rb') if line.startswith(b'>')) for x in files)
 refs=sorted(a.refogs.glob('RefOG*.txt'))
 if (len(files),n,len(refs)) != (12,251378,70): raise SystemExit(f'Unexpected Open Orthobench inputs: files={len(files)}, sequences={n}, refogs={len(refs)}')
 rows=[(str(x),x.stem,sum(1 for line in x.open('rb') if line.startswith(b'>')),sha(x)) for x in files]
 a.out.parent.mkdir(parents=True,exist_ok=True)
 with a.out.open('w') as f:
  f.write('file\tspecies\tn_sequences\tsha256\n')
  for row in rows:f.write('\t'.join(map(str,row))+'\n')
 digest=hashlib.sha256(''.join(f'{r[3]}  {r[0]}\n' for r in rows).encode()).hexdigest()
 (a.out.parent/'OPEN_ORTHOBENCH_INPUT_DIGEST').write_text(digest+'\n')
if __name__=='__main__':main()
