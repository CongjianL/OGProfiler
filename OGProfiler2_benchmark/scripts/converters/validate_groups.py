#!/usr/bin/env python3
"""Strict terminal-family partition validator."""
from __future__ import annotations
import argparse, csv
from collections import defaultdict
from pathlib import Path

def read_groups(path: Path) -> dict[str, list[str]]:
    groups: dict[str, list[str]] = defaultdict(list)
    with path.open(newline='', encoding='utf-8') as f:
        for n,row in enumerate(csv.reader(f, delimiter='\t'), 1):
            if not row or row[0] == 'group_id': continue
            if len(row) != 2 or not row[0] or not row[1]: raise ValueError(f'{path}:{n}: expected group_id<TAB>protein_id')
            groups[row[0]].append(row[1])
    if not groups: raise ValueError(f'{path}: no groups')
    return dict(groups)
def fasta_ids(path: Path) -> set[str]:
    ids=set()
    for p in sorted(path.iterdir() if path.is_dir() else [path]):
        if p.suffix.lower() not in {'.fa','.faa','.fasta','.fna'}: continue
        for line in p.read_text(encoding='utf-8').splitlines():
            if line.startswith('>'): ids.add(line[1:].split()[0])
    return ids
def validate(groups: dict[str,list[str]], expected: set[str]|None=None) -> list[str]:
    errors=[]; assigned={}
    for group,members in groups.items():
        if not members: errors.append(f'empty group: {group}')
        for p in members:
            if p in assigned: errors.append(f'duplicate assignment: {p}: {assigned[p]}, {group}')
            assigned[p]=group
            if expected is not None and p not in expected: errors.append(f'unknown input protein: {p}')
    if expected is not None:
        missing=sorted(expected-set(assigned))
        if missing: errors.append(f'missing input proteins ({len(missing)}): {", ".join(missing[:20])}')
    return errors
def main() -> int:
    ap=argparse.ArgumentParser(); ap.add_argument('groups',type=Path); ap.add_argument('--fasta',type=Path); args=ap.parse_args()
    errors=validate(read_groups(args.groups), fasta_ids(args.fasta) if args.fasta else None)
    if errors: print('\n'.join(errors)); return 1
    print('OK'); return 0
if __name__=='__main__': raise SystemExit(main())
