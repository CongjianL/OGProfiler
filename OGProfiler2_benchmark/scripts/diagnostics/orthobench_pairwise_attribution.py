#!/usr/bin/env python3
"""Official-like Open Orthobench pairwise FP attribution.

The combinatorics reproduce ``BENCHMARKS/benchmark.py:calculate_benchmarks_pairwise``
without modifying the supplied scorer.  Its per-RefOG low-certainty exclusion and
``N = len(RefOG)-1`` normalization are retained.  Contributions are then grouped
by predicted family, a diagnostic decomposition not emitted by the official script.
"""
from __future__ import annotations
import argparse, csv, math
from collections import Counter, defaultdict
from pathlib import Path


def read_groups(path: Path) -> dict[str, set[str]]:
    groups: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            groups[row["group_id"]].add(row["protein_id"])
    return dict(groups)


def read_refogs(root: Path) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    truth = {p.stem: set(p.read_text().splitlines()) for p in sorted(root.glob("RefOG*.txt"))}
    low_root = root / "low_certainty_assignments"
    low = {name: set() for name in truth}
    for p in sorted(low_root.glob("RefOG*.txt")):
        low[p.stem] = set(p.read_text().splitlines())
    return truth, low


def official_family_contribution(reference: set[str], uncertain: set[str], predicted: set[str]) -> tuple[float,float,float,int]:
    """Return normalized TP, FP, family-attributable FN, and overlap.

    Formula source: official Open Orthobench ``calculate_benchmarks_pairwise``.
    The missing-gene FN term has no predicted-family owner and is intentionally
    absent from the family-level FN value.
    """
    ref = reference - uncertain
    pred = predicted - uncertain
    overlap = len(ref & pred)
    if overlap == 0:
        return 0.0, 0.0, 0.0, 0
    denominator = float(len(ref) - 1)
    tp = (overlap * (overlap - 1) / 2.0) / denominator
    fp = (overlap * (len(pred) - overlap)) / denominator
    fn = (overlap * (len(ref) - overlap) / 2.0) / denominator
    return tp, fp, fn, overlap


def attribute(groups: dict[str,set[str]], truth: dict[str,set[str]], low: dict[str,set[str]]) -> list[dict[str,object]]:
    universe = set().union(*truth.values())
    rows: list[dict[str,object]] = []
    for family_id, members in groups.items():
        contributions=[]
        for refog, ref in truth.items():
            tp,fp,fn,overlap=official_family_contribution(ref,low.get(refog,set()),members)
            if overlap: contributions.append((refog,overlap,tp,fp,fn))
        ranked=sorted(contributions,key=lambda x:(-x[1],x[0]))
        rows.append({
            "family_id":family_id,"family_size_total":len(members),
            "n_refog_members":len(members&universe),"n_non_refog_members":len(members-universe),
            "n_distinct_refogs":len(contributions),
            "official_like_TP":sum(x[2] for x in contributions),
            "official_like_FP":sum(x[3] for x in contributions),
            "official_like_FN_if_applicable":sum(x[4] for x in contributions),
            "largest_refog":ranked[0][0] if ranked else "NA",
            "largest_refog_n":ranked[0][1] if ranked else 0,
            "second_refog":ranked[1][0] if len(ranked)>1 else "NA",
            "second_refog_n":ranked[1][1] if len(ranked)>1 else 0,
        })
    rows.sort(key=lambda r:(-float(r["official_like_FP"]),str(r["family_id"])))
    total=sum(float(r["official_like_FP"]) for r in rows)
    cumulative=0.0
    for r in rows:
        fraction=float(r["official_like_FP"])/total if total else 0.0
        cumulative+=fraction;r["fraction_of_all_FP"]=fraction;r["cumulative_FP_fraction"]=cumulative
    return rows


def family_size_rows(groups: dict[str,set[str]]) -> list[dict[str,object]]:
    counts=Counter(map(len,groups.values())); total_f=len(groups);total_p=sum(k*v for k,v in counts.items())
    return [{"family_size":size,"n_families":n,"n_proteins":size*n,
             "fraction_families":n/total_f,"fraction_proteins":size*n/total_p}
            for size,n in sorted(counts.items())]


def write_tsv(path: Path, rows: list[dict[str,object]]) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument("--groups",type=Path,required=True);p.add_argument("--refogs",type=Path,required=True);p.add_argument("--outdir",type=Path,required=True);a=p.parse_args()
    groups=read_groups(a.groups);truth,low=read_refogs(a.refogs);rows=attribute(groups,truth,low)
    write_tsv(a.outdir/"predicted_family_pairwise_burden.tsv",rows)
    write_tsv(a.outdir/"family_size_distribution.tsv",family_size_rows(groups))
    total_fp=sum(float(r["official_like_FP"]) for r in rows)
    with (a.outdir/"pairwise_fp_concentration.tsv").open("w",encoding="utf-8") as h:
        h.write("metric\tvalue\n")
        for n in (1,5,10,20): h.write(f"top_{n}_family_FP_fraction\t{sum(float(r['official_like_FP']) for r in rows[:n])/total_fp}\n")
    sizes=list(map(len,groups.values()))
    bins={"n_terminal_families":len(sizes),"n_singletons":sum(x==1 for x in sizes),"fraction_singletons":sum(x==1 for x in sizes)/len(sizes),"n_size_2":sum(x==2 for x in sizes),"n_size_3_5":sum(3<=x<=5 for x in sizes),"n_size_6_10":sum(6<=x<=10 for x in sizes),"n_size_11_50":sum(11<=x<=50 for x in sizes),"n_size_51_100":sum(51<=x<=100 for x in sizes),"n_size_gt100":sum(x>100 for x in sizes),"n_size_gt500":sum(x>500 for x in sizes),"n_size_gt1000":sum(x>1000 for x in sizes),"largest_family_size":max(sizes)}
    with (a.outdir/"family_size_summary.tsv").open("w",encoding="utf-8") as h:
        h.write("metric\tvalue\n");[h.write(f"{k}\t{v}\n") for k,v in bins.items()]

if __name__=="__main__":main()
