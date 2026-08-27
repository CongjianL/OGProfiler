#!/usr/bin/env python3
"""Read-only stage attribution for the immutable Open Orthobench run."""
from __future__ import annotations
import argparse,csv,math,random,statistics
from collections import Counter,defaultdict
from pathlib import Path
from typing import Iterable
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

NA="NA"
class DSU:
 def __init__(self,items:Iterable[int]): self.p={x:x for x in items}
 def find(self,x:int)->int:
  while self.p[x]!=x:self.p[x]=self.p[self.p[x]];x=self.p[x]
  return x
 def union(self,a:int,b:int)->None:
  a,b=self.find(a),self.find(b)
  if a!=b:self.p[b]=a
 def components(self)->list[set[int]]:
  out:dict[int,set[int]]=defaultdict(set)
  for x in self.p:out[self.find(x)].add(x)
  return list(out.values())

def rows(path:Path)->list[dict[str,object]]: return pq.read_table(path).to_pylist()
def tsv(path:Path,data:list[dict[str,object]])->None:
 path.parent.mkdir(parents=True,exist_ok=True)
 if not data:return
 with path.open("w",newline="",encoding="utf-8") as h:
  w=csv.DictWriter(h,fieldnames=list(data[0]),delimiter="\t");w.writeheader();w.writerows(data)
def read_table(path:Path)->list[dict[str,str]]:
 with path.open(newline="",encoding="utf-8") as h:return list(csv.DictReader(h,delimiter="\t"))

def classify(*,f1:float,recall:float,contamination:float,split_count:int,n_global_components:int)->str:
 catastrophic=recall>=0.8 and contamination>=0.5
 split=split_count>0
 if f1>=0.85:return "WELL_RECOVERED"
 if catastrophic and split:return "MIXED"
 if catastrophic:return "CATASTROPHIC_FUSION"
 if split and n_global_components>1:return "SEARCH_OR_EDGE_DISCONNECTED"
 if split and n_global_components==1:return "HIERARCHY_SPLIT"
 return "MIXED"

def summarize_values(values:list[float])->dict[str,object]:
 return {"n":len(values),"min":min(values) if values else NA,"median":statistics.median(values) if values else NA,"mean":statistics.fmean(values) if values else NA,"max":max(values) if values else NA}

def ranks(values:list[float])->list[float]:
 order=sorted(range(len(values)),key=values.__getitem__);out=[0.0]*len(values);i=0
 while i<len(order):
  j=i+1
  while j<len(order) and values[order[j]]==values[order[i]]:j+=1
  rank=(i+j-1)/2+1
  for k in order[i:j]:out[k]=rank
  i=j
 return out

def spearman(x:list[float],y:list[float])->float:
 if len(x)<3 or len(set(x))<2 or len(set(y))<2:return math.nan
 rx,ry=ranks(x),ranks(y);mx,my=statistics.fmean(rx),statistics.fmean(ry)
 num=sum((a-mx)*(b-my) for a,b in zip(rx,ry,strict=True));den=math.sqrt(sum((a-mx)**2 for a in rx)*sum((b-my)**2 for b in ry))
 return num/den if den else math.nan

def bootstrap_spearman(x:list[float],y:list[float],n:int=2000,seed:int=42)->tuple[float,float,float,int]:
 rho=spearman(x,y);rng=random.Random(seed);boot=[]
 for _ in range(n):
  idx=[rng.randrange(len(x)) for _ in x];v=spearman([x[i] for i in idx],[y[i] for i in idx])
  if not math.isnan(v):boot.append(v)
 boot.sort()
 if not boot:return rho,math.nan,math.nan,0
 return rho,boot[int(.025*(len(boot)-1))],boot[int(.975*(len(boot)-1))],len(boot)

def hierarchy_path(run:Path,component:int,ref_pids:set[int])->tuple[list[dict[str,object]],dict[str,object]|None]:
 d=run/"hierarchy/components"/f"component={component:08d}"
 if not d.is_dir():return [],None
 nodes={int(r["cluster_id"]):r for r in rows(d/"nodes.parquet")}; membership=rows(d/"members.parquet")
 term_by_pid={int(r["protein_id"]):int(r["terminal_cluster_id"]) for r in membership}
 children:dict[int,list[int]]=defaultdict(list)
 for nid,r in nodes.items():
  if r["parent_id"] is not None:children[int(r["parent_id"])].append(nid)
 event_path=run/"evolution/components"/f"component={component:08d}"/"events.parquet"
 events={int(r["cluster_id"]):r for r in rows(event_path)} if event_path.is_file() else {}
 counts=Counter(); selected=ref_pids&set(term_by_pid)
 for pid in selected:
  n=term_by_pid[pid]
  while True:
   counts[n]+=1;parent=nodes[n]["parent_id"]
   if parent is None:break
   n=int(parent)
 first=None
 for nid in sorted(nodes,key=lambda x:(int(nodes[x]["depth"]),x)):
  child_counts=[counts[c] for c in children.get(nid,[]) if counts[c]>0]
  if len(child_counts)>=2:
   first={"first_split_node":f"{component}:{nid}","first_split_depth":nodes[nid]["depth"],"parent_size":nodes[nid]["n_genes"],"children_sizes":",".join(str(nodes[c]["n_genes"]) for c in children[nid]),"members_sent_to_each_child":",".join(str(counts[c]) for c in children[nid])};break
 out=[]
 for nid in sorted(nodes,key=lambda x:(int(nodes[x]["depth"]),x)):
  if counts[nid]==0:continue
  r=nodes[nid];e=events.get(nid,{})
  out.append({"node_id":f"{component}:{nid}","parent_id":NA if r["parent_id"] is None else f"{component}:{r['parent_id']}","depth":r["depth"],"node_size":r["n_genes"],"n_refog_members":counts[nid],"refog_fraction":counts[nid]/len(selected),"children_count":r["child_count"],"event_label_if_available":e.get("network_event",NA),"split_decision_if_available":r.get("split_status",NA) or NA,"stopping_reason_if_available":r.get("terminal_reason",NA) or NA,**(first or {"first_split_node":NA,"first_split_depth":NA,"parent_size":NA,"children_sizes":NA,"members_sent_to_each_child":NA})})
 return out,first

def main()->None:
 p=argparse.ArgumentParser();p.add_argument("--run",type=Path,required=True);p.add_argument("--refogs",type=Path,required=True);p.add_argument("--metrics",type=Path,required=True);p.add_argument("--burden",type=Path,required=True);p.add_argument("--outdir",type=Path,required=True);a=p.parse_args();a.outdir.mkdir(parents=True,exist_ok=True)
 protein_rows=rows(a.run/"input/proteins.parquet");by_name={str(r["original_id"]):int(r["protein_id"]) for r in protein_rows};meta={int(r["protein_id"]):r for r in protein_rows}
 species={int(r["species_id"]):str(r["species_name"]) for r in rows(a.run/"input/species.parquet")}
 truth_names={f.stem:set(f.read_text().splitlines()) for f in sorted(a.refogs.glob("RefOG*.txt"))};truth={k:{by_name[x] for x in v} for k,v in truth_names.items()};pid_ref={pid:k for k,v in truth.items() for pid in v};truth_pids=set(pid_ref)
 metric={r["refog"]:r for r in read_table(a.metrics)}
 component_by_pid={int(r["protein_id"]):int(r["component_id"]) for r in rows(a.run/"components/index.parquet")}
 # Raw-search diagnostic: scan stored hits only; no search is rerun.
 raw_dsu={k:DSU(v) for k,v in truth.items()};raw_hits=Counter();any_hit:dict[str,set[int]]=defaultdict(set)
 truth_array=pa.array(sorted(truth_pids),type=pa.int64());pf=pq.ParquetFile(a.run/"search/hits.parquet")
 for batch in pf.iter_batches(columns=["query_id","target_id"],batch_size=1_000_000):
  mask=pc.or_(pc.is_in(batch.column(0),value_set=truth_array),pc.is_in(batch.column(1),value_set=truth_array));selected=batch.filter(mask)
  qs=selected.column(0).to_pylist();ts=selected.column(1).to_pylist()
  for q,t in zip(qs,ts,strict=True):
   q=int(q);t=int(t)
   if q==t:continue
   if q in truth_pids:any_hit[pid_ref[q]].add(q)
   if t in truth_pids:any_hit[pid_ref[t]].add(t)
   if q in pid_ref and t in pid_ref and pid_ref[q]==pid_ref[t]:
    ref=pid_ref[q];raw_hits[ref]+=1;raw_dsu[ref].union(q,t)
 retained_dsu={k:DSU(v) for k,v in truth.items()};retained_hits=Counter()
 edge_pf=pq.ParquetFile(a.run/"edges/retained_edges.parquet")
 for batch in edge_pf.iter_batches(columns=["u","v"],batch_size=500_000):
  mask=pc.and_(pc.is_in(batch.column(0),value_set=truth_array),pc.is_in(batch.column(1),value_set=truth_array));selected=batch.filter(mask)
  for u,v in zip(selected.column(0).to_pylist(),selected.column(1).to_pylist(),strict=True):
   u=int(u);v=int(v)
   if u in pid_ref and v in pid_ref and pid_ref[u]==pid_ref[v]:
    ref=pid_ref[u];retained_hits[ref]+=1;retained_dsu[ref].union(u,v)
 stage=[];split_origin=[];classes=Counter()
 for ref in sorted(truth):
  members=truth[ref];m=metric[ref];global_counts=Counter(component_by_pid[x] for x in members);rawc=raw_dsu[ref].components();retc=retained_dsu[ref].components()
  terminal_groups=int(m["split_count"])+1;split_count=int(m["split_count"]);recall=float(m["recall"]);cont=float(m["contamination"]);f1=float(m["F1"])
  cls=classify(f1=f1,recall=recall,contamination=cont,split_count=split_count,n_global_components=len(global_counts));classes[cls]+=1
  stage.append({"refog":ref,"true_size":len(members),"n_input":len(members),"n_members_with_any_search_hit":len(any_hit[ref]),"n_raw_within_refog_hits":raw_hits[ref],"raw_search_connected_components":len(rawc),"raw_search_largest_component_fraction":max(map(len,rawc))/len(members),"n_retained_within_refog_edges":retained_hits[ref],"retained_induced_components":len(retc),"retained_largest_induced_component_fraction":max(map(len,retc))/len(members),"n_global_components_containing_refog":len(global_counts),"largest_global_component_refog_fraction":max(global_counts.values())/len(members),"n_terminal_families":terminal_groups,"best_terminal_family":m["best_predicted_group"],"best_terminal_intersection":m["intersection"],"best_terminal_F1":f1,"split_count":split_count,"contamination":cont,"missing_fraction":m["missing_fraction"],"diagnostic_classification":cls})
  if split_count>0:
   hier_split=False
   for comp in global_counts:
    d=a.run/"hierarchy/components"/f"component={comp:08d}"/"members.parquet"
    if d.is_file():
     terms={int(r["terminal_cluster_id"]) for r in rows(d) if int(r["protein_id"]) in members}
     hier_split|=len(terms)>1
   graph_split=len(global_counts)>1
   origin="MIXED_PRE_AND_HIERARCHY" if graph_split and hier_split else "PRE_HIERARCHY_DISCONNECTED" if graph_split else "HIERARCHY_GENERATED" if hier_split else "UNCERTAIN"
   split_origin.append({"refog":ref,"true_size":len(members),"final_split_count":split_count,"n_global_components":len(global_counts),"largest_component_fraction":max(global_counts.values())/len(members),"hierarchy_generated_split":str(hier_split).lower(),"search_or_graph_generated_split":str(graph_split).lower(),"classification":origin})
 tsv(a.outdir/"refog_stage_attribution.tsv",stage);tsv(a.outdir/"refog_error_classification.tsv",[{"classification":k,"n":v,"fraction":v/70} for k,v in sorted(classes.items())]);tsv(a.outdir/"split_origin.tsv",split_origin)
 # Paths for requested and high-split RefOGs.
 requested={"RefOG020","RefOG014","RefOG006","RefOG054","RefOG002","RefOG061"}|{r for r,m in metric.items() if int(m["split_count"])>=5}
 path_root=a.outdir/"hierarchy_paths";path_root.mkdir(exist_ok=True)
 for ref in sorted(requested):
  allrows=[]
  for comp in sorted({component_by_pid[x] for x in truth[ref]}):allrows.extend(hierarchy_path(a.run,comp,truth[ref])[0])
  tsv(path_root/f"{ref}.tsv",allrows)
 # Family/fusion diagnostics.
 result_members=read_table(a.run/"results/members.tsv");fam_pid:dict[str,set[int]]=defaultdict(set)
 for r in result_members:fam_pid[r["family_id"]].add(int(r["protein_id"]))
 family_rows={r["family_id"]:r for r in read_table(a.run/"results/families.tsv")};burden=read_table(a.burden);top=burden[:10]
 catastrophic=[];bridge_rows=[]
 for b in top:
  fid=b["family_id"];members=fam_pid[fid];fr=family_rows[fid];comp=int(fr["component_id"]);cluster=int(fr["cluster_id"])
  comp_edges=a.run/"components/edges"/f"component={comp:08d}"
  erows=pq.read_table(comp_edges).to_pylist() if comp_edges.is_dir() else []
  induced=[r for r in erows if int(r["u"]) in members and int(r["v"]) in members]
  refcounts=Counter(pid_ref[x] for x in members if x in pid_ref);ranked=refcounts.most_common();is_cat=False
  for ref,n in ranked:
   recall=n/len(truth[ref]);contam=(len(members)-n)/len(members);is_cat|=recall>=.8 and contam>=.5
  catastrophic.append({"family_id":fid,"total_size":len(members),"n_species":len({meta[x]["species_id"] for x in members}),"n_refog_members":sum(refcounts.values()),"n_non_refog_members":len(members)-sum(refcounts.values()),"n_refogs":len(refcounts),"top_refog":ranked[0][0] if ranked else NA,"top_refog_fraction":ranked[0][1]/len(members) if ranked else 0,"second_refog":ranked[1][0] if len(ranked)>1 else NA,"official_like_FP":b["official_like_FP"],"fraction_total_FP":b["fraction_of_all_FP"],"hierarchy_depth":next((r["depth"] for r in read_table(a.run/"results/hierarchy.tsv") if r["component_id"]==str(comp) and r["cluster_id"]==str(cluster)),NA),"parent_node":next((r["parent_id"] or NA for r in read_table(a.run/"results/hierarchy.tsv") if r["component_id"]==str(comp) and r["cluster_id"]==str(cluster)),NA),"terminal_reason_if_available":fr["terminal_reason"],"catastrophic_definition_met":str(is_cat).lower()})
  if induced:
   import igraph as ig
   ids=sorted(members);idx={x:i for i,x in enumerate(ids)};g=ig.Graph(n=len(ids),edges=[(idx[int(r["u"])],idx[int(r["v"])]) for r in induced],directed=False);weights=[float(r["weight"]) for r in induced];degree=g.degree();wdegree=g.strength(weights=weights);arts=set(g.articulation_points()) if len(ids)<=10000 else set();between=g.betweenness(cutoff=4,directed=False) if len(ids)<=10000 else [math.nan]*len(ids)
   candidates=sorted(range(len(ids)),key=lambda i:(-(between[i] if not math.isnan(between[i]) else -1),-degree[i],ids[i]))[:20]
   for i in candidates:
    pid=ids[i];bridge_rows.append({"family_id":fid,"protein_id":meta[pid]["original_id"],"species":species[int(meta[pid]["species_id"])],"degree":degree[i],"weighted_degree":wdegree[i],"refog_membership":pid_ref.get(pid,"NON_REFOG"),"bridge_metric":between[i],"articulation_point":str(i in arts).lower()})
 tsv(a.outdir/"catastrophic_fusions.tsv",catastrophic);tsv(a.outdir/"fusion_bridge_candidates.tsv",bridge_rows)
 # Detailed focal family.
 fid="OG000000491";members=fam_pid[fid];detail=a.outdir/fid;detail.mkdir(exist_ok=True);fr=family_rows[fid];comp=int(fr["component_id"])
 tsv(detail/"members.tsv",[{"protein_id":meta[x]["original_id"],"internal_protein_id":x,"species":species[int(meta[x]["species_id"])],"length":meta[x]["length"],"refog":pid_ref.get(x,"NON_REFOG")} for x in sorted(members)])
 sc=Counter(species[int(meta[x]["species_id"])] for x in members);tsv(detail/"species_composition.tsv",[{"species":k,"n":v,"fraction":v/len(members)} for k,v in sorted(sc.items())])
 rc=Counter(pid_ref.get(x,"NON_REFOG") for x in members);tsv(detail/"refog_composition.tsv",[{"category":k,"n":v,"fraction":v/len(members)} for k,v in rc.most_common()])
 lc=Counter(int(meta[x]["length"])//50*50 for x in members);tsv(detail/"sequence_length_distribution.tsv",[{"length_bin_start":k,"n":v} for k,v in sorted(lc.items())])
 path,_=hierarchy_path(a.run,comp,members);tsv(detail/"hierarchy_path.tsv",path)
 erows=pq.read_table(a.run/"components/edges"/f"component={comp:08d}").to_pylist();induced=[r for r in erows if int(r["u"]) in members and int(r["v"]) in members];cats:dict[str,list[dict[str,object]]]=defaultdict(list)
 def cat(pid:int)->str:
  r=pid_ref.get(pid);return r if r in {"RefOG032","RefOG026"} else "OTHER_REFOG" if r else "NON_REFOG"
 for r in induced:
  a1,a2=cat(int(r["u"])),cat(int(r["v"]));key=f"{a1}__{a2}" if a1<=a2 else f"{a2}__{a1}";cats[key].append(r)
 edge_summary=[]
 for k,v in sorted(cats.items()):
  row={"edge_category":k,"n_edges":len(v),"bitscore":NA,"evalue":NA,"identity":NA}
  for field in ("weight","coverage","score_uv","score_vu"):
   vals=[float(x[field]) for x in v];row[f"{field}_mean"]=statistics.fmean(vals);row[f"{field}_median"]=statistics.median(vals)
  edge_summary.append(row)
 tsv(detail/"edge_summary.tsv",edge_summary)
 # Challenge metadata: bundled RefOG membership size is reliable; no new annotation.
 challenge=[{"refog":r,"refog_size":len(truth[r]),"mean_sequence_identity":NA,"evolutionary_rate":NA,"alignment_quality":NA,"domain_number_category":NA,"F1":metric[r]["F1"],"split_count":metric[r]["split_count"],"contamination":metric[r]["contamination"],"missing_fraction":metric[r]["missing_fraction"]} for r in sorted(truth)]
 tsv(a.outdir/"refog_challenge_metadata.tsv",challenge)
 correlations=[]
 for response in ("F1","split_count","contamination","missing_fraction"):
  x=[float(r["refog_size"]) for r in challenge];y=[float(r[response]) for r in challenge];rho,lo,hi,nboot=bootstrap_spearman(x,y)
  correlations.append({"challenge_variable":"refog_size","response":response,"n":len(x),"spearman_rho":rho,"bootstrap_ci_2.5":lo,"bootstrap_ci_97.5":hi,"bootstrap_replicates_valid":nboot,"interpretation":"exploratory"})
 for unavailable in ("mean_sequence_identity","evolutionary_rate","alignment_quality","domain_number_category"):
  for response in ("F1","split_count","contamination","missing_fraction"):
   correlations.append({"challenge_variable":unavailable,"response":response,"n":0,"spearman_rho":NA,"bootstrap_ci_2.5":NA,"bootstrap_ci_97.5":NA,"bootstrap_replicates_valid":0,"interpretation":"NA: no reliable structured field located in bundled Supporting_Data"})
 tsv(a.outdir/"refog_challenge_correlations.tsv",correlations)

if __name__=="__main__":main()
