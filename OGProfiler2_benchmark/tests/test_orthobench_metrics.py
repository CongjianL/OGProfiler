import subprocess,sys,csv
from pathlib import Path
def invoke(tmp_path,groups):
 ref=tmp_path/'ref';ref.mkdir();(ref/'RefOG001.txt').write_text('a\nb\n'); (ref/'RefOG002.txt').write_text('c\n')
 g=tmp_path/'g.tsv';g.write_text('group_id\tprotein_id\n'+groups);out=tmp_path/'out.tsv';summary=tmp_path/'s.tsv';s=Path(__file__).parents[1]/'scripts/metrics/orthobench_metrics.py'
 subprocess.run([sys.executable,s,'--refogs',ref,'--groups',g,'--out',out,'--summary',summary],check=True)
 return list(csv.DictReader(out.open(),delimiter='\t'))
def test_perfect(tmp_path): assert float(invoke(tmp_path,'x\ta\nx\tb\ny\tc\n')[0]['F1'])==1
def test_split(tmp_path): assert float(invoke(tmp_path,'x\ta\ny\tb\nz\tc\n')[0]['F1'])==2/3
def test_contamination(tmp_path): assert float(invoke(tmp_path,'x\ta\nx\tb\nx\tc\n')[0]['contamination'])==1/3
