import subprocess,sys
from pathlib import Path
def test_converter(tmp_path):
 run=tmp_path/'run/results';run.mkdir(parents=True)
 (run/'members.tsv').write_text('family_id\toriginal_id\nOG1\tp1\nOG1\tp2\n')
 (run/'hierarchy.tsv').write_text('component_id\tcluster_id\tparent_id\tdepth\tn_genes\n0\t1\t\t0\t2\n')
 out=tmp_path/'groups.tsv'; script=Path(__file__).parents[1]/'scripts/converters/convert_ogprofiler.py'
 subprocess.run([sys.executable,script,'--run',tmp_path/'run','--groups-out',out],check=True)
 assert out.read_text().splitlines()==['group_id\tprotein_id','OG1\tp1','OG1\tp2']
