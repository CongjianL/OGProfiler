import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parents[1]/'scripts/converters'))
from validate_groups import validate
def test_duplicate_and_missing():
 assert any('duplicate assignment' in x for x in validate({'A':['p1'],'B':['p1']},{'p1','p2'}))
def test_valid_partition(): assert validate({'A':['p1'],'B':['p2']},{'p1','p2'}) == []
