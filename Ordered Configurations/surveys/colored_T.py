"""colored_T.py -- Young-subgroup quotients of Conf_n(T(a,b,c)).
Usage:  python3 colored_T.py 2,2,2 aabcde,aabbcd,aaabcd     (appends to colored_T.jsonl)
Result: every proper quotient of Conf_6(T(2,2,2)) is torsion free; for T(1,3,3) the quotient by one
transposition has 3-torsion exactly Z/3 (data/survey/colored_T.jsonl).
"""
import _common  # sets up the module path and ENGINE (see _common.py)
import json, sys
from gconf import homology
def T(a,b,c):
    """Edge list of the triangle multigraph T(a, b, c)."""
    return [(0,1)]*a+[(0,2)]*b+[(1,2)]*c
abc = tuple(map(int, sys.argv[1].split(","))); 
cols_list = sys.argv[2].split(",")
out = open("colored_T.jsonl", "a")
for cols in cols_list:
    r = homology(T(*abc), len(cols), colors=cols, tag="colT")
    rec = dict(abc=abc, colors=cols, betti=r["betti"], torsion=r["torsion"], chi_ok=r["chi"] == r["gal"])
    out.write(json.dumps(rec) + "\n"); out.flush(); print(rec, flush=True)
