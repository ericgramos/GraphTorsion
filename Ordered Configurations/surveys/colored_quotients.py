"""colored_quotients.py -- quotients Conf_n(G)/S_lambda by Young subgroups (positive controls).

A word like "aabc" makes the particles of equal letters indistinguishable; "aa..a" is the unordered space.
Recovers the known torsion of unordered configuration spaces of non-planar graphs (for example
H_1(UConf_2(K_5)) = Z^6 + Z/2) and records the coloured quotients of K_4, K_5, K_{3,3}, K_6.
Usage:  python3 colored_quotients.py               (writes colored.jsonl)
For the triangle graphs (Young quotients of Conf_6 T(2,2,2) and T(1,3,3)) see colored_T.py.
"""
import _common  # sets up the module path and ENGINE (see _common.py)
import itertools, json
from gconf import homology
K=lambda m: list(itertools.combinations(range(m),2)); Kb=lambda a,b:[(i,a+j) for i in range(a) for j in range(b)]
out=open('colored.jsonl','w')
cases=[('K4',K(4),c) for c in ['ab','aab','abc','aabb','aaab','aabc','abcd']]+\
      [('K5',K(5),c) for c in ['ab','aab','abc','aabb','aaab','aabc','abcd','aabcd','abcde']]+\
      [('K33',Kb(3,3),c) for c in ['ab','aab','abc','aabb','aaab','aabc','abcd','aabcd','abcde']]+\
      [('K6',K(6),c) for c in ['aa','ab','aabb','aaab','aabc','abcd','aaaa']]
for name,g,cols in cases:
    r=homology(g,len(cols),colors=cols,tag='col')
    rec=dict(name=name,colors=cols,betti=r['betti'],torsion=r['torsion'],chi_ok=r['chi']==r['gal'])
    out.write(json.dumps(rec)+'\n'); out.flush(); print(rec,flush=True)
