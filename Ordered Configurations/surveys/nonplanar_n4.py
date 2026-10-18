"""nonplanar_n4.py -- first non-planar tests: doubled K_4 (planar control), doubled K_{3,3}, doubled K_5 and
K_5 with a doubled K_4, at n = 4, by comparing dimensions over F_2, F_3, F_5 with a large prime.
Result: data/nonplanar/earlier/nonplanar4.jsonl (no excess, so no torsion, in every case).
"""
import _common  # sets up the module path and ENGINE (see _common.py)
import itertools, json, time, sys
from gconf import homology
K=lambda m: list(itertools.combinations(range(m),2)); Kb=lambda a,b:[(i,a+j) for i in range(a) for j in range(b)]
d=lambda g: g+g
cases=[('K4 doubled (planar control)', d(K(4))), ('K33 doubled', d(Kb(3,3))), ('K5 doubled', d(K(5))),
       ('K5 with a doubled K4', d(K(4))+[(i,4) for i in range(4)])]
out=open('nonplanar4.jsonl','a')
for name,g in cases:
    t=time.time(); rec={'name':name,'n':4}
    b={p: homology(g,4,prime=p,tag='np4')['betti'] for p in (1000003,2,3,5)}
    rec['betti']=b[1000003]; rec['excess']={p:[x-y for x,y in zip(b[p],b[1000003])] for p in (2,3,5)}; rec['secs']=round(time.time()-t)
    out.write(json.dumps(rec)+'\n'); out.flush(); print(rec, flush=True)
