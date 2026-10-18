"""unordered_scan.py -- UConf_n of standard graphs with the independent unordered complex (python/unordered.py).
Usage:  python3 unordered_scan.py    (appends to unordered_scan.jsonl; needs an existing file, possibly empty)
"""
import _common  # sets up the module path and ENGINE (see _common.py)
import itertools, json, sys
from unordered import uhomology
K=lambda m: list(itertools.combinations(range(m),2)); Kb=lambda a,b:[(i,a+j) for i in range(a) for j in range(b)]
pet=[(i,(i+1)%5) for i in range(5)]+[(i,i+5) for i in range(5)]+[(5+i,5+(i+2)%5) for i in range(5)]
done={(d['name'],d['n']) for d in map(json.loads,open('unordered_scan.jsonl'))}
out=open('unordered_scan.jsonl','a')
for name,g,N in [('K5',K(5),8),('K33',Kb(3,3),8),('K6',K(6),7),('Pet',pet,6),('K34',Kb(3,4),7),('K4',K(4),8),('K7',K(7),5),('K44',Kb(4,4),6)]:
    for n in range(2,N+1):
        if (name,n) in done: continue
        r=uhomology(g,n)
        t={i:x for i,x in enumerate(r['torsion']) if x}
        out.write(json.dumps(dict(name=name,n=n,betti=r['betti'],torsion=r['torsion']))+'\n'); out.flush()
        print(name,n,r['betti'],{i:dict((p,x.count(p)) for p in set(x)) for i,x in t.items()},flush=True)
