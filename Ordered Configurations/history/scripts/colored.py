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
