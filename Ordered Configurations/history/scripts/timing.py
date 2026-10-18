import itertools, time, sys
from gconf import *
K = lambda m: list(itertools.combinations(range(m), 2))
Kb = lambda a, b: [(i, a + j) for i in range(a) for j in range(b)]
pet = [(i,(i+1)%5) for i in range(5)] + [(i, i+5) for i in range(5)] + [(5+i, 5+(i+2)%5) for i in range(5)]
for name, g, n in [("K5",K(5),3),("K33",Kb(3,3),3),("K5",K(5),4),("K33",Kb(3,3),4),("Pet",pet,3),("K6",K(6),3),("K4",K(4),5)]:
    t=time.time(); a=homology(g,n)
    print(name,n,a["betti"],a["torsion"],a["cells"],"chi",a["chi"]==a["gal"],"resid",a["residual"],f"{time.time()-t:.1f}s",flush=True)
