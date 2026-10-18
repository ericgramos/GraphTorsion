import sys, time
from gconf import homology, subdivide, abrams_homology
g = {"dt": [(0,1),(0,1),(0,2),(0,2),(1,2),(1,2)], "g26": [(0,1),(0,2),(0,2),(1,2),(1,2),(1,2)]}[sys.argv[1]]
n = int(sys.argv[2]); m = int(sys.argv[3])
t = time.time()
a = homology(g, n, tag="adt")
b = abrams_homology(subdivide(g, m), n)
ok = b["betti"][:len(a["betti"])] == a["betti"] and all(x == 0 for x in b["betti"][len(a["betti"]):]) and b["torsion"][:len(a["torsion"])] == a["torsion"]
print(sys.argv[1], n, "Swiatkowski", a["betti"], a["torsion"], "| Abrams (subdiv %d)" % m, b["betti"], b["torsion"], b["cells"], "AGREE" if ok else "DISAGREE", f"{time.time()-t:.0f}s", flush=True)
