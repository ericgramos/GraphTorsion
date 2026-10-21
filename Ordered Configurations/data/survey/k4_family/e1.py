"""e1.py -- which doubling patterns of K_4 carry torsion in H_2(Conf_5)?

Exact integral homology at n = 5 (python/gconf.homology, i.e. the engine in exact mode) of K_4 with each of the eleven
patterns of doubled edges (up to symmetry), and of the doubled K_4 minus an edge.  Appends one JSON line per graph to
e1.jsonl (name, edges, Betti numbers, torsion by degree, seconds); the saved results are doubling_patterns_n5.jsonl.
Run from this folder:  python3 e1.py        (about 5 hours in total; the largest case needs about 5 GB)
"""
import itertools, json, sys, time
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "python"))
from gconf import homology
K4 = list(itertools.combinations(range(4), 2))
pats = [(), ((0, 1),), ((0, 1), (0, 2)), ((0, 1), (2, 3)), ((0, 1), (0, 2), (0, 3)), ((0, 1), (0, 2), (1, 2)),
        ((0, 1), (0, 2), (1, 3)), ((0, 1), (0, 2), (0, 3), (1, 2)), ((0, 1), (0, 2), (1, 3), (2, 3)),
        ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3)), tuple(K4)]
jobs = [("K4 + doubled " + str(list(p)), K4 + list(p)) for p in pats]
jobs.append(("doubled (K4 - e)", [e for e in K4 if e != (2, 3)] * 2))
for name, g in jobs:
    t = time.time()
    r = homology(g, 5, timeout=3600)
    rec = dict(name=name, edges=g, betti=r["betti"], torsion=r["torsion"], secs=round(time.time() - t))
    open("e1.jsonl", "a").write(json.dumps(rec) + "\n")
    print(name, r["betti"], r["torsion"], rec["secs"], flush=True)
print("E1 DONE")
