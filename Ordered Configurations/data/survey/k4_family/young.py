"""young.py -- how many labels does the torsion need, and how do odd label symmetries act?

For the five torsion graphs of the K_4 family at n = 5: exact homology of every partially labelled space
Conf_5(G)/S_lambda (python/gconf.homology with colours; "aabcd" means particles 0 and 1 are indistinguishable), and
the local Smith form mod 16 of the quotients by a 3-cycle and a 5-cycle of the labels (whose coinvariants are the fixed
part of the 2-primary torsion, as the groups have odd order).  Appends to young.jsonl; the saved results are
young_and_odd_symmetry_quotients_n5.jsonl.
Run from this folder:  python3 young.py    (several hours)
"""
import itertools, json, sys, time
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "python"))
from gconf import homology, local_torsion
K4 = list(itertools.combinations(range(4), 2))
G = {"doubled K4": K4 * 2,
     "T222+apex (triangle 012 doubled)": K4 + [(0, 1), (0, 2), (1, 2)],
     "doubled C4 + diagonals": K4 + [(0, 1), (0, 2), (1, 3), (2, 3)],
     "Z/8 pattern (012 and 03 doubled)": K4 + [(0, 1), (0, 2), (0, 3), (1, 2)],
     "Z/16 pattern (all but 23 doubled)": K4 + [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]}
cols = ["aabcd", "aabbc", "aaabc", "aaabb", "aaaab", "aaaaa"]
out = open("young.jsonl", "a")
for name, g in G.items():
    for c in cols:
        t = time.time()
        try:
            r = homology(g, 5, colors=c, timeout=5400)
            rec = dict(graph=name, quotient=c, betti=r["betti"], torsion=r["torsion"], secs=round(time.time() - t))
        except Exception as e:
            rec = dict(graph=name, quotient=c, error=str(e)[-200:])
        out.write(json.dumps(rec) + "\n"); out.flush(); print(rec, flush=True)
for name, g in G.items():
    for lab, gen in [("3-cycle", "1,2,0,3,4"), ("5-cycle", "1,2,3,4,0")]:
        t = time.time()
        r = local_torsion(g, 5, 2, k=4, extra=("-g", gen, "-F", "-D", "70000"))
        rec = dict(graph=name, quotient=lab, torsion_mod16=r["torsion"][2], rank_local=r["rank_local"], secs=round(time.time() - t))
        out.write(json.dumps(rec) + "\n"); out.flush(); print(rec, flush=True)
print("YOUNG DONE")
