"""v.py -- independent checks of the Z/8 and Z/16 patterns (universal coefficients and the local Smith form).

For K_4 with {01,02,03,12} doubled and with {01,02,03,12,13} doubled, at n = 5: the dimensions of the homology over
Q (a large prime) and over F_2, their difference (which counts the cyclic 2-primary summands in degrees 2 and 3), the
local Smith form mod 16 (orders 2, 4, 8 seen exactly) and the rank deficit of each map against the rational rank
(summands of order >= 16, invisible mod 16).  The printed lines are saved in z8_z16_checks.txt.
Run from this folder:  python3 v.py        (about 30 minutes)
"""
import itertools, json, sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "python"))
from gconf import homology, local_torsion
K4 = list(itertools.combinations(range(4), 2))
for name, extra in [("four doubled {01,02,03,12}", [(0,1),(0,2),(0,3),(1,2)]), ("five doubled", [(0,1),(0,2),(0,3),(1,2),(1,3)])]:
    g = K4 + extra
    q = homology(g, 5, prime=1000003); f2 = homology(g, 5, prime=2)
    l = local_torsion(g, 5, 2, k=4, extra=("-D", "70000"))
    print(name, "| Q dims", q["betti"], "| F2 dims", f2["betti"], "| F2 excess", [a - b for a, b in zip(f2["betti"], q["betti"])],
          "| mod 16 torsion", l["torsion"][2], "| rank deficit (summands of order >= 16)",
          [a - b for a, b in zip(q["rank"], l["rank_local"])], flush=True)
print("V DONE")
