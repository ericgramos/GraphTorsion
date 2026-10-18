"""headline.py -- large single integral computations, run one at a time (K_5, K_{3,3}, ... at small n).

usage: python3 headline.py NAME n [maxdeg]
  NAME     one of the keys of GRAPHS below
  n        number of particles
  maxdeg   optional: build only d_1 .. d_maxdeg.  Betti numbers and torsion are then exact for
           H_0 .. H_{maxdeg-1}; the torsion of H_{maxdeg-1} needs d_maxdeg, which is built.
Output: one JSON line on stdout, also appended to headline.jsonl in the current directory, with
cells, betti, torsion (per degree, the elementary divisors > 1), the Euler characteristic of the complex
and Gal's formula value (they must agree for a complete run), time and peak memory.
The torsion is computed exactly over Z: unit-pivot elimination in the engine, then the Smith normal form
of the residual block with FLINT (python/gconf.py: snf_divisors).
"""
import json, os, sys, time
from collections import Counter
import _common
from _common import K, Kb, PETERSEN
from gconf import normalize, snf_divisors, gal_euler

GRAPHS = {"K4": K(4), "K5": K(5), "K6": K(6), "K7": K(7), "K33": Kb(3, 3), "K34": Kb(3, 4), "K44": Kb(4, 4),
          "Petersen": PETERSEN, "theta4": [(0, 1)] * 4, "theta5": [(0, 1)] * 5}
name, n = sys.argv[1], int(sys.argv[2])
maxdeg = int(sys.argv[3]) if len(sys.argv) > 3 else -1
V, edges = normalize(GRAPHS[name])
t0 = time.time()
info, secs, peak, rc = _common.run_engine(["-d", str(maxdeg)] if maxdeg >= 0 else [], edges, n, tag=f"hl_{name}_{n}_")
if rc != 0:
    print("FAILED", rc)
    sys.exit(1)
cells = info["cells"]
top = len(cells) - 1
rank = {0: 0, top + 1: 0}
divs = {}
for m in info["maps"]:
    rr, d = snf_divisors(m["res_file"])        # rank and elementary divisors of the residual block
    rank[m["k"]] = m["pivots"] + rr
    divs[m["k"]] = [x for x in d if x != 1]
    os.remove(m["res_file"])
deg = Counter(x for e in edges for x in e)
true_dim = min(n, sum(1 for v in range(V) if deg[v] >= 2))
complete = maxdeg < 0 or maxdeg >= true_dim
last = top if complete else top - 1
betti = [cells[k] - rank[k] - rank[k + 1] for k in range(last + 1)]
torsion = [divs.get(k + 1, []) for k in range(last + 1)]
res = dict(name=name, n=n, edges=edges, cells=cells, betti=betti, torsion=torsion, complete=complete,
           chi=sum((-1) ** k * c for k, c in enumerate(cells)) if complete else None,
           gal=gal_euler(V, edges, n), secs=round(time.time() - t0, 1), peak_mb=peak)
print(json.dumps(res))
with open("headline.jsonl", "a") as f:
    f.write(json.dumps(res) + "\n")
