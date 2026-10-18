# Second batch: the two survey cases whose first run reported spurious 2^32 divisors, checked exactly and over fields.
import json, os, sys, time
PKG = sys.argv[1]; OUT = sys.argv[2]
sys.path.insert(0, os.path.join(PKG, "python"))
from gconf import homology
G = {"g4_78": [[0, 1], [0, 1], [0, 2], [1, 4], [2, 4], [2, 4], [3, 3], [3, 4]],
     "g4_80": [[0, 1], [0, 1], [0, 4], [1, 4], [2, 3], [2, 4], [2, 4], [3, 3]]}
for name, g in G.items():
    for p in (0, 2, 3, 1000003):
        t = time.time(); r = homology(g, 5, prime=p)
        rec = dict(key=f"{'exact' if p == 0 else 'field'}:{name}:5:{p}", betti=r["betti"], torsion=r.get("torsion"),
                   residual=r.get("residual"), secs=round(time.time() - t, 1))
        open(OUT, "a").write(json.dumps(rec) + "\n"); print(rec, flush=True)
