"""z16.py -- the Z/8 and Z/16 patterns recomputed on different complexes for the same spaces.

Each graph is computed exactly over Z (python/gconf.homology) in its original form, after two random relabellings
with random edge orientations and edge orders, and after subdividing an edge (a homeomorphic graph with a different,
larger complex).  All must give the same Betti numbers and torsion.  Appends to z16.jsonl; the saved results are
z8_z16_different_complexes.jsonl.
Run from this folder:  python3 z16.py      (several hours)
"""
import itertools, json, random, sys, time
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "python"))
from gconf import homology
K4 = list(itertools.combinations(range(4), 2))
base = {"Z16": K4 + [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)], "Z8": K4 + [(0, 1), (0, 2), (0, 3), (1, 2)]}
def subdivide(g, idx, new):            # replace edge number idx by a path through a new vertex
    """Replace edge number idx by a path of two edges through a new vertex (a homeomorphic graph)."""
    a, b = g[idx]; return g[:idx] + g[idx + 1:] + [(a, new), (new, b)]
def scramble(g, seed):
    """Random relabelling of the vertices, random edge orientations and a random edge order."""
    rnd = random.Random(seed); V = 1 + max(max(e) for e in g); p = list(range(V)); rnd.shuffle(p)
    h = [(p[a], p[b]) if rnd.random() < 0.5 else (p[b], p[a]) for a, b in g]; rnd.shuffle(h); return h
jobs = []
for nm, g in base.items():
    jobs += [(nm, "original", g), (nm, "relabelled+reoriented seed 1", scramble(g, 1)), (nm, "relabelled+reoriented seed 2", scramble(g, 2)),
             (nm, "single edge 23 subdivided", subdivide(g, 5, 4)), (nm, "one copy of a doubled edge subdivided", subdivide(g, 6, 4))]
for nm, how, g in jobs:
    t = time.time()
    try:
        r = homology(g, 5, timeout=7200)
        rec = dict(pattern=nm, variant=how, edges=g, betti=r["betti"], torsion=r["torsion"], chi_ok=r["chi"] == r["gal"], cells=r["cells"], secs=round(time.time() - t))
    except Exception as e:
        rec = dict(pattern=nm, variant=how, edges=g, error=str(e)[-200:])
    open("z16.jsonl", "a").write(json.dumps(rec) + "\n"); print(rec["pattern"], rec["variant"], rec.get("betti"), rec.get("torsion", rec.get("error")), rec.get("secs"), flush=True)
print("Z16 DONE")
