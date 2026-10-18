"""triangle_family.py -- the triangle multigraphs T(a,b,c) (vertices 0, 1, 2; a, b, c parallel edges
on the sides 01, 02, 12), exact integral homology for the listed (a,b,c) and n.
This is where the first torsion was found: T(2,2,2) and T(1,2,3) at n = 6 (H_2), then T(1,3,3) (odd torsion)
and T(2,2,3).  Output appended to family.jsonl in the current directory (data/survey/family.jsonl).
Usage:  python3 triangle_family.py
Cases already in family.jsonl are skipped.  The first seven jobs are the recorded ones.  The remaining five
(T(2,2,3), T(1,2,4), T(2,3,3), T(0,4,4) at n = 6 and T(1,2,2) at n = 7) were never completed in exact mode:
exact elimination stalls on entry growth (T(2,2,3)) or the complexes are too large for 8 GB.  T(2,2,3) was done
with local Smith forms (local_survey.py) and T(1,2,4) with the shortcut (triangle_shortcut.py).
"""
import _common  # sets up the module path and ENGINE (see _common.py)
import json, sys, time, os
from gconf import homology
def T(a, b, c): return [(0, 1)] * a + [(0, 2)] * b + [(1, 2)] * c
jobs = [((1,2,2),5),((1,2,2),6),((0,3,3),6),((0,2,4),6),((1,2,3),5),((2,2,2),5),
        ((1,3,3),6),((2,2,3),6),((1,2,4),6),((1,2,2),7),((0,4,4),6),((2,3,3),6)]
done = set()
if os.path.exists("family.jsonl"):
    done = {(tuple(json.loads(l)["abc"]), json.loads(l)["n"]) for l in open("family.jsonl")}
out = open("family.jsonl", "a")
for abc, n in jobs:
    if (abc, n) in done: continue
    t = time.time()
    r = homology(T(*abc), n, tag="fam")
    rec = dict(abc=abc, n=n, betti=r["betti"], torsion=r["torsion"], chi_ok=r["chi"] == r["gal"],
               residual={str(k): v for k, v in r["residual"].items()}, secs=round(time.time() - t))
    out.write(json.dumps(rec) + "\n"); out.flush()
    print(rec, flush=True)
