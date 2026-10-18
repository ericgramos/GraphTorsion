# Re-verification of the major computations with the delivered package (fresh build).
import json, os, sys, time, subprocess, itertools
PKG = sys.argv[1]; OUT = sys.argv[2]
sys.path.insert(0, os.path.join(PKG, "python"))
os.environ.setdefault("GCONF_SCRATCH", os.path.dirname(OUT))
import gconf
from gconf import homology, local_torsion, gal_euler, normalize
K4 = list(itertools.combinations(range(4), 2))
C4 = [(0, 1), (1, 2), (2, 3), (3, 0)]
T = lambda a, b, c: [(0, 1)] * a + [(0, 2)] * b + [(1, 2)] * c

def log(rec):
    with open(OUT, "a") as f: f.write(json.dumps(rec) + "\n")
    print(json.dumps(rec)[:400], flush=True)

def peakrun(fn, *a, **k):
    t = time.time(); r = fn(*a, **k); return r, round(time.time() - t, 1)

jobs = []
# exact integral homology
for name, g, n in [("K4", K4, 2), ("K4", K4, 3), ("K4", K4, 4), ("K4", K4, 5), ("K4", K4, 6),
                   ("T(1,2,2)", T(1,2,2), 6), ("T(2,2,2)", T(2,2,2), 5), ("T(1,2,3)", T(1,2,3), 5),
                   ("T(2,2,2)", T(2,2,2), 6), ("T(1,2,3)", T(1,2,3), 6), ("doubled C4", C4 * 2, 4),
                   ("T(1,3,3)", T(1,3,3), 6)]:
    jobs.append(("exact", name, g, n, None))
# field ranks (universal coefficients cross-check)
for name, g, n in [("T(2,2,2)", T(2,2,2), 6), ("T(1,3,3)", T(1,3,3), 6)]:
    for p in (2, 3, 1000003):
        jobs.append(("field", name, g, n, p))
# local Smith forms
for name, g, n, ps in [("doubled K4", K4 * 2, 5, (2, 3, 5)), ("tripled C4", C4 * 3, 5, (2, 3, 5)),
                       ("doubled C4", C4 * 2, 5, (2, 3)), ("doubled K4", K4 * 2, 4, (2, 3)),
                       ("tripled C4", C4 * 3, 4, (2, 3)), ("T(2,2,3)", T(2,2,3), 6, (2, 3))]:
    for p in ps:
        jobs.append(("local", name, g, n, p))
# large-prime Betti numbers for the local cases
for name, g, n in [("doubled K4", K4 * 2, 5), ("tripled C4", C4 * 3, 5), ("doubled C4", C4 * 2, 5)]:
    jobs.append(("field", name, g, n, 1000003))

only = set(sys.argv[3:])
for kind, name, g, n, p in jobs:
    key = f"{kind}:{name}:{n}:{p}"
    if only and key not in only: continue
    try:
        if kind == "exact":
            r, s = peakrun(homology, g, n)
            log(dict(key=key, betti=r["betti"], torsion=r["torsion"], chi=r["chi"], gal=r["gal"], cells=r["cells"], secs=s))
        elif kind == "field":
            r, s = peakrun(homology, g, n, prime=p)
            log(dict(key=key, dims=r["betti"], cells=r["cells"], secs=s))
        else:
            r, s = peakrun(local_torsion, g, n, p, extra=("-D", "70000"))
            lo = [m.get("leftover_nnz", 0) for m in r["maps"]]
            log(dict(key=key, torsion=r["torsion"], leftover=lo, cells=r["cells"], secs=s))
    except Exception as e:
        log(dict(key=key, error=str(e)[-300:]))
