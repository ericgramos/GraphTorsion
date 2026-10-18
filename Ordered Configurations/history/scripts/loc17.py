import json, subprocess, sys, itertools, time
from loc import S
def loc(binary, edges, n, p, k, extra=()):
    V = 1 + max(max(e) for e in edges)
    inp = f"{V} {len(edges)} {n}\n" + "".join(f"{a} {b}\n" for a, b in edges)
    out = subprocess.run([binary, "-q", str(p), str(k), "-o", S + "/loc17", *extra], input=inp, capture_output=True, text=True)
    if out.returncode: return "FAILED %d %s" % (out.returncode, out.stderr[-300:])
    d = json.loads(out.stdout)
    return [(m["k"], m["pivots"], m["stages"], m["leftover_nnz"]) for m in d["maps"]], out.stderr.count("[dense] switching")
K = lambda m: list(itertools.combinations(range(m), 2))
T = lambda a, b, c: [(0,1)]*a + [(0,2)]*b + [(1,2)]*c
C4 = [(0,1),(1,2),(2,3),(3,0)]
tests = [("UConf_2(K5)", K(5), 2, 2, 3, ("-k", "aa")),
         ("UConf_4(K6)", K(6), 4, 2, 3, ("-k", "aaaa")),
         ("UConf_4(K6) mod 27", K(6), 4, 3, 3, ("-k", "aaaa")),
         ("K4 n=4 mod 16", K(4), 4, 2, 4, ()),
         ("T(2,2,2) n=4 mod 9", T(2,2,2), 4, 3, 2, ()),
         ("doubled K4 n=5 /Z5", K(4)*2, 5, 2, 4, ("-g", "1,2,3,4,0")),
         ("tripled C4 n=5 /Z5", C4*3, 5, 3, 3, ("-g", "1,2,3,4,0")),
         ("T(1,2,3) n=6 /Z5", T(1,2,3), 6, 2, 4, ("-g", "1,2,3,4,0,5"))]
modes = [m.split() for m in sys.argv[1:]]
for name, g, n, p, k, extra in tests:
    ref = loc("./oconf15", g, n, p, k, extra)
    for md in modes:
        t = time.time(); r = loc("./oconf17", g, n, p, k, tuple(extra) + tuple(md))
        ok = isinstance(r, tuple) and r[0] == ref[0]
        print(f"{name:22s} {' '.join(md):14s} {'MATCH' if ok else 'DIFF'} dense_switches={r[1] if isinstance(r, tuple) else '-'} [{time.time()-t:.0f}s]", flush=True)
        if not ok: print("   ref", ref[0], "\n   new", r)
