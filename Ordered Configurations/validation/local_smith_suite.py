"""local_smith_suite.py -- the local Smith form mode (-q p K) on every known torsion example.

Each line runs the engine with -q p K (optionally with a quotient) and prints the p-primary torsion found
per degree, next to the expected answer.  The expected answers come from independent computations:
UConf_2(K5) and UConf_4(K6) are classical (Abrams; Ko-Park); the T(a,b,c) answers come from exact
integral Smith normal forms (python/gconf.homology) and certificates; doubled K4 and tripled C4 from full
local runs that agree with field ranks.  The T(2,2,3) runs take 25-35 minutes and about 3 GB each;
pass --quick to skip them.
"""
import itertools, json, os, subprocess, sys, tempfile, time
HERE = os.path.dirname(os.path.abspath(__file__))
ENGINE = os.environ.get("OCONF_BIN", os.path.join(HERE, "..", "engine", "oconf"))
S = tempfile.mkdtemp()

def loc(edges, n, p, k, extra=()):
    V = 1 + max(max(e) for e in edges)
    inp = f"{V} {len(edges)} {n}\n" + "".join(f"{a} {b}\n" for a, b in edges)
    out = subprocess.run([ENGINE, "-q", str(p), str(k), "-o", S + "/loc", *extra], input=inp, capture_output=True, text=True)
    if out.returncode:
        return "FAILED %d" % out.returncode
    d = json.loads(out.stdout)
    tors = {}
    for m in d["maps"]:
        st = m["stages"]
        t = [p ** s for s, c in enumerate(st) for _ in range(c) if s > 0]
        if t or m["leftover_nnz"]:
            tors[f"H{m['k'] - 1}"] = {"orders": dict((x, t.count(x)) for x in sorted(set(t))), "leftover": m["leftover_nnz"]}
    return tors

K = lambda m: list(itertools.combinations(range(m), 2))
T = lambda a, b, c: [(0, 1)] * a + [(0, 2)] * b + [(1, 2)] * c
C4 = [(0, 1), (1, 2), (2, 3), (3, 0)]
tests = [("UConf_2(K5)            expect H1: Z/2", K(5), 2, 2, 3, ("-k", "aa")),
         ("UConf_4(K6)            expect H1, H2: Z/2", K(6), 4, 2, 3, ("-k", "aaaa")),
         ("UConf_4(K6) mod 27     expect none", K(6), 4, 3, 3, ("-k", "aaaa")),
         ("doubled K4 n=5 /Z5     expect H2: (Z/2)^8 + Z/4", K(4) * 2, 5, 2, 4, ("-g", "1,2,3,4,0", "-F")),
         ("tripled C4 n=5 /Z5     expect H2: (Z/3)^3", C4 * 3, 5, 3, 3, ("-g", "1,2,3,4,0", "-F")),
         ("T(1,2,3) n=6 /Z5       expect H2: Z/2", T(1, 2, 3), 6, 2, 4, ("-g", "1,2,3,4,0", "-F")),
         ("T(2,2,2) n=6 full      expect H2: (Z/2)^2", T(2, 2, 2), 6, 2, 4, ()),
         ("T(1,2,3) n=6 full      expect H2: Z/2", T(1, 2, 3), 6, 2, 4, ()),
         ("T(1,3,3) n=6 full      expect H2: (Z/3)^6", T(1, 3, 3), 6, 3, 3, ())]
if "--quick" not in sys.argv:
    tests += [("T(2,2,3) n=6 mod 2^4   expect H2: (Z/2)^12", T(2, 2, 3), 6, 2, 4, ()),
              ("T(2,2,3) n=6 mod 3^3   expect H2: (Z/3)^8", T(2, 2, 3), 6, 3, 3, ())]
for name, g, n, p, k, extra in tests:
    t = time.time()
    print(name, "->", loc(g, n, p, k, extra), f"[{time.time() - t:.0f}s]", flush=True)
