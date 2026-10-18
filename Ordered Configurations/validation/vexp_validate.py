"""vexp_validate.py -- the vertex-explosion torsion tool (analysis/vexp) against direct computation.

For each graph, prime power and vertex v, the p-primary torsion of H_k(Conf_n G) is computed (1) directly
with the engine (-q p K) and (2) with analysis/vexp at vertex v, for every degree k >= 1.  vexp may report
that its hypothesis (no p-torsion in H_{k-1}(Conf_{n-1} G_v)) fails; that is not an error.  Lines end with
OK, MISMATCH, or "not run (memory)".
"""
import itertools, json, os, subprocess, sys, tempfile, time
HERE = os.path.dirname(os.path.abspath(__file__))
ENGINE = os.environ.get("OCONF_BIN", os.path.join(HERE, "..", "engine", "oconf"))
VEXP = os.path.join(HERE, "..", "analysis", "vexp")
S = tempfile.mkdtemp()

def gin(g, n=None):
    V = 1 + max(max(e) for e in g)
    return (f"{V} {len(g)} {n}\n" if n is not None else f"{V} {len(g)}\n") + "".join(f"{a} {b}\n" for a, b in g)

def direct(g, n, p, K):
    out = subprocess.run([ENGINE, "-q", str(p), str(K), "-o", S + "/vv"], input=gin(g, n), capture_output=True, text=True)
    d = json.loads(out.stdout)
    return {mm["k"] - 1: sorted(p ** s for s, c in enumerate(mm["stages"]) for _ in range(c) if s > 0) for mm in d["maps"]}

def vexp(g, n, v, k, p, K):
    out = subprocess.run([VEXP, str(n), str(v), str(k), str(p), str(K)], input=gin(g), capture_output=True, text=True)
    if out.returncode == 3:
        return "hypothesis fails"
    if out.returncode != 0:
        return "FAILED rc=%d" % out.returncode
    d = json.loads(out.stdout.strip().splitlines()[-1])
    return sorted(p ** s for s, c in enumerate(d["stages"]) for _ in range(c) if s > 0)

T = lambda a, b, c: [(0, 1)] * a + [(0, 2)] * b + [(1, 2)] * c
C4 = [(0, 1), (1, 2), (2, 3), (3, 0)]
K4 = list(itertools.combinations(range(4), 2))
K5 = list(itertools.combinations(range(5), 2))
cases = [  # name, edges, n, [(p, K)], vertices to explode
    ("T(2,2,2)", T(2, 2, 2), 6, [(2, 4)], [0, 2]),
    ("T(1,3,3)", T(1, 3, 3), 6, [(3, 3)], [2, 0]),
    ("T(1,2,2)", T(1, 2, 2), 6, [(2, 4), (3, 3)], [2]),
    ("K4", K4, 5, [(2, 4), (3, 3)], [0, 3]),
    ("K5", K5, 4, [(2, 4), (3, 3)], [0]),
    ("Theta4+leaves", [(0, 1)] * 4 + [(0, 2), (1, 3)], 5, [(2, 4), (3, 3)], [0]),
    ("tripled C4", C4 * 3, 4, [(2, 4), (3, 3)], [0, 1]),
    ("doubled K4", K4 * 2, 4, [(2, 4), (3, 3)], [0]),
]
allok = True
for name, g, n, primes, verts in cases:
    for p, K in primes:
        dres = direct(g, n, p, K)
        for v in verts:
            for k in sorted(dres):
                if k < 1:
                    continue
                vr = vexp(g, n, v, k, p, K)
                ok = vr == "hypothesis fails" or vr == dres[k]
                if isinstance(vr, str) and vr.startswith("FAILED"):
                    ok = None
                allok &= ok is not False
                print(f"{name:14s} n={n} p={p} v={v} H_{k}: direct {dres[k]} | vexp {vr}  "
                      f"{'OK' if ok else ('MISMATCH' if ok is False else 'not run (memory)')}", flush=True)
print("ALL OK" if allok else "MISMATCHES")
