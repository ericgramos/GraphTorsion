"""dense_phase_suite.py -- the dense finishing phases agree with pure sparse elimination.

For each test the engine is run twice with the same coefficients: once with the dense phase disabled
(-D 0: sparse elimination only) and once with each of several dense settings that force an early switch:
  "-D -1"           dense from the very first pivot (only feasible for small matrices),
  "-D 3000 -C 2"    switch as soon as <= 3000 columns are alive,
  "-D 20000 -C 8"   switch as soon as <= 20000 columns are alive.
With -q 2 2 the bit-sliced Z/4 phase is used, with -q 3 2 the Z/9 phase, otherwise the byte-per-entry
phase.  All maps (ranks and local Smith stages) must agree.  Prints MATCH/DIFF per test and setting.
"""
import itertools, json, os, subprocess, sys, tempfile, time
HERE = os.path.dirname(os.path.abspath(__file__))
ENGINE = os.environ.get("OCONF_BIN", os.path.join(HERE, "..", "engine", "oconf"))
S = tempfile.mkdtemp()

def run(edges, n, p, k, extra):
    """Run the engine in local mode with extra options (forcing a dense phase) and return its answer and time."""
    V = 1 + max(max(e) for e in edges)
    inp = f"{V} {len(edges)} {n}\n" + "".join(f"{a} {b}\n" for a, b in edges)
    out = subprocess.run([ENGINE, "-q", str(p), str(k), "-o", S + "/dp", *extra], input=inp, capture_output=True, text=True)
    if out.returncode:
        return "FAILED %d" % out.returncode, 0
    d = json.loads(out.stdout)
    return [(m["k"], m["pivots"], m["stages"], m["leftover_nnz"]) for m in d["maps"]], out.stderr.count("switching")

K = lambda m: list(itertools.combinations(range(m), 2))
T = lambda a, b, c: [(0, 1)] * a + [(0, 2)] * b + [(1, 2)] * c
C4 = [(0, 1), (1, 2), (2, 3), (3, 0)]
tests = [("UConf_4(K6) mod 8", K(6), 4, 2, 3, ("-k", "aaaa")),
         ("UConf_4(K6) mod 27", K(6), 4, 3, 3, ("-k", "aaaa")),
         ("K4 n=4 mod 16", K(4), 4, 2, 4, ()),
         ("T(2,2,2) n=4 mod 9", T(2, 2, 2), 4, 3, 2, ()),
         ("doubled K4 n=5 /Z5 mod 4 (Z/4 phase)", K(4) * 2, 5, 2, 2, ("-g", "1,2,3,4,0", "-F")),
         ("doubled K4 n=5 /Z5 mod 16", K(4) * 2, 5, 2, 4, ("-g", "1,2,3,4,0", "-F")),
         ("tripled C4 n=5 /Z5 mod 9 (Z/9 phase)", C4 * 3, 5, 3, 2, ("-g", "1,2,3,4,0", "-F")),
         ("T(1,2,3) n=6 /Z5 mod 4", T(1, 2, 3), 6, 2, 2, ("-g", "1,2,3,4,0", "-F"))]
settings = [("-D", "-1"), ("-D", "3000", "-C", "2"), ("-D", "20000", "-C", "8")]
allok = True
for name, g, n, p, k, extra in tests:
    ref, _ = run(g, n, p, k, tuple(extra) + ("-D", "0"))
    for st in settings:
        if st[1] == "-1" and "/Z5" in name:
            continue                      # dense from the start needs too much memory here
        t = time.time()
        r, sw = run(g, n, p, k, tuple(extra) + st)
        ok = r == ref
        allok &= ok
        print(f"{name:40s} {' '.join(st):14s} {'MATCH' if ok else 'DIFF'} (dense switches: {sw}) [{time.time() - t:.0f}s]", flush=True)
        if not ok:
            print("   sparse:", ref, "\n   dense: ", r)
print("ALL OK" if allok else "FAILURES")
