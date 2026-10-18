"""Regression test: compare the release engine with the reference binary it was derived from.

Usage:   python3 regression_vs_reference.py NEW_BINARY REFERENCE_BINARY

Both binaries are run on the same small inputs with many option combinations.  The JSON outputs
(cell counts, ranks, local Smith stages, residual sizes) and, in exact mode, the residual matrices
must agree.  Prints one line per test and "ALL OK" at the end if every test agrees.
"""
import json, os, subprocess, sys, tempfile, itertools

NEW, REF = sys.argv[1], sys.argv[2]
TMP = tempfile.mkdtemp()

def graph(edges, n):
    V = 1 + max(max(e) for e in edges)
    return f"{V} {len(edges)} {n}\n" + "".join(f"{a} {b}\n" for a, b in edges)

def run(binary, args, inp, tag):
    out = subprocess.run([binary, *args, "-o", os.path.join(TMP, tag)], input=inp, capture_output=True, text=True)
    if out.returncode != 0:
        return ("rc", out.returncode)
    d = json.loads(out.stdout)
    # residual files differ only in name; compare their contents instead
    for m in d.get("maps", []):
        fn = m.get("res_file")
        if fn and fn != "none" and os.path.exists(fn):
            m["res_content"] = open(fn).read()
        m.pop("res_file", None)
    if "res_file" in d:
        d["res_content"] = open(d.pop("res_file")).read()
    return d

K4 = list(itertools.combinations(range(4), 2))
T = lambda a, b, c: [(0, 1)] * a + [(0, 2)] * b + [(1, 2)] * c
C4 = [(0, 1), (1, 2), (2, 3), (3, 0)]
K33 = [(i, 3 + j) for i in range(3) for j in range(3)]
T222r = [(0, 1), (1, 2), (2, 0)] * 2
fig8 = [(0, 0), (0, 0)]
lolli = [(0, 1), (1, 1)]
tests = [
    ("K4 n=3 exact", K4, 3, []),
    ("K4 n=4 exact", K4, 4, []),
    ("K4 n=4 p=2", K4, 4, ["-p", "2"]),
    ("K4 n=4 p=3", K4, 4, ["-p", "3"]),
    ("K4 n=4 p=1000003", K4, 4, ["-p", "1000003"]),
    ("T(1,2,3) n=4 exact", T(1, 2, 3), 4, []),
    ("T(1,2,3) n=5 q 2 4", T(1, 2, 3), 5, ["-q", "2", "4"]),
    ("T(1,2,3) n=5 q 3 3", T(1, 2, 3), 5, ["-q", "3", "3"]),
    ("figure eight n=3 exact", fig8, 3, []),
    ("lollipop n=4 exact", lolli, 4, []),
    ("K4 n=3 colours aab", K4, 3, ["-k", "aab"]),
    ("K4 n=4 colours aaaa (unordered)", K4, 4, ["-k", "aaaa"]),
    ("K4 n=4 -d 2", K4, 4, ["-d", "2"]),
    ("K4 n=4 -c", K4, 4, ["-c"]),
    ("K4 n=4 -z", K4, 4, ["-z"]),
    ("K33 n=3 -z with automorphisms", K33, 3, ["-z", "-a", "3,4,5,0,1,2;1,0,2,3,4,5"]),
    ("K4 n=4 -T -p 2", K4, 4, ["-T", "-p", "2"]),
    ("K4 n=4 labels Z3 slow", K4, 4, ["-g", "1,2,0", "-p", "1000003"]),
    ("K4 n=4 labels Z3 fast", K4, 4, ["-g", "1,2,0", "-F", "-p", "1000003"]),
    ("tripled C4 n=4 labels S3 tw=l", C4 * 3, 4, ["-g", "1,0,2,3;1,2,0,3", "-X", "l", "-F", "-q", "5", "2"]),
    ("K33 n=3 side swap + (01), q 3 3", K33, 3, ["-a", "3,4,5,0,1,2;1,0,2,3,4,5", "-F", "-q", "3", "3"]),
    ("K33 n=3 same, tw=le", K33, 3, ["-a", "3,4,5,0,1,2;1,0,2,3,4,5", "-X", "le", "-F", "-q", "3", "3"]),
    ("K33 n=3 same, -Y -1,1", K33, 3, ["-a", "3,4,5,0,1,2;1,0,2,3,4,5", "-Y", "-1,1", "-F", "-q", "3", "3"]),
    ("K33 n=4 diag Z3 x labels Z3, q 2 4", K33, 4, ["-a", "1,2,0,4,5,3", "-g", "1,2,0,3", "-F", "-q", "2", "4"]),
    ("T222 n=5 rotation omega part", T222r, 5, ["-g", "1,2,0,3,4", "-a", "1,2,0", "-W", "1", "-F"]),
    ("T222 n=5 rotation omega part, dense", T222r, 5, ["-g", "1,2,0,3,4", "-a", "1,2,0", "-W", "1", "-F", "-G", "2000"]),
    ("doubled K4 n=4 q 2 2 dense z4", K4 * 2, 4, ["-q", "2", "2", "-D", "3000", "-C", "2"]),
    ("doubled K4 n=4 q 3 2 dense z9", K4 * 2, 4, ["-q", "3", "2", "-D", "3000", "-C", "2"]),
    ("doubled K4 n=4 q 2 4 dense generic", K4 * 2, 4, ["-q", "2", "4", "-D", "3000", "-C", "2"]),
    ("K4 n=4 q 3 3 dense from start", K4, 4, ["-q", "3", "3", "-D", "-1"]),
    ("tripled C4 n=5 /Z5 q 3 2 dense (has Z/3^3)", C4 * 3, 5, ["-g", "1,2,3,4,0", "-F", "-q", "3", "2", "-D", "20000", "-C", "4"]),
    ("T(1,2,3) n=6 /Z5 q 2 2 dense (has Z/2)", T(1, 2, 3), 6, ["-g", "1,2,3,4,0", "-F", "-q", "2", "2", "-D", "20000", "-C", "4"]),
]
allok = True
for name, edges, n, args in tests:
    inp = graph(edges, n)
    a, b = run(NEW, args, inp, "new"), run(REF, args, inp, "ref")
    ok = a == b
    allok &= ok
    print(f"{'OK  ' if ok else 'DIFF'} {name}", flush=True)
    if not ok:
        print("   new:", str(a)[:400]); print("   ref:", str(b)[:400])
# generic matrix mode (-m)
mat = "3 3\n0 0 2\n0 1 1\n1 1 3\n2 2 1\n2 0 4\n"
a = run(NEW, ["-m"], mat, "mn"); b = run(REF, ["-m"], mat, "mr")
print(f"{'OK  ' if a == b else 'DIFF'} generic matrix mode"); allok &= a == b
print("ALL OK" if allok else "FAILURES")
