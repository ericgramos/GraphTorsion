"""k5_partial_doublings.py -- K_5 with a subset of its edges doubled, n = 5.

Step 1 (--classes): the subsets S of the 10 edges of K_5, up to the action of S_5 on the vertices, give 34
classes; for each, the graph K_5 + S (the edges of S doubled) and the size of its complex in degrees <= 3
(the closed cell-count formula) are written to k5sub.json, sorted by size.
Step 2 (--run JOBS OUT): for each job [S, p, transpose] compute ranks of d_1, d_2, d_3 over F_p (degrees
<= 3, i.e. H_0, H_1, H_2 exactly), on the full complex if it has at most 15.5M cells in degrees <= 3 and
otherwise on the quotient by the label rotation Z/5.  Comparing the ranks for p = 2, 3 with a large prime
(1000003) detects 2- and 3-torsion in H_1 and H_2.
Result (data/nonplanar/earlier/k5res*.jsonl): every class that fit in memory (up to 4 doubled edges) has
identical ranks for p = 2, 3 and 1000003: no 2- or 3-torsion in H_1 or H_2 of Conf_5.
Usage:  python3 k5_partial_doublings.py --classes
        python3 k5_partial_doublings.py --run '[[[[0,1]], 2, false], [[[0,1]], 1000003, false]]' out.jsonl
"""
import itertools, json, math, sys
from collections import Counter
import _common

K5 = list(itertools.combinations(range(5), 2))


def canon(S):
    """Canonical form of a set of doubled edges of K_5 under the symmetric group on the five vertices."""
    return min(tuple(sorted(tuple(sorted((p[a], p[b]))) for a, b in S)) for p in itertools.permutations(range(5)))


def counts(edges, n):            # |C_k| for k = 0..5 by the closed formula
    """Number of cells |C_k| of the complex for k = 0..5, by the closed formula n! e_k(d_v - 1) C(n-k+E-1, E-1)."""
    deg = Counter(x for e in edges for x in e)
    E = len(edges)
    es = [1]
    for v in range(5):
        a = deg[v] - 1
        es = [(es[i] if i < len(es) else 0) + (a * es[i - 1] if i >= 1 else 0) for i in range(len(es) + 1)]
    return [math.factorial(n) * es[k] * math.comb(n - k + E - 1, E - 1) for k in range(min(n, 5) + 1)]


if sys.argv[1] == "--classes":
    classes = {}
    for r in range(0, 11):
        for S in itertools.combinations(K5, r):
            classes.setdefault(canon(S), S)
    rows = []
    for c, S in classes.items():
        g = K5 + list(S)
        cc = counts(g, 5)
        rows.append(dict(doubled=list(S), ndoubled=len(S), deg_le3=sum(cc[:4]), q=sum(cc[:4]) // 5, edges=g))
    rows.sort(key=lambda r: r["q"])
    json.dump(rows, open("k5sub.json", "w"))
    for r in rows:
        print(r["ndoubled"], f"{r['deg_le3'] / 1e6:6.1f}M full  {r['q'] / 1e6:5.1f}M /Z5", r["doubled"])
    print(len(rows), "classes")
else:
    rows = {tuple(map(tuple, r["doubled"])): r for r in json.load(open("k5sub.json"))}
    jobs, outfn = json.loads(sys.argv[2]), sys.argv[3]
    for (Dd, p, Tr) in jobs:
        r = rows[tuple(map(tuple, Dd))]
        mode = "full" if r["deg_le3"] <= 15_500_000 else "Z5"
        opts = ["-d", "3", "-p", str(p)] + (["-T"] if Tr else []) + (["-g", "1,2,3,4,0", "-F"] if mode == "Z5" else [])
        data, secs, peak, rc = _common.run_engine(opts, [tuple(e) for e in r["edges"]], 5, tag="k5")
        if rc != 0:
            print("FAILED", Dd, p, Tr, rc, "peak", peak, flush=True)
            continue
        rec = dict(doubled=Dd, mode=mode, prime=p, cells=data["cells"], ranks=[m["pivots"] for m in data["maps"]],
                   secs=round(secs), peak_mb=peak, transpose=Tr)
        with open(outfn, "a") as f:
            f.write(json.dumps(rec) + "\n")
        print(rec, flush=True)
