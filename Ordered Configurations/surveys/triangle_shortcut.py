r"""triangle_shortcut.py -- torsion of H_2(Conf_n T(a,b,c)) through the vertex-explosion shortcut.

For each (a, b, c, n): choose the vertex x opposite the side of smallest multiplicity; Gamma \ x is a banana
graph with hairs (two essential vertices).  The script checks the two hypotheses of the shortcut
  (1) H_1(Conf_{n-1}(Gamma \ x)) is torsion free (python/gconf.homology), and
  (2) star x star tori span H_2(Conf_{n-1}(Gamma \ x)) mod p and over Q (analysis/prod2),
and then computes the rank of the relation matrix M with analysis/g2tors for a large prime and for
p = 2, 3, 5, 7; the number of p-torsion summands of H_2(Conf_n Gamma) is rank_Q(M) - rank_p(M).
Usage:  python3 triangle_shortcut.py '[[1,2,3,6],[2,2,2,6]]' out.jsonl
Result: data/survey/triangle_family_shortcut.jsonl and family2.jsonl.
"""
import _common
import json, subprocess, sys, os, time
from _common import T

def g2(edges, x):
    """Gamma minus x: vertex x is removed, and each of its edges becomes a hair ending at a new leaf.
    Returns (number of vertices, edge list); non-x vertices keep their relative order."""
    verts = sorted({v for e in edges for v in e} - {x})
    idx = {v: i for i, v in enumerate(verts)}
    nxt = len(verts)
    out = []
    for (s, t) in edges:
        if s == x:
            out.append((nxt, idx[t])); nxt += 1
        elif t == x:
            out.append((idx[s], nxt)); nxt += 1
        else:
            out.append((idx[s], idx[t]))
    return nxt, out
from gconf import homology
PRIMES = [2, 3, 5, 7]
BIG = 1000003
def gin(edges, V=None):
    V = V or 1 + max(max(e) for e in edges)
    return f"{V} {len(edges)}\n" + "".join(f"{a} {b}\n" for a, b in edges)
def run(cmd, inp):
    out = subprocess.run(cmd, input=inp, capture_output=True, text=True)
    if out.returncode != 0: raise RuntimeError(out.stderr[-300:])
    return json.loads(out.stdout.strip().splitlines()[-1])
cases = json.loads(sys.argv[1]); outfn = sys.argv[2]
done = set()
if os.path.exists(outfn):
    done = {(tuple(json.loads(l)["abc"]), json.loads(l)["n"]) for l in open(outfn)}
for (a, b, c, n) in cases:
    if ((a, b, c), n) in done: continue
    t0 = time.time()
    sides = {2: a, 1: b, 0: c}          # vertex opposite each side: x=2 opp. side a (01), x=1 opp. b (02), x=0 opp. c (12)
    x = min(sides, key=lambda v: sides[v])
    G = T(a, b, c); V2, E2 = g2(G, x)
    rec = dict(abc=[a, b, c], n=n, x=x)
    # validity: H_1(Conf_{n-1}(G\x)) torsion-free; star products span H_2(Conf_{n-1}(G\x)) mod p and over Q
    h = homology(E2, n - 1, tag="f2h")
    rec["h1_torsion_free"] = not any(h["torsion"][1]) if len(h["torsion"]) > 1 else True
    span = {}
    for p in [BIG] + PRIMES:
        r = run([os.path.join(_common.HERE, "..", "analysis", "prod2"), str(n - 1), str(p)], gin(E2, V2))
        span[p] = r["rankP"] == r["b2"]
    rec["stars_span"] = span
    ranks = {}
    for p in [BIG] + PRIMES:
        ranks[p] = run([os.path.join(_common.HERE, "..", "analysis", "g2tors"), str(n), str(x), str(p)], gin(G, 3))["rankM"]
    rec["rankM"] = ranks
    rec["torsion"] = {p: ranks[BIG] - ranks[p] for p in PRIMES}
    rec["valid"] = rec["h1_torsion_free"] and all(span.values())
    rec["secs"] = round(time.time() - t0)
    with open(outfn, "a") as f: f.write(json.dumps(rec) + "\n")
    print(rec, flush=True)
