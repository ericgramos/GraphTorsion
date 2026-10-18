import json
from collections import defaultdict
res = defaultdict(dict)
for f in ("k5res.jsonl", "k5res2.jsonl"):
    try:
        for l in open(f):
            r = json.loads(l); res[(tuple(map(tuple, r["doubled"])), r["mode"])][r["prime"]] = r
    except FileNotFoundError: pass
for (S, mode), d in sorted(res.items(), key=lambda kv: (len(kv[0][0]), kv[0][0])):
    refs = [p for p in (1000003, 13, 11, 7) if p in d]
    if not refs: print(len(S), mode, list(S), "no reference yet", sorted(d)); continue
    ref = refs[0]; q = d[ref]["ranks"]; c = d[ref]["cells"]
    b = [c[k] - ([0]+q)[k] - (q+[0])[k] for k in range(3)]
    tors = {}
    for p in d:
        if p == ref: continue
        diff = [x - y for x, y in zip(q, d[p]["ranks"])]
        if any(diff): tors[p] = {f"H{k}": v for k, v in enumerate(diff) if v}
    print(len(S), mode, list(S), "b0..b2", b, "ref", ref, "primes", sorted(p for p in d if p != ref), ("TORSION " + str(tors)) if tors else "")
