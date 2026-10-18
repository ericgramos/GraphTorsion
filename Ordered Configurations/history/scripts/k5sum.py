import json
from collections import defaultdict
res = defaultdict(dict)
for l in open("k5res.jsonl"):
    r = json.loads(l); res[(tuple(map(tuple, r["doubled"])), r["mode"])][r["prime"]] = r
for (S, mode), d in res.items():
    if 1000003 not in d: continue
    q = d[1000003]["ranks"]; c = d[1000003]["cells"]
    b = [c[k] - ([0]+q)[k] - (q+[0])[k] for k in range(3)]
    flags = {p: [x - y for x, y in zip(q, d[p]["ranks"])] for p in d if p != 1000003}
    tors = {p: {f"H{k}": v for k, v in enumerate(f) if v} for p, f in flags.items()}
    tors = {p: v for p, v in tors.items() if v}
    print(len(S), mode, list(S), "b0..b2", b, "primes", sorted(flags), "TORSION " + str(tors) if tors else "")
