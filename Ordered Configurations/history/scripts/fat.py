import json, subprocess, sys, time, threading, os
S = "/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad"
C4 = [(0,1),(1,2),(2,3),(3,0)]; K23 = [(i,2+j) for i in range(2) for j in range(3)]
import itertools; K4 = list(itertools.combinations(range(4),2))
jobs = [("C4 tripled", C4*3, 5, True), ("K23 doubled", K23*2, 5, True), ("C4 doubled", C4*2, 5, False),
        ("C4 doubled + diagonal", C4*2+[(0,2)], 5, False), ("C4 doubled", C4*2, 6, True), ("C4 doubled + diagonal", C4*2+[(0,2)], 6, True)]
outfn = "fat.jsonl"
for name, g, n, quot in jobs:
    V = 1 + max(max(e) for e in g)
    inp = f"{V} {len(g)} {n}\n" + "".join(f"{a} {b}\n" for a, b in g)
    perm = ",".join(str(x) for x in ([1,2,3,4,0] + list(range(5, n))))
    rec = dict(name=name, n=n, quotient="Z5" if quot else "none", ranks={})
    for p in ([13, 2, 3] + ([] if quot else [5])):
        args = ["./oconf12", "-p", str(p), "-o", f"{S}/fat"] + (["-g", perm] if quot else [])
        t = time.time(); out = subprocess.run(args, input=inp, capture_output=True, text=True)
        if out.returncode != 0: rec["ranks"][p] = "FAILED"; continue
        d = json.loads(out.stdout); rec["cells"] = d["cells"]; rec["ranks"][p] = [m["pivots"] for m in d["maps"]]
    ref = rec["ranks"].get(13)
    if isinstance(ref, list):
        c = rec["cells"]; r = [0] + ref + [0]; rec["betti"] = [c[k] - r[k] - r[k+1] for k in range(len(c))]
        rec["torsion_signal"] = {p: {f"H{k}": a - b for k, (a, b) in enumerate(zip(ref, v)) if a != b} for p, v in rec["ranks"].items() if p != 13 and isinstance(v, list)}
    with open(outfn, "a") as f: f.write(json.dumps(rec) + "\n")
    print(json.dumps({k: rec[k] for k in rec if k != "ranks"}), flush=True)
