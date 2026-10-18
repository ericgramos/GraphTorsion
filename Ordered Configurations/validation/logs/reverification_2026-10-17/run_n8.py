# K_{3,3} at n = 8: quotient by Sylow-2 subgroup of S_8 (labels, order 128) x order-8 automorphisms, degrees <= 3.
import json, os, sys, time
sys.path.insert(0, "/home/claude/release/gconf-torsion/surveys")
import _common
OUT = "/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad/verify/n8/k33n8.jsonl"
edges = [tuple(map(int, l.split())) for l in open("/home/claude/release/gconf-torsion/examples/graphs/k33_7.in").read().split("\n")[1:] if l.strip()]
L = "1,0,2,3,4,5,6,7;2,3,0,1,4,5,6,7;4,5,6,7,0,1,2,3"
A = "1,0,2,3,4,5;0,1,2,4,3,5;3,4,5,0,1,2"
jobs = [("k33n8_p3_n", ["-q", "3", "2", "-X", "n"]), ("k33n8_p3_l", ["-q", "3", "2", "-X", "l"]),
        ("k33n8_p7_n", ["-q", "7", "1", "-X", "n"]), ("k33n8_p7_l", ["-q", "7", "1", "-X", "l"])]
for name, co in jobs:
    opts = co + ["-g", L, "-a", A, "-F", "-d", "3", "-D", "70000"]
    data, secs, peak, rc = _common.run_engine(opts, edges, 8, tag=name)
    rec = dict(name=name, V=6, n=8, edges=edges, options=opts, rc=rc, secs=round(secs), peak_MB=peak)
    if data:
        rec["cells"] = data["cells"]; rec["maps"] = [dict(k=m["k"], stages=m.get("stages"), pivots=m["pivots"]) for m in data["maps"]]
        if "group" in data: rec["group"] = data["group"]
    with open(OUT, "a") as f: f.write(json.dumps(rec) + "\n")
    print(name, rc, rec.get("maps"), round(secs), peak, flush=True)
print("N8 DONE", flush=True)
