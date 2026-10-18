import json, subprocess, os, sys, time
S = "/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad"
rows = json.load(open("k5sub.json"))
primes = [int(x) for x in sys.argv[1].split(",")]; outfn = sys.argv[2]; maxdoubled = int(sys.argv[3])
done = set()
if os.path.exists(outfn):
    for l in open(outfn): r = json.loads(l); done.add((tuple(map(tuple, r["doubled"])), r["prime"], r["mode"]))
old = {}
for l in open("k5res.jsonl"):
    q = json.loads(l); old.setdefault(tuple(map(tuple, q["doubled"])), set()).add(q["prime"])
for r in rows:
    have = old.get(tuple(map(tuple, r["doubled"])), set())
    if r["ndoubled"] > maxdoubled: continue
    mode = "full" if r["deg_le3"] <= 15_500_000 else "Z5"
    inp = f"5 {len(r['edges'])} 5\n" + "".join(f"{a} {b}\n" for a, b in r["edges"])
    for p in primes:
        key = (tuple(map(tuple, r["doubled"])), p, mode)
        if key in done or p in have: continue
        if p == 13 and 1000003 in have: continue
        t = time.time()
        args = ["./oconf11", "-d", "3", "-p", str(p), "-o", f"{S}/k5r"] + (["-g", "1,2,3,4,0"] if mode == "Z5" else [])
        out = subprocess.run(args, input=inp, capture_output=True, text=True)
        if out.returncode != 0:
            print("FAILED", r["doubled"], p, out.returncode, flush=True); continue
        d = json.loads(out.stdout); ranks = [m["pivots"] for m in d["maps"]]
        for m in d["maps"]:
            try: os.remove(m["res_file"])
            except OSError: pass
        rec = dict(doubled=r["doubled"], mode=mode, prime=p, cells=d["cells"], ranks=ranks, secs=round(time.time() - t))
        with open(outfn, "a") as f: f.write(json.dumps(rec) + "\n")
        print(rec, flush=True)
