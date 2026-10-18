import json, subprocess, os, sys, time, threading
S = "/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad"
rows = {tuple(map(tuple, r["doubled"])): r for r in json.load(open("k5sub.json"))}
jobs = json.loads(sys.argv[1]); outfn = sys.argv[2]
for (D, p, T) in jobs:
    r = rows[tuple(map(tuple, D))]
    mode = "full" if r["deg_le3"] <= 15_500_000 else "Z5"
    inp = f"5 {len(r['edges'])} 5\n" + "".join(f"{a} {b}\n" for a, b in r["edges"])
    args = ["./oconf12", "-d", "3", "-p", str(p), "-o", f"{S}/k5r3"] + (["-T"] if T else []) + (["-g", "1,2,3,4,0"] if mode == "Z5" else [])
    t = time.time(); proc = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    peak = [0]
    def watch():
        while proc.poll() is None:
            try: peak[0] = max(peak[0], int(open(f"/proc/{proc.pid}/status").read().split("VmRSS:")[1].split()[0]) // 1024)
            except Exception: pass
            time.sleep(2)
    threading.Thread(target=watch, daemon=True).start()
    out, _ = proc.communicate(inp)
    if proc.returncode != 0: print("FAILED", D, p, T, proc.returncode, "peak", peak[0], flush=True); continue
    d = json.loads(out)
    rec = dict(doubled=D, mode=mode, prime=p, cells=d["cells"], ranks=[m["pivots"] for m in d["maps"]], secs=round(time.time()-t), peak_mb=peak[0], transpose=T)
    with open(outfn, "a") as f: f.write(json.dumps(rec) + "\n")
    print(rec, flush=True)
