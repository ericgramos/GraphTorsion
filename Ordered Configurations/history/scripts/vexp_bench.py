import subprocess, time, threading, json, itertools
S = "/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad"
def run(args, inp):
    t = time.time(); p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    peak = [0]
    def w():
        while p.poll() is None:
            try: peak[0] = max(peak[0], int(open(f"/proc/{p.pid}/status").read().split("VmRSS:")[1].split()[0]) // 1024)
            except Exception: pass
            time.sleep(0.5)
    threading.Thread(target=w, daemon=True).start()
    out, err = p.communicate(inp); return time.time() - t, peak[0], p.returncode, out, err
def gin(g, n=None):
    V = 1 + max(max(e) for e in g)
    return (f"{V} {len(g)} {n}\n" if n is not None else f"{V} {len(g)}\n") + "".join(f"{a} {b}\n" for a, b in g)
C4 = [(0,1),(1,2),(2,3),(3,0)]; K4 = list(itertools.combinations(range(4), 2))
T = lambda a, b, c: [(0,1)]*a + [(0,2)]*b + [(1,2)]*c
cases = [("tripled C4", C4*3, 5, 3, 3), ("doubled K4", K4*2, 5, 2, 4), ("T(2,2,3)", T(2,2,3), 6, 2, 4)]
for name, g, n, p, K in cases:
    t, m, rc, out, _ = run(["./oconf14", "-q", str(p), str(K), "-d", "3", "-o", S + "/vb"], gin(g, n))
    if rc == 0:
        d = json.loads(out); st = [mm["stages"] for mm in d["maps"] if mm["k"] == 3][0]
        dres = sorted(p ** s for s, c in enumerate(st) for _ in range(c) if s > 0)
    else: dres = f"FAILED rc={rc}"
    print(json.dumps({"case": name, "n": n, "p": p, "direct": {"H2_torsion": dres, "secs": round(t), "peak_MB": m}}), flush=True)
    for v in ([0] if name != "T(2,2,3)" else [2]):
        t, m, rc, out, err = run(["./vexp", str(n), str(v), "2", str(p), str(K)], gin(g))
        if rc == 0:
            d = json.loads(out.strip().splitlines()[-1]); vres = sorted(p ** s for s, c in enumerate(d["stages"]) for _ in range(c) if s > 0)
        else: vres = f"FAILED rc={rc}"
        print(json.dumps({"case": name, "v": v, "vexp": {"H2_torsion": vres, "secs": round(t), "peak_MB": m}, "log": err.strip().splitlines()}), flush=True)
