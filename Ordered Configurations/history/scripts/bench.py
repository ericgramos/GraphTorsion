import subprocess, time, threading, json, itertools, sys
S = "/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad"
def run(args, inp):
    t = time.time(); p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    peak = [0]
    def w():
        while p.poll() is None:
            try: peak[0] = max(peak[0], int(open(f"/proc/{p.pid}/status").read().split("VmRSS:")[1].split()[0]) // 1024)
            except Exception: pass
            time.sleep(0.5)
    threading.Thread(target=w, daemon=True).start()
    out, _ = p.communicate(inp); return time.time() - t, peak[0], p.returncode, out
def gin(g, n):
    V = 1 + max(max(e) for e in g); return f"{V} {len(g)} {n}\n" + "".join(f"{a} {b}\n" for a, b in g)
T = lambda a, b, c: [(0,1)]*a + [(0,2)]*b + [(1,2)]*c
C4 = [(0,1),(1,2),(2,3),(3,0)]; K4 = list(itertools.combinations(range(4), 2))
cases = [("T(2,2,2) n=6", T(2,2,2), 6, 2, []), ("T(1,2,3) n=6", T(1,2,3), 6, 2, []),
         ("tripled C4 n=5 /Z5", C4*3, 5, 3, ["-g", "1,2,3,4,0"]), ("doubled K4 n=5 /Z5", K4*2, 5, 2, ["-g", "1,2,3,4,0"]),
         ("T(1,3,3) n=6", T(1,3,3), 6, 3, [])]
for name, g, n, p, grp in cases:
    inp = gin(g, n); res = {}
    if not grp:   # exact integer computation (original engine, residual SNF negligible here)
        t, m, rc, _ = run(["./oconf3", "-o", S + "/bA"], inp); res["A exact Z"] = (t, m, rc)
    t1, m1, rc1, _ = run(["./oconf9", "-p", "1000003", "-o", S + "/bB"] + grp, inp)
    t2, m2, rc2, _ = run(["./oconf9", "-p", str(p), "-o", S + "/bB"] + grp, inp)
    res["B two field runs"] = (t1 + t2, max(m1, m2), rc1 or rc2)
    oc7 = ["./oconf8"] if grp else ["./oconf7"]
    t, m, rc, _ = run(oc7 + ["-q", str(p), "3" if p == 3 else "4", "-o", S + "/bC"] + grp, inp); res["C old local 16-byte"] = (t, m, rc)
    t, m, rc, _ = run(["./oconf14", "-q", str(p), "3" if p == 3 else "4", "-o", S + "/bD"] + grp, inp); res["D new packed local"] = (t, m, rc)
    print(json.dumps({"case": name, "p": p, **{k: {"secs": round(v[0]), "peak_MB": v[1], "rc": v[2]} for k, v in res.items()}}), flush=True)
