import subprocess, json, time, threading, itertools
S = "/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad"
K33 = [(i,3+j) for i in range(3) for j in range(3)]; K5 = list(itertools.combinations(range(5), 2))
G = {"K33": K33, "K33+de": K33 + [(0,3)], "K5": K5, "K5+dtri": K5 + [(0,1),(1,2),(0,2)]}
def run(args, inp):
    t = time.time(); p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    peak = [0]
    def w():
        while p.poll() is None:
            try: peak[0] = max(peak[0], int(open(f"/proc/{p.pid}/status").read().split("VmRSS:")[1].split()[0]) // 1024)
            except Exception: pass
            time.sleep(1)
    threading.Thread(target=w, daemon=True).start()
    out, err = p.communicate(inp); return round(time.time() - t), peak[0], p.returncode, out, err
def gin(g, n=None):
    V = 1 + max(max(e) for e in g)
    return (f"{V} {len(g)} {n}\n" if n is not None else f"{V} {len(g)}\n") + "".join(f"{a} {b}\n" for a, b in g)
jobs = [("direct", "K5", 6, None, None), ("direct", "K33+de", 6, None, None)]
for kind, name, n, v, k in jobs:
    for p, K in ((2, 4), (3, 3)):
        if kind == "vexp":
            t, m, rc, out, err = run(["./vexp", str(n), str(v), str(k), str(p), str(K)], gin(G[name]))
            if rc == 0:
                d = json.loads(out.strip().splitlines()[-1]); st = d["stages"]
                res = {"rank": st[0], "torsion": {p ** s: c for s, c in enumerate(st) if s > 0 and c}, "leftover": d["leftover_nnz"]}
            else: res = "hypothesis fails (p-torsion in H_{k-1} of Gamma_v at n-1)" if rc == 3 else f"FAILED rc={rc}"
            print(json.dumps({"job": f"{name} n={n} H_{k} via explosion at v={v}", "p": p, "result": res, "secs": t, "peak_MB": m}), flush=True)
        else:
            t, m, rc, out, err = run(["./oconf14", "-q", str(p), str(K), "-g", "1,2,3,4,0,5", "-o", S + "/npq"], gin(G[name], n))
            if rc == 0:
                d = json.loads(out); res = {}
                for mm in d["maps"]:
                    st = mm["stages"]; tt = {p ** s: c for s, c in enumerate(st) if s > 0 and c}
                    if tt or mm["leftover_nnz"]: res[f"H{mm['k']-1}"] = {"torsion": tt, "leftover": mm["leftover_nnz"]}
                res = res or "no p-torsion in any degree"
            else: res = f"FAILED rc={rc}"
            print(json.dumps({"job": f"{name} n={n} all degrees, Z/5 quotient", "p": p, "cells": json.loads(out)["cells"] if rc == 0 else None, "result": res, "secs": t, "peak_MB": m}), flush=True)
