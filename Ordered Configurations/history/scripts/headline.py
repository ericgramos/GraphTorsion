"""Large single cases, run one at a time. Usage: python3 headline.py name n [maxdeg]
With maxdeg m only d_1..d_m are built: Betti numbers and torsion are then exact for
H_0..H_{m-1} (torsion of H_{m-1} needs d_m)."""
import itertools, json, os, subprocess, sys, time, threading
from gconf import normalize, snf_divisors, gal_euler

K = lambda m: list(itertools.combinations(range(m), 2))
Kb = lambda a, b: [(i, a + j) for i in range(a) for j in range(b)]
GRAPHS = {
    "K4": K(4), "K5": K(5), "K6": K(6), "K7": K(7), "K33": Kb(3, 3), "K34": Kb(3, 4), "K44": Kb(4, 4),
    "Petersen": [(i, (i + 1) % 5) for i in range(5)] + [(i, i + 5) for i in range(5)] + [(5 + i, 5 + (i + 2) % 5) for i in range(5)],
    "theta4": [(0, 1)] * 4, "theta5": [(0, 1)] * 5,
}
name, n = sys.argv[1], int(sys.argv[2])
maxdeg = int(sys.argv[3]) if len(sys.argv) > 3 else -1
V, edges = normalize(GRAPHS[name])
scratch = os.environ["GCONF_SCRATCH"]
prefix = os.path.join(scratch, f"hl_{name}_{n}")
inp = f"{V} {len(edges)} {n}\n" + "".join(f"{u} {v}\n" for u, v in edges)
args = [os.path.join(os.path.dirname(os.path.abspath(__file__)), "oconf4"), "-o", prefix] + (["-d", str(maxdeg)] if maxdeg >= 0 else [])
t0 = time.time()
proc = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
peak = [0]
def watch():
    while proc.poll() is None:
        try:
            with open(f"/proc/{proc.pid}/status") as f:
                for line in f:
                    if line.startswith("VmRSS"):
                        peak[0] = max(peak[0], int(line.split()[1]) // 1024)
        except Exception:
            pass
        time.sleep(1)
threading.Thread(target=watch, daemon=True).start()
out, err = proc.communicate(inp)
if proc.returncode != 0:
    print("FAILED", proc.returncode, err[-500:]); sys.exit(1)
info = json.loads(out)
cells = info["cells"]
top = len(cells) - 1
rank = {0: 0, top + 1: 0}
divs = {}
for m in info["maps"]:
    rr, d = snf_divisors(m["res_file"])
    rank[m["k"]] = m["pivots"] + rr
    divs[m["k"]] = [x for x in d if x != 1]
    os.remove(m["res_file"])
from collections import Counter
_deg = Counter(x for e in edges for x in e)
true_dim = min(n, sum(1 for v in range(V) if _deg[v] >= 2))
complete = (maxdeg < 0) or (maxdeg >= true_dim)
last = top if complete else top - 1          # highest degree whose homology is fully known
betti = [cells[k] - rank[k] - rank[k + 1] for k in range(last + 1)]
torsion = [divs.get(k + 1, []) for k in range(last + 1)]
res = dict(name=name, n=n, edges=edges, cells=cells, betti=betti, torsion=torsion, complete=complete,
           residuals={m["k"]: [m["res_rows"], m["res_cols"], m["maxabs"]] for m in info["maps"]},
           chi=sum((-1) ** k * c for k, c in enumerate(cells)) if complete else None,
           gal=gal_euler(V, edges, n), secs=round(time.time() - t0, 1), peak_mb=peak[0])
print(json.dumps(res))
with open("headline.jsonl", "a") as f:
    f.write(json.dumps(res) + "\n")
