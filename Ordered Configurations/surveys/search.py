"""search.py -- the exact integral survey (Betti numbers and torsion of every degree).

For every job (graph, n) in a job list whose complex has at most maxcells cells, compute H_*(Conf_n; Z)
exactly with python/gconf.homology (unit-pivot elimination plus exact Smith form of residual blocks), and
append the result as one JSON line.  Torsion, or a disagreement between the Euler characteristic of the
complex and Gal's formula, is reported with "!!!".  Already finished jobs are skipped, so the script can be
restarted.  cell_count(edges, n) is the closed formula
     |S_n(G)| = n! * sum_k e_k(d_v - 1) * C(n - k + E - 1, E - 1)
(e_k the elementary symmetric polynomial in the numbers d_v - 1 over all vertices).
Usage:  python3 search.py jobs.json out.jsonl [maxcells, default 3000000]
Result: data/survey/survey_exact.jsonl (556 cases).
"""
import _common  # sets up the module path and ENGINE (see _common.py)
import json, math, sys, time, os
from gconf import homology, smooth, normalize


def cell_count(edges, n):
    V, edges = normalize(edges)
    deg = [0] * V
    for u, v in edges:
        deg[u] += 1; deg[v] += 1
    E = len(edges)
    # elementary symmetric polynomials of (d_v - 1)
    es = [1]
    for d in deg:
        a = d - 1
        es = [(es[i] if i < len(es) else 0) + (a * es[i - 1] if i >= 1 else 0) for i in range(len(es) + 1)]
    return math.factorial(n) * sum(es[k] * math.comb(n - k + E - 1, E - 1) for k in range(min(n, len(es) - 1) + 1))


def ess(edges):
    V, edges = normalize(edges)
    deg = [0] * V
    for u, v in edges:
        deg[u] += 1; deg[v] += 1
    return sum(1 for d in deg if d >= 3)


if __name__ == "__main__":
    jobs = json.load(open(sys.argv[1]))
    outfn = sys.argv[2]
    maxcells = int(sys.argv[3]) if len(sys.argv) > 3 else 3_000_000
    done = set()
    if os.path.exists(outfn):
        for line in open(outfn):
            r = json.loads(line); done.add((r["name"], r["n"]))
    with open(outfn, "a") as f:
        for job in jobs:
            name, edges, n = job["name"], job["edges"], job["n"]
            if (name, n) in done:
                continue
            cc = cell_count(edges, n)
            if cc > maxcells:
                continue
            t = time.time()
            try:
                r = homology(edges, n, tag="s" + str(os.getpid()))
            except Exception as ex:
                print("FAIL", name, n, ex, flush=True); continue
            r.update(name=name, edges=edges, secs=round(time.time() - t, 2))
            r["residual"] = {str(k): v for k, v in r["residual"].items()}
            tors = any(r["torsion"])
            if tors or r["chi"] != r["gal"]:
                print("!!!", name, n, r["torsion"], r["chi"], r["gal"], flush=True)
            f.write(json.dumps(r) + "\n"); f.flush()
            print(name, n, r["betti"], "tors" if tors else "free", f"{r['secs']}s", flush=True)
