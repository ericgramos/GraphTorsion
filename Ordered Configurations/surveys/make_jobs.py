"""make_jobs.py -- build the job list of the torsion survey (writes jobs.json).

Jobs (graph, n) for n = 3..8: named graphs (K_4..K_7, K_{3,3}, K_{3,4}, K_{4,4}, K_{3,5}, Petersen, Q_3,
prism, Wagner, K_5 - e, Heawood, wheels W_4, W_5); every multigraph of genus 3 and 4 with all degrees >= 3
and at least three essential vertices (from multigraphs.json, produced with python/graphs.py; the file is in
data/survey/); sun graphs (a cycle with one leaf per vertex); and genus 2 and 3 graphs with one leaf added.
Jobs are sorted by the number of cells of the reduced complex (search.cell_count).
Usage:  cd surveys && python3 make_jobs.py      (needs multigraphs.json in the current directory)
"""
import _common  # sets up the module path and ENGINE (see _common.py)
import json, itertools
from search import ess, cell_count
from gconf import smooth, normalize
K = lambda m: list(itertools.combinations(range(m), 2))
Kb = lambda a, b: [(i, a + j) for i in range(a) for j in range(b)]
named = {
 "K4": K(4), "K5": K(5), "K6": K(6), "K7": K(7), "K33": Kb(3,3), "K34": Kb(3,4), "K44": Kb(4,4), "K35": Kb(3,5),
 "Petersen": [(i,(i+1)%5) for i in range(5)] + [(i,i+5) for i in range(5)] + [(5+i,5+(i+2)%5) for i in range(5)],
 "Q3": [(a,b) for a in range(8) for b in range(a+1,8) if bin(a^b).count("1")==1],
 "Prism3": [(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5)],
 "Wagner": [(i,(i+1)%8) for i in range(8)] + [(i,i+4) for i in range(4)],
 "K5-e": [e for e in K(5) if e != (0,1)],
 "Heawood": [(i,(i+1)%14) for i in range(14)] + [(i,(i+5)%14) for i in range(0,14,2)],
 "W4": [(0,i) for i in range(1,5)] + [(i, i%4+1) for i in range(1,5)],
 "W5": [(0,i) for i in range(1,6)] + [(i, i%5+1) for i in range(1,6)],
}
jobs = []
for name, g in named.items():
    for n in range(3, 9):
        jobs.append(dict(name=name, edges=g, n=n))
mg = json.load(open("multigraphs.json"))
for g, lst in mg.items():
    for i, (V, edges) in enumerate(lst):
        if ess(edges) < 3: continue
        for n in range(3, 9):
            jobs.append(dict(name=f"g{g}_{i}", edges=edges, n=n))
# genus 1 and 2 graphs with leaves ("hairs"): attach leaves to a cycle / small graphs
for k in range(3, 7):
    cyc = [(i,(i+1)%k) for i in range(k)] + [(i, k+i) for i in range(k)]
    for n in range(3, 9): jobs.append(dict(name=f"sun{k}", edges=cyc, n=n))
for g in ("2", "3"):
    for i, (V, edges) in enumerate(mg[g]):
        for v in range(V):
            e2 = edges + [(v, V)]
            if ess(e2) < 3: continue
            for n in range(3, 8): jobs.append(dict(name=f"g{g}_{i}+leaf{v}", edges=e2, n=n))
# dedupe and sort by size
jobs.sort(key=lambda j: cell_count(j["edges"], j["n"]))
json.dump(jobs, open("jobs.json", "w"))
print(len(jobs), "jobs;", sum(1 for j in jobs if cell_count(j["edges"], j["n"]) <= 1_000_000), "under 1M cells;",
      sum(1 for j in jobs if cell_count(j["edges"], j["n"]) <= 3_000_000), "under 3M")
