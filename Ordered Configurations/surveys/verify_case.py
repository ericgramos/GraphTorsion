"""verify_case.py -- independent re-checks of a torsion case (used for T(2,2,2), T(1,2,3), T(1,3,3)).

usage: python3 verify_case.py GRAPH MODE [n]
  GRAPH   one of the names in GRAPHS below (g4_26 = T(1,2,3), g4_27 = T(2,2,2), T133, and homeomorphic
          variants: one edge subdivided, or every parallel class broken by subdivision into a simple graph;
          these give different cell complexes for the same space)
  MODE    exact      : exact integral computation
          p<prime>   : dimensions over F_prime
          perm<seed> : exact computation after a random relabelling of vertices and edges and random edge
                       orientations (a different complex for the same space)
"""
import random, sys, time
import _common
from gconf import homology

GRAPHS = {"g4_26": [[0, 1], [0, 2], [0, 2], [1, 2], [1, 2], [1, 2]],
          "g4_27": [[0, 1], [0, 1], [0, 2], [0, 2], [1, 2], [1, 2]],
          "T133": [[0, 1], [0, 2], [0, 2], [0, 2], [1, 2], [1, 2], [1, 2]],
          "g4_27sub": [[0, 1], [0, 1], [0, 2], [0, 2], [1, 2], [1, 3], [3, 2]],
          "g4_26sub": [[0, 1], [0, 2], [0, 2], [1, 2], [1, 2], [1, 3], [3, 2]],
          "g4_27simple": [[0, 1], [0, 3], [3, 1], [0, 2], [0, 4], [4, 2], [1, 2], [1, 5], [5, 2]],
          "g4_26simple": [[0, 1], [0, 2], [0, 3], [3, 2], [1, 2], [1, 4], [4, 2], [1, 5], [5, 2]]}
name, mode = sys.argv[1], sys.argv[2]
n = int(sys.argv[3]) if len(sys.argv) > 3 else 6
edges = GRAPHS[name]
t = time.time()
if mode.startswith("p") and mode[1:].isdigit():
    r = homology(edges, n, prime=int(mode[1:]), tag="vc")
    print(name, n, f"F_{mode[1:]}", r["betti"], f"{time.time() - t:.0f}s", flush=True)
elif mode.startswith("perm"):
    rnd = random.Random(int(mode[4:]))
    V = 1 + max(max(e) for e in edges)
    perm = list(range(V))
    rnd.shuffle(perm)
    e2 = [[perm[a], perm[b]] for a, b in edges]
    e2 = [e[::-1] if rnd.random() < 0.5 else e for e in e2]
    rnd.shuffle(e2)
    r = homology(e2, n, tag="vc")
    print(name, n, mode, e2, r["betti"], r["torsion"], r["residual"], f"{time.time() - t:.0f}s", flush=True)
else:
    r = homology(edges, n, tag="vc")
    print(name, n, "exact", r["betti"], r["torsion"], r["residual"], f"{time.time() - t:.0f}s", flush=True)
