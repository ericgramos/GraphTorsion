"""Independent re-checks of a suspected torsion case.
usage: verify.py <graph> <mode> [n]
  mode = oconf | oconf2 | oconf3 | oconf4 : exact Z computation with that engine build
  mode = p<prime>                          : ranks over F_p
  mode = perm<seed>                        : Z computation after random relabelling of
                                             vertices/edges and random edge orientations"""
import json, sys, os, time, random
from gconf import homology

GRAPHS = {"g4_26": [[0, 1], [0, 2], [0, 2], [1, 2], [1, 2], [1, 2]],
          "g4_27": [[0, 1], [0, 1], [0, 2], [0, 2], [1, 2], [1, 2]],
          "T133": [[0, 1], [0, 2], [0, 2], [0, 2], [1, 2], [1, 2], [1, 2]],
          # homeomorphic copies with one edge subdivided (different cell complex, same space)
          "g4_27sub": [[0, 1], [0, 1], [0, 2], [0, 2], [1, 2], [1, 3], [3, 2]],
          "g4_26sub": [[0, 1], [0, 2], [0, 2], [1, 2], [1, 2], [1, 3], [3, 2]],
          # all parallel classes broken: simple graphs
          "g4_27simple": [[0, 1], [0, 3], [3, 1], [0, 2], [0, 4], [4, 2], [1, 2], [1, 5], [5, 2]],
          "g4_26simple": [[0, 1], [0, 2], [0, 3], [3, 2], [1, 2], [1, 4], [4, 2], [1, 5], [5, 2]]}
name, mode = sys.argv[1], sys.argv[2]
n = int(sys.argv[3]) if len(sys.argv) > 3 else 6
edges = GRAPHS[name]
t = time.time()
tag = f"v{name}{mode}{n}"
if mode.startswith("nodrop-p"):
    os.environ["OCONF_BIN"] = "oconf2"
    r = homology(edges, n, prime=int(mode[8:]), tag=tag)
    print(name, n, f"F_{mode[8:]} (no column dropping)", r["betti"], f"{time.time()-t:.0f}s", flush=True)
elif mode.startswith("p") and mode[1:].isdigit():
    os.environ["OCONF_BIN"] = "oconf4"
    r = homology(edges, n, prime=int(mode[1:]), tag=tag)
    print(name, n, f"F_{mode[1:]}", r["betti"], f"{time.time()-t:.0f}s", flush=True)
elif mode.startswith("perm"):
    rnd = random.Random(int(mode[4:]))
    V = 1 + max(max(e) for e in edges)
    perm = list(range(V)); rnd.shuffle(perm)
    e2 = [[perm[a], perm[b]] for a, b in edges]
    e2 = [e[::-1] if rnd.random() < 0.5 else e for e in e2]
    rnd.shuffle(e2)
    os.environ["OCONF_BIN"] = "oconf4"
    r = homology(e2, n, tag=tag)
    print(name, n, mode, e2, r["betti"], r["torsion"], r["residual"], f"{time.time()-t:.0f}s", flush=True)
else:
    os.environ["OCONF_BIN"] = mode
    r = homology(edges, n, tag=tag)
    print(name, n, mode, r["betti"], r["torsion"], r["residual"], f"{time.time()-t:.0f}s", flush=True)
