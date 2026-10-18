#!/usr/bin/env python3
"""
Is every bit of planar torsion 2-torsion?  A search for a counterexample.

THE CONJECTURE UNDER TEST
-------------------------
Conjecture: if Gamma is planar then every torsion summand of H_*(B_k Gamma) is Z/2.  In
degree one this is a theorem -- Ko--Park: planar graphs have torsion-free H_1, and non-planar
ones have exactly one Z/2.  Above degree one it is an empirical claim, and this script is the
attempt to break it.

WHERE ODD TORSION COULD COME FROM
---------------------------------
The only mechanism we know that produces odd torsion is the inclusion-matrix one: the top
differential of K_{n+1,b} contains the unsigned inclusion matrix W_{n-1,n}(b), and Wilson's
diagonal form supplies the divisor n.  The smallest instance with an odd divisor is K_{4,5},
which gives Z/3 -- and K_{4,5} contains K_{3,3}, so it is not planar.  More generally that
mechanism needs an independent set of size b >= 2n-1 >= 5 attached to n+1 >= 4 vertices, and
K_{4,5} is a minor of any such configuration, so the mechanism is unavailable to planar
graphs.  If this is the only source of odd torsion (prime_census.py), the planar conjecture
follows.

Products cannot help either: Tor(Z/2, Z/2) = Z/2, so Kunneth classes built from 2-torsion stay
2-primary.

So a planar counterexample would need a genuinely new mechanism.  The families below are
chosen as the places where the known mechanism comes closest to being available in a planar
graph -- that is, planar graphs with one or two vertices of large degree and a large
independent set:

    generalised thetas K_{2,m}   two hubs of degree m, m independent branch vertices.  This is
                                 the planar graph that most resembles K_{n+1,b} with a sparse
                                 large side: it IS K_{a,b} with a = 2.  The vanishing for
                                 a <= n applies, so the top homology vanishes, but the lower
                                 degrees are exactly where a new mechanism might hide.
    wheels W_m                   one hub of degree m, attached to an m-cycle.
    bipyramids                   two hubs of degree m, attached to an m-cycle (bipyramid_4 is
                                 the octahedron).
    grids                        the other extreme: planar, 4-regular in the interior, no
                                 high-degree vertex at all.
    triangulations               maximal planar graphs, as dense as planarity allows.
    Platonic solids at k=5       pushing the known carriers to more particles.

RESULT
------
No odd torsion anywhere.  Every torsion summand found in a planar graph, in any degree and at
any particle number reached, is Z/2 -- and the invariant factors are all literally 2, never
4 or 8 either.  The conjecture survives.

    python3 experiments/planar_odd_search.py [--kmax 4] [--cap CELLS]

Exits non-zero if a non-2-primary summand is found in a planar graph, so it doubles as a
regression test on the conjecture.
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import (Graph, antiprism, complete_bip, from_networkx, paths,        # noqa: E402
                 platonic_solids, prism_graph, reduced_cell_count, wheel)
from gbt import fast, lean                                                   # noqa: E402
from gbt.chain import odd_part                                               # noqa: E402
from gbt.report import homology_line, torsion_line                           # noqa: E402


DATA_FILE = "planar_odd_search.json"
FAST_LIMIT = 300_000


def bipyramid(n):
    """Two apexes joined to every vertex of an n-cycle: planar, apex degree n."""
    edges = ([(i, (i + 1) % n) for i in range(n)]
             + [(n, i) for i in range(n)]
             + [(n + 1, i) for i in range(n)])
    return Graph(n + 2, edges)


def triangulation(n, seed):
    """A maximal planar graph on n vertices, grown by stacking into a random face."""
    import random
    import networkx as nx
    G = nx.complete_graph(4)
    rng = random.Random(seed)
    faces = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    for v in range(4, n):
        f = faces.pop(rng.randrange(len(faces)))
        for x in f:
            G.add_edge(v, x)
        faces += [(v, f[0], f[1]), (v, f[0], f[2]), (v, f[1], f[2])]
    return from_networkx(G)


def families():
    """[(name, Graph, note)] -- all planar, chosen as explained in the docstring."""
    import networkx as nx
    out = []
    for m in range(4, 9):
        out.append(("K_{2,%d}" % m, complete_bip(2, m),
                    "generalised theta: two hubs of degree %d" % m))
    for m in range(5, 11):
        out.append(("W_%d" % m, wheel(m), "wheel: one hub of degree %d" % m))
    for m in range(4, 8):
        note = "bipyramid: two hubs of degree %d" % m
        if m == 4:
            note += " (= the octahedron)"
        out.append(("bipyramid_%d" % m, bipyramid(m), note))
    for (a, b) in [(3, 3), (3, 4)]:
        out.append(("grid %dx%d" % (a, b), from_networkx(nx.grid_2d_graph(a, b)),
                    "no high-degree vertex"))
    for n, seed in [(8, 1), (9, 2), (10, 3), (11, 4)]:
        out.append(("triangulation_%d" % n, triangulation(n, seed),
                    "maximal planar, %d vertices" % n))
    for name, G in platonic_solids():
        out.append((name, G, "Platonic"))
    for n in (5, 6):
        out.append(("prism_%d" % n, prism_graph(n), "planar 3-regular"))
        out.append(("antiprism_%d" % n, antiprism(n), "planar 4-regular"))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kmax", type=int, default=4)
    ap.add_argument("--cap", type=int, default=600_000)
    args = ap.parse_args()

    import networkx as nx
    res = paths.load_json(DATA_FILE)
    offenders = []
    n_tested = 0
    exponents = set()

    hdr = ("%-18s %3s %3s %3s %11s %-46s %s"
           % ("graph", "N", "M", "k", "cells", "homology", "torsion"))
    print(hdr)
    print("-" * len(hdr))
    for name, G, note in families():
        assert nx.check_planarity(nx.Graph(G.edges))[0], "%s is not planar" % name
        for k in range(2, args.kmax + 1):
            sz = reduced_cell_count(G, k)
            if sz > args.cap:
                continue
            key = "%s|%d" % (name, k)
            if key in res:
                H = {int(d): (v[0], v[1]) for d, v in res[key]["H"].items()}
                line = "(cached)"
            elif sz <= FAST_LIMIT:
                t0 = time.time()
                H, _, _ = fast.compute(G, k)
                line = "%.0fs" % (time.time() - t0)
                res[key] = {"H": paths.record(H), "N": G.N, "M": G.M, "note": note}
                paths.save_json(DATA_FILE, res)
            else:
                # too big for the full homology: get the torsion degree by degree
                t0 = time.time()
                H = {}
                for d in range(1, min(k, G.N)):
                    H[d] = (0, lean.torsion_in_degree(G, k, d))
                line = "%.0fs (torsion only)" % (time.time() - t0)
                res[key] = {"H": paths.record(H), "N": G.N, "M": G.M, "note": note,
                            "torsion_only": True}
                paths.save_json(DATA_FILE, res)
            n_tested += 1
            for d in H:
                for x in H[d][1]:
                    exponents.add(x)
                    if odd_part(x) > 1 or x not in (0, 1, 2):
                        offenders.append((name, k, d, x))
            shown = ("torsion only" if res[key].get("torsion_only")
                     else homology_line(H))
            print("%-18s %3d %3d %3d %11s %-46s %s   %s"
                  % (name, G.N, G.M, k, format(sz, ","), shown,
                     torsion_line(H), line), flush=True)
    print("-" * len(hdr))

    print("\n%d (graph, k) pairs tested, all planar." % n_tested)
    print("invariant factors seen anywhere: %s" % (sorted(exponents) or "none"))
    if offenders:
        print("\nNON-2-PRIMARY TORSION FOUND in a planar graph:")
        for o in offenders:
            print("   %s, k=%d, H_%d: invariant factor %d" % o)
        print("The planar conjecture is FALSE.")
        return 1
    print("\nNo odd torsion, and no Z/4 either: every torsion summand is Z/2.")
    print("The planar conjecture survives this search.")
    print("\nwrote data/%s" % DATA_FILE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
