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
diagonal form supplies the divisors 2, 3, ..., min(n, b-n+1).  An odd divisor therefore needs
n >= 3 and b >= n+2, i.e. at least four vertices on one side attached to an independent set of
at least five on the other.  The smallest instance is K_{4,5}, which gives Z/3 -- and K_{4,5}
contains K_{3,3}, so it is not planar; every larger instance contains K_{4,5}.  So the
mechanism is unavailable to planar graphs.  If it is the only source of odd torsion
(`prime_census.py`), the planar conjecture follows.

Products cannot help either: Tor(Z/2, Z/2) = Z/2 and Z/2 (x) Z = Z/2, so Kunneth classes built
from 2-torsion stay 2-primary.

So a planar counterexample would need a genuinely new mechanism.  The families below are
chosen as the places where the known mechanism comes closest to being available in a planar
graph -- planar graphs with one or two vertices of large degree and a large independent set --
together with dense and regular controls:

    generalised thetas K_{2,m}   m = 4..8.  Two hubs of degree m, m independent branch
                                 vertices.  This is the planar graph that most resembles
                                 K_{n+1,b} with a sparse large side: it IS K_{a,b} with a = 2.
                                 The vanishing for a <= n applies, so the top homology
                                 vanishes, but the lower degrees are where a new mechanism
                                 might hide.
    wheels W_m                   m = 5..10.  One hub of degree m, attached to an m-cycle.
    bipyramids                   m = 4..7.  Two hubs of degree m, attached to an m-cycle
                                 (bipyramid_4 is the octahedron).
    grids 3x3, 3x4               the other extreme: no high-degree vertex at all.
    triangulations               maximal planar graphs on 8, 9, 10, 11 vertices, grown by
                                 stacking into a random face with a fixed seed.
    Platonic solids              the tetrahedron, octahedron, dodecahedron and icosahedron.
    prisms, antiprisms           n = 5, 6; planar and 3- resp. 4-regular.

That is 29 entries, 28 graphs up to isomorphism (bipyramid_4 and the octahedron coincide);
each is checked to be planar (networkx) before it is used.

METHOD
------
For every graph and every k = 2, ..., --kmax whose fully reduced complex has at most --cap cells:
if the complex has at most FAST_LIMIT = 300,000 cells the full homology is computed with
`gbt.fast`; otherwise only the torsion of H_1, ..., H_{k-1} is computed, degree by degree, with
`gbt.lean.torsion_in_degree` (the free ranks are then not known and are stored as 0, with the
flag "torsion_only").  Every invariant factor found is then tested: anything other than 2 (an
odd prime factor, or 4, 8, ...) is an offender.

WHAT IS PRINTED
---------------
A table with one row per (graph, k) tested, columns

    graph | N (vertices) | M (edges) | k | cells (fully reduced complex) |
    homology ("H_1=... H_2=...", or "torsion only") | torsion | time ("(cached)" if stored)

then the number of pairs tested, the set of all invariant factors seen, and the verdict.
Excerpt (`--kmax 3`, everything cached, a few seconds):

    graph                N   M   k       cells homology                          torsion
    K_{2,4}              6   8   2         155 H_1=Z^6  H_2=0                    none   (cached)
    ...
    58 (graph, k) pairs tested, all planar.
    invariant factors seen anywhere: none
    No odd torsion, and no Z/4 either: every torsion summand is Z/2.
    The planar conjecture survives this search.

At the default --kmax 4 the invariant factors seen are exactly {2}: every torsion summand found
in a planar graph, in any degree, is literally Z/2.

USAGE
-----
    python3 experiments/planar_odd_search.py [--kmax 4] [--cap 600000]

    --kmax K    largest particle number tried (default 4); k runs from 2 to K.
    --cap N     skip a (graph, k) whose fully reduced complex has more than N cells
                (default 600,000).  Skipped pairs are not printed.

DATA FILE
---------
data/planar_odd_search.json, keyed "<name>|<k>", each record
{"H": {degree: [free rank, [invariant factors > 1]]}, "N": vertices, "M": edges, "note": text}
and, for pairs too big for the full homology, "torsion_only": true.  Stored pairs are reused.

COST AND EXIT CODE
------------------
A few minutes at the defaults from scratch; seconds when everything is stored.  Requires
networkx.  Exits with 1 if a non-2-primary summand is found in a planar graph (the data file
is still written), and 0 otherwise, so it doubles as a regression test on the conjecture.

PITFALL
-------
The cap is applied before the stored data are consulted, so lowering --cap hides stored pairs
from the table and the count; it does not delete them.
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
#: above this many reduced cells only the torsion is computed (gbt.lean.torsion_in_degree)
FAST_LIMIT = 300_000


def bipyramid(n):
    """The n-bipyramid: two apexes joined to every vertex of an n-cycle.

    Planar, n+2 vertices, 3n edges; the cycle is 0..n-1 and the apexes are n and n+1, each of
    degree n.  bipyramid(4) is the octahedron."""
    edges = ([(i, (i + 1) % n) for i in range(n)]
             + [(n, i) for i in range(n)]
             + [(n + 1, i) for i in range(n)])
    return Graph(n + 2, edges)


def triangulation(n, seed):
    """A maximal planar graph on n vertices (3n - 6 edges), grown by stacking.

    Start from K_4 with its four triangular faces; for v = 4, ..., n-1 pick a face uniformly at
    random (random.Random(seed)), join v to its three corners, and replace the face by the
    three new ones.  Deterministic for a given seed.  Returns a gbt Graph (needs networkx)."""
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
    """The 29 planar test graphs (28 up to isomorphism), as a list of (name, Graph, note).

    `name` is the key used in the data file, `note` a one-line description stored with each
    record.  Chosen as explained in the module docstring.  Needs networkx (for the grids)."""
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
    """Run the search over families() and print the table and verdict.

    Returns 1 if some planar graph has an invariant factor other than 2, else 0."""
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
                # reuse the stored answer; JSON keys are strings, so convert degrees to int
                H = {int(d): (v[0], v[1]) for d, v in res[key]["H"].items()}
                line = "(cached)"
            elif sz <= FAST_LIMIT:
                t0 = time.time()
                H, _, _ = fast.compute(G, k)
                line = "%.0fs" % (time.time() - t0)
                res[key] = {"H": paths.record(H), "N": G.N, "M": G.M, "note": note}
                paths.save_json(DATA_FILE, res)
            else:
                # too big for the full homology: get the torsion degree by degree, using
                # tors H_d = tors coker d_{d+1}.  H_d for d >= k is zero (the complex has
                # no cells above degree min(k, N)), and the top degree is torsion-free.
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
                    # an offender is any invariant factor other than 2: an odd prime
                    # factor, or a higher power of 2 such as 4
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
