#!/usr/bin/env python3
"""
Do the Platonic solids have torsion in H_*(B_k)?  Four of the five do, and all at k = 4.

THE RESULT
----------
    solid          N   M   degree   k<=3        k=4
    tetrahedron    4   6   3        free        free  (and free for every k <= 10)
    cube           8  12   3        free        Z/2 in H_2
    octahedron     6  12   4        free        Z/2 in H_2
    dodecahedron  20  30   3        free        Z/2 in H_2
    icosahedron   12  30   5        free        Z/2 in H_2

So the tetrahedron is the exception, and it is not a matter of waiting for more particles:
H_*(B_k K_4) is torsion-free for every k up to 10, with rank H_2 = 6k - 15 growing linearly.
Nor is it only the Platonic solids: the triangular prism, which is the 3-prism just as the
cube is the 4-prism, is torsion-free for every k up to 7 (k <= 4 is stored in data/; k = 5, 6, 7
were checked separately with fast.compute and are not stored).

WHY k = 4
---------
An--Knudsen (J. Topol. 15 (2022) 666-691) prove that for a PLANAR graph, H_2 of the
configuration spaces is generated -- under embedding, disjoint union and edge stabilization --
by three atoms: the cycle, the star and the theta graph.  Of those three, only the theta has
nonzero H_2 at all:

    cycle  UConf_k(S^1) is a circle, so H_2 = 0 for every k;
    star   configuration spaces of trees have free fundamental group (Farley--Sabalka), so
           H_i = 0 for i >= 2 and every k;
    theta  H_2(B_k theta) = 0 for k <= 3, and Z for k = 4.

The threshold in the last line is the answer to "why four particles": below k = 4 a planar
graph has no second homology to speak of beyond products of circles, and at k = 4 the theta
class switches on.  The script verifies all three rows.

WHAT THE TORSION IS
-------------------
It cannot be a new generator, since a direct sum of torsion-free groups has no torsion.  It is
the failure of a free theta class to be integral in the ambient graph: there is a cycle z,
pushed forward from an embedded theta, with 2z a boundary and z not.  `theta_class.py`
exhibits this for the octahedron, where one embedded theta accounts for the entire Z/2.

NO CLEAN CRITERION
------------------
The mechanism is common to all four; a combinatorial criterion is not, and the data says
fairly firmly that none exists.  All five cubic graphs on eight vertices have twelve edges,
and exactly two of them have H_2 torsion at k = 4 -- the cube, and one non-planar graph.
Among the eighteen graphs on seven vertices with twelve edges and minimum degree three, three
have the torsion and fifteen do not, and the three are not separated from the fifteen by
degree sequence, girth, planarity, connectivity or triangle count.  What IS visible is a
threshold: no graph of minimum degree three with at most eleven edges has H_2 torsion at
k = 4, in any case computed.  The property is not minor-closed (density_interpolation.py), so no forbidden-minor
characterisation can exist.

    python3 experiments/platonic_solids.py [--kmax 3] [--include-k4] [--extras]

    --kmax N        compute the full homology for k = 2..N (cheap for N <= 3)
    --include-k4    also do k = 4, where the torsion is.  The dodecahedron and the
                    icosahedron need `gbt.lean.torsion_in_degree`, which uses the
                    "torsion is a cokernel" lemma to avoid the lower differentials
                    entirely; expect a few minutes for the dodecahedron and considerably
                    longer, with several GB, for the icosahedron.
    --extras        the prism and antiprism families, the three An--Knudsen atoms, and the
                    tetrahedron pushed to high k
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import (Graph, antiprism, fast, lean, paths, platonic_solids,     # noqa: E402
                 prism_graph, primary, reduced_cell_count, tetrahedron, theta)
from gbt.report import homology_line, torsion_line                        # noqa: E402

DATA_FILE = "platonic_solids.json"
FAST_LIMIT = 300_000         # above this, use the lean torsion-only path


def full_homology(G, k):
    """{degree: (free, invariant factors)} by the ordinary pipeline."""
    H, _, _ = fast.compute(G, k)
    return H


def torsion_only(G, k, d=2, verbose=False):
    """The torsion of H_d alone, via the identity "torsion is a cokernel" (cheap for big cases)."""
    return lean.torsion_in_degree(G, k, d, verbose=verbose)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kmax", type=int, default=3)
    ap.add_argument("--include-k4", action="store_true")
    ap.add_argument("--extras", action="store_true")
    args = ap.parse_args()

    res = paths.load_json(DATA_FILE)

    # ------------------------------------------------------------------ k <= kmax
    print("== the five Platonic solids, full homology for k <= %d ==\n" % args.kmax)
    hdr = "%-14s %3s %3s %4s %10s %-46s %s" % ("solid", "N", "M", "k", "cells",
                                               "homology", "torsion")
    print(hdr)
    print("-" * len(hdr))
    for name, G in platonic_solids():
        for k in range(2, args.kmax + 1):
            sz = reduced_cell_count(G, k)
            if sz > FAST_LIMIT:
                print("%-14s %3d %3d %4d %10s  (skipped: too large for the full pipeline)"
                      % (name, G.N, G.M, k, format(sz, ",")))
                continue
            t0 = time.time()
            H = full_homology(G, k)
            print("%-14s %3d %3d %4d %10s %-46s %s   [%.0fs]"
                  % (name, G.N, G.M, k, format(sz, ","), homology_line(H),
                     torsion_line(H), time.time() - t0), flush=True)
            res["%s|%d" % (name, k)] = {"H": paths.record(H), "N": G.N, "M": G.M}
            paths.save_json(DATA_FILE, res)
    print("-" * len(hdr))
    print("No Platonic solid has any torsion at k <= 3.  H_1 is torsion-free throughout,")
    print("as Ko--Park requires: all five are planar.")

    # ------------------------------------------------------------------ k = 4
    if args.include_k4:
        print("\n== k = 4, where the torsion is ==\n")
        hdr = "%-14s %3s %3s %12s %-26s %s" % ("solid", "N", "M", "cells",
                                               "torsion of H_2", "route")
        print(hdr)
        print("-" * len(hdr))
        for name, G in platonic_solids():
            key = "%s|4" % name
            sz = reduced_cell_count(G, 4)
            if key in res and "H2_torsion" in res[key]:
                t = res[key]["H2_torsion"]
                route = "cached"
            elif sz <= FAST_LIMIT:
                H = full_homology(G, 4)
                t = list(H.get(2, (0, []))[1])
                route = "full homology"
                res[key] = {"H": paths.record(H), "N": G.N, "M": G.M, "H2_torsion": t}
            else:
                t0 = time.time()
                t = torsion_only(G, 4, 2, verbose=True)
                route = "lean, torsion only, %.0fs" % (time.time() - t0)
                res[key] = {"N": G.N, "M": G.M, "H2_torsion": t}
            paths.save_json(DATA_FILE, res)
            print("%-14s %3d %3d %12s %-26s %s"
                  % (name, G.N, G.M, format(sz, ","), primary(t) if t else "none", route),
                  flush=True)
        print("-" * len(hdr))

    # ------------------------------------------------------------------ extras
    if args.extras:
        print("\n== why k = 4: the three An--Knudsen planar atoms ==\n")
        cycle = Graph(4, [(0, 1), (1, 2), (2, 3), (0, 3)])
        star = Graph(5, [(0, 1), (0, 2), (0, 3), (0, 4)])
        print("%-16s %4s %-44s %s" % ("atom", "k", "homology", "H_2"))
        for nm, G in [("cycle C_4", cycle), ("star K_{1,4}", star),
                      ("theta (3 paths)", theta(3, 2)), ("theta (4 paths)", theta(4, 2))]:
            for k in range(2, 6):
                H = full_homology(G, k)
                print("%-16s %4d %-44s %s"
                      % (nm, k, homology_line(H),
                         "Z^%d" % H[2][0] if H.get(2, (0, []))[0] else "0"), flush=True)
            print()
        print("Only the theta contributes to H_2 at all, and the three-path theta does so")
        print("first at k = 4 -- exactly where the planar torsion appears.")

        print("\n== the tetrahedron at high k: does it ever acquire torsion? ==\n")
        print("%4s %10s %-54s %s" % ("k", "cells", "homology", "torsion"))
        for k in range(2, 11):
            G = tetrahedron()
            H = full_homology(G, k)
            print("%4d %10s %-54s %s"
                  % (k, format(reduced_cell_count(G, k), ","), homology_line(H),
                     torsion_line(H)), flush=True)
        print("rank H_2(B_k K_4) = 6k - 15 for k >= 3: linear, and torsion-free throughout.")

        print("\n== the prism and antiprism families at k = 4 ==")
        print("   (the cube is the 4-prism and the octahedron is the 3-antiprism)\n")
        print("%-16s %3s %3s %10s %-20s %s" % ("graph", "N", "M", "cells",
                                               "torsion of H_2", ""))
        for n in range(3, 7):
            for fam, make in (("prism", prism_graph), ("antiprism", antiprism)):
                G = make(n)
                sz = reduced_cell_count(G, 4)
                if sz > FAST_LIMIT:
                    print("%-16s %3d %3d %10s  (skipped: too large)"
                          % ("%s_%d" % (fam, n), G.N, G.M, format(sz, ",")))
                    continue
                H = full_homology(G, 4)
                t = H.get(2, (0, []))[1]
                note = ""
                if fam == "prism" and n == 4:
                    note = "= the cube"
                if fam == "antiprism" and n == 3:
                    note = "= the octahedron"
                print("%-16s %3d %3d %10s %-20s %s"
                      % ("%s_%d" % (fam, n), G.N, G.M, format(sz, ","),
                         primary(t) if t else "none", note), flush=True)
                res["%s_%d|4" % (fam, n)] = {"H": paths.record(H), "N": G.N, "M": G.M}
                paths.save_json(DATA_FILE, res)

    print("\nwrote data/%s" % DATA_FILE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
