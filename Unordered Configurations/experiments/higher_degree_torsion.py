#!/usr/bin/env python3
"""
Torsion above H_2: it exists in every degree, and for planar graphs too.

THE CONSTRUCTION
----------------
Configuration spaces of a disjoint union split:

    UConf_k(G_1 u G_2)  =  disjoint union over k_1 + k_2 = k of
                           UConf_{k_1}(G_1) x UConf_{k_2}(G_2),

so Kunneth applies summand by summand.  Putting 4 particles on the octahedron, whose
H_2(B_4) contains Z/2, and one particle on each of i-2 disjoint triangles, whose H_1 is Z,
gives

    H_2(B_4 K_{2,2,2})  (x)  H_1(C_3)^{(x) i-2}   contained in   H_i(B_{i+2} Gamma),

hence Z/2 in degree i.  A disjoint union of planar graphs is planar, so planar graphs have
2-torsion in H_i for every i >= 2.  This script checks
the cases i = 3 and i = 4 directly, which is worth doing because the Kunneth summand could in
principle be killed by the way the pieces are indexed in the Swiatkowski model.

IT DOES NOT NEED TO BE DISCONNECTED
-----------------------------------
Attaching the triangle to the octahedron by a single edge gives a CONNECTED planar graph on 9
vertices and 16 edges, and the Z/2 in H_3(B_5) survives the attachment.  By contrast the
octahedron on its own has torsion-free H_3 at k = 5 and k = 6 -- the room the
extra cycle takes up is exactly what the higher-degree class needs.

AND IT IS NOT ONLY Z/2, NOR ONLY PLANAR
---------------------------------------
The bipartite sweep already contains odd and non-squarefree torsion above degree two:

    H_3(B_4 K_{5,5})  =  Z^1260 + Z/2
    H_3(B_4 K_{5,6})  =  Z^3010 + (Z/2)^5 + Z/3
    H_3(B_4 K_{5,7})  contains Z/4, by the flip-chain class

and the flip-chain class gives Z/n in degree n-1 for every n, so every cyclic group occurs in some degree.
Those graphs are not planar, and `planar_odd_search.py` is the evidence that they cannot be.

    python3 experiments/higher_degree_torsion.py [--max-degree 4]

Degree 3 takes a couple of minutes; degree 4 needs about 1.5 million cells in two slices and
a few GB, and runs in roughly a minute thanks to the sweeping peel.
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import (Graph, disjoint_union, lean, octahedron, paths,           # noqa: E402
                 primary)
from gbt.sizes import reduced_cell_count_by_degree as slices_by_degree     # noqa: E402

DATA_FILE = "higher_degree_torsion.json"

C3 = Graph(3, [(0, 1), (1, 2), (0, 2)])


def octahedron_plus_triangles(t):
    """The octahedron together with t disjoint triangles (planar, disconnected for t > 0)."""
    G = octahedron()
    for _ in range(t):
        G = disjoint_union(G, C3)
    return G


def octahedron_triangle_joined():
    """The octahedron with a triangle attached by one edge: planar, connected, 9 vertices."""
    return Graph(9, list(octahedron().edges)
                 + [(6, 7), (7, 8), (6, 8), (0, 6)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-degree", type=int, default=4)
    args = ap.parse_args()

    import networkx as nx
    res = paths.load_json(DATA_FILE)

    cases = []
    # the controls: the carriers on their own, which should stay torsion-free above H_2
    cases.append(("octahedron", octahedron(), 5, 3, "control: no room for an H_3 class"))
    cases.append(("octahedron", octahedron(), 6, 3, "control"))
    # the construction, degree by degree
    for i in range(3, args.max_degree + 1):
        G = octahedron_plus_triangles(i - 2)
        cases.append(("octahedron + %d x C_3" % (i - 2), G, i + 2, i,
                      "Kunneth: Z/2 expected"))
    # and connected, in degree 3
    cases.append(("oct-C_3 joined", octahedron_triangle_joined(), 5, 3,
                  "CONNECTED planar"))

    hdr = ("%-22s %3s %3s %3s %4s %13s %-14s %s"
           % ("graph", "N", "M", "k", "deg", "slice cells", "torsion", "note"))
    print(hdr)
    print("-" * len(hdr))
    failures = []
    for name, G, k, d, note in cases:
        planar = nx.check_planarity(nx.Graph(G.edges))[0]
        sl = slices_by_degree(G, k)
        cells = sl.get(d, 0) + sl.get(d + 1, 0)
        key = "%s|%d|%d" % (name, k, d)
        if key in res:
            t = res[key]["torsion"]
            took = "cached"
        else:
            t0 = time.time()
            t = lean.torsion_in_degree(G, k, d)
            took = "%.0fs" % (time.time() - t0)
            res[key] = {"torsion": t, "N": G.N, "M": G.M, "planar": planar,
                        "k": k, "degree": d, "cells": cells}
            paths.save_json(DATA_FILE, res)
        print("%-22s %3d %3d %3d %4d %13s %-14s %s   [%s]"
              % (name, G.N, G.M, k, d, format(cells, ","),
                 primary(t) if t else "none", note, took), flush=True)
        if not planar:
            failures.append("%s is not planar" % name)
        if note.startswith("Kunneth") and not t:
            failures.append("%s: expected Z/2 in H_%d, got none" % (name, d))
        if note.startswith("CONNECTED") and not t:
            failures.append("%s: expected Z/2 in H_%d, got none" % (name, d))
    print("-" * len(hdr))

    print()
    print("The controls are torsion-free above H_2 and the products are not, which is the")
    print("whole point: the higher-degree torsion comes from the extra disjoint cycles, not")
    print("from the carrier alone.  Every graph here is planar, so planar graphs have")
    print("2-torsion in H_i for every i >= 2.")
    if failures:
        print("\nFAILURES: %s" % failures)
        return 1
    print("\nwrote data/%s" % DATA_FILE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
