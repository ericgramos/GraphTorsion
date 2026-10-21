#!/usr/bin/env python3
"""
Planar graphs with torsion in H_2, and the three planar generators of An--Knudsen.

THE QUESTION
------------
In degree one, Ko--Park says torsion happens exactly for non-planar graphs.  In degree two that
is false: the octahedron K_{2,2,2} is planar and carries Z/2 in H_2(B_4).  This script
establishes that, and then puts it next to the theorem that makes it surprising.

An and Knudsen (J. Topol. 15 (2022) 666-691) prove that for a PLANAR graph, H_2 of the
configuration spaces is generated, under embedding, disjoint union and edge stabilization, by
three atoms: the cycle, the star, and the theta graph.  This script computes the homology of
all three and finds it torsion-free.  Since torsion cannot appear in a direct sum of
torsion-free groups, the 2-torsion in the octahedron cannot be a "new generator" -- it has to
be produced by the RELATIONS, i.e. a free class on an embedded atom must map to a torsion
class.  `theta_class.py` verifies exactly that.

The script also computes the graphs around the octahedron -- the triangular prism (same vertex
count, no torsion), the wheel W_5, K_4, K_{2,3}, K_{2,4} -- so the finding is not a one-off.

METHOD
------
The nine graphs of cases() are computed in full (H_*(B_k Gamma; Z) in every degree) with
`gbt.fast`, or with `gbt.lean` above LEAN_THRESHOLD = 600,000 reduced cells.  A graph "carries
torsion above H_1" if some H_i with i >= 2 has a non-trivial torsion subgroup.

WHAT IS PRINTED
---------------
1.  For every (graph, k) NOT already stored: one progress line
        name  n=<vertices> m=<edges> k=<k>  <raw cells> cells  [time]  H_1=... H_2=...
    followed by the torsion.  Stored cases are not recomputed and print nothing here.
2.  A table, one row per (graph, k), columns
        graph | n (vertices) | m (edges) | k | torsion (all degrees) |
        above H_1 (YES / no) | note
3.  The verdict: whether the three atoms are torsion-free, and the list of planar carriers of
    torsion above H_1, with a short explanation.

Excerpt (default `--k 3,4`, everything stored, about a second):

    graph            n   m    k torsion                        above H_1  note
    cycle C_4        4   4    3 none                           no         An--Knudsen atom
    ...
    octahedron       6  12    4 H_2: Z/2                       YES        planar, 6 vertices ...
    An--Knudsen's three planar atoms are torsion-free in H_2: True
    Planar graphs carrying torsion above H_1: [('octahedron', 4)]

USAGE
-----
    python3 experiments/planar_torsion.py [--k 3,4]

    --k LIST    comma-separated particle numbers (default "3,4").  At k = 3 nothing planar has
                torsion; the octahedron needs k = 4.

DATA FILE
---------
data/planar_torsion.json, keyed "<name>|<k>" (e.g. "octahedron|4"), each record
{"H": {degree: [free rank, [invariant factors > 1]]}, "vertices": N, "edges": M, "note": text}.

COST AND EXIT CODE
------------------
About a second from scratch at k <= 4.  Always exits 0.  To recompute a stored case delete its
key from the data file.
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import (Graph, complete, complete_bip, octahedron, prism,          # noqa: E402
                 theta, wheel, reduced_cell_count, paths)
from gbt import fast, lean                                                 # noqa: E402
from gbt.report import has_torsion, homology_line, torsion_line             # noqa: E402

#: above this many reduced cells, use the memory-lean pipeline instead of gbt.fast
LEAN_THRESHOLD = 600_000
DATA_FILE = "planar_torsion.json"


def cases():
    """The nine test graphs, as a list of (name, Graph, note).

    The first three are An--Knudsen's planar atoms (their note contains the word "atom", which
    main() uses to recognise them): the 4-cycle, the star K_{1,4} and the theta graph with
    three paths of length 2.  Then K_4, the triangular prism, the wheel W_5, K_{2,3}, K_{2,4}
    and the octahedron.  All are planar."""
    star = Graph(5, [(0, 1), (0, 2), (0, 3), (0, 4)])
    cyc = Graph(4, [(0, 1), (1, 2), (2, 3), (0, 3)])
    return [
        ("cycle C_4", cyc, "An--Knudsen atom"),
        ("star K_{1,4}", star, "An--Knudsen atom"),
        ("theta", theta(3, 2), "An--Knudsen atom"),
        ("K_4", complete(4), "planar, smallest non-trivial"),
        ("prism", prism(), "planar, 6 vertices 9 edges -- control"),
        ("wheel W_5", wheel(5), "planar"),
        ("K_{2,3}", complete_bip(2, 3), "planar"),
        ("K_{2,4}", complete_bip(2, 4), "planar"),
        ("octahedron", octahedron(), "planar, 6 vertices 12 edges"),
    ]


def main():
    """Compute or reload every (graph, k) of cases(), print the table and the verdict.

    Returns 0 (the process exit code)."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", default="3,4", help="comma-separated particle numbers")
    args = ap.parse_args()
    ks = [int(x) for x in args.k.split(",")]

    res = paths.load_json(DATA_FILE)
    rows = []
    for name, G, note in cases():
        for k in ks:
            key = "%s|%d" % (name, k)
            if key not in res:
                sz = reduced_cell_count(G, k)
                mod = lean if sz > LEAN_THRESHOLD else fast
                t0 = time.time()
                H, raw, _ = mod.compute(G, k)
                print("%-14s n=%-2d m=%-2d k=%d  %8s cells  [%.0fs]  %s"
                      % (name, G.N, G.M, k, format(sum(raw.values()), ","),
                         time.time() - t0, homology_line(H)))
                print("%-14s %s torsion: %s" % ("", " " * 18, torsion_line(H)), flush=True)
                res[key] = {"H": paths.record(H), "vertices": G.N, "edges": G.M,
                            "note": note}
                paths.save_json(DATA_FILE, res)
            # metadata from the live graph, so records cached by an earlier version
            # of this script (which may hold only the homology) still work
            rows.append((name, k, G.N, G.M,
                         {int(d): tuple(v) for d, v in res[key]["H"].items()}, note))

    print("\n" + "=" * 90)
    print("planar graphs: torsion in H_i for i >= 2?")
    print("=" * 90)
    print("%-14s %3s %3s %4s %-30s %-10s %s"
          % ("graph", "n", "m", "k", "torsion", "above H_1", "note"))
    carriers = []
    atoms_ok = True
    for name, k, n, m, H, note in rows:
        above = has_torsion(H, min_degree=2)       # torsion in some H_i with i >= 2
        if above:
            carriers.append((name, k))
            if "atom" in note:
                atoms_ok = False
        print("%-14s %3d %3d %4d %-30s %-10s %s"
              % (name, n, m, k, torsion_line(H), "YES" if above else "no", note))

    print()
    print("An--Knudsen's three planar atoms are torsion-free in H_2: %s" % atoms_ok)
    if carriers:
        print("Planar graphs carrying torsion above H_1: %s" % carriers)
        print()
        print("These two facts together say something definite.  H_2 of a planar graph's")
        print("configuration spaces is generated by the three atoms under embedding,")
        print("disjoint union and stabilization; all three atoms are torsion-free; and yet")
        print("the octahedron has Z/2.  So the torsion is created by relations, not carried")
        print("by a generator -- a free class on an embedded atom maps to a torsion class.")
        print("theta_class.py exhibits that map explicitly.")
    else:
        print("No planar carrier found at these k.  The octahedron needs k = 4.")
    print("\nwrote data/%s" % DATA_FILE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
