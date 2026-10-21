#!/usr/bin/env python3
"""
The other half of the Ko--Park question: linklessly embeddable graphs WITH H_2 torsion.

THE QUESTION
------------
Ko--Park: H_1(B_k Gamma) has torsion (for k >= 2) exactly when Gamma is not planar.  Is there
an analogue one degree up, with planarity replaced by linkless embeddability (whose obstruction
set is the Petersen family)?  `petersen_family.py` shows that being intrinsically linked does
not force torsion above H_1 at k = 3.  This script tests the converse, which is the more
interesting direction: do linklessly embeddable graphs stay free of torsion above H_1?  They
do not.  The apex graph apex(planar8) already has (Z/2)^4 in H_2 at k = 3, and the planar
octahedron has Z/2 in H_2 at k = 4.

THE CONTROLS
------------
All of these are linklessly embeddable.

    K_5, K_{3,3}        the Kuratowski graphs: non-planar, but linklessly embeddable.
                        They carry the Ko--Park Z/2 in H_1 and nothing above it at k = 3.
    K_{3,3}+e           K_{3,3} plus one edge inside a side: a slightly denser non-planar
                        control.
    apex(planar7),      a maximal planar graph on 7 (resp. 8) vertices plus one new vertex
    apex(planar8)       joined to every vertex.  Apex graphs are linklessly embeddable (a
                        classical fact), and are the natural non-planar controls.
    planar9             a maximal planar graph on 9 vertices: as many edges as planarity
                        allows.
    K_{2,2,2}           the octahedron: planar, 6 vertices, 12 edges -- Z/2 in H_2 at k = 4.
    K_6 - PM            K_6 with the perfect matching {0,1}, {2,3}, {4,5} removed.  This graph
                        IS the octahedron, so this row is a consistency check on the graph
                        construction (it must agree with the previous row), not a new control.

The maximal planar graphs are grown by stacking a vertex into a random face of a
triangulation, from fixed seeds (7: seed 1, 8: seed 2, 9: seed 3), so the run is reproducible.

WHAT IS PRINTED
---------------
For each (graph, k) not yet stored, a progress block (name, n vertices, m edges, planarity,
k, raw cell count, time) with the full homology and the torsion.  Then a table:

    graph        the control's name (also the data key prefix)
    k            particles
    planar       True / False
    edges        number of edges
    torsion      all torsion, by degree, as primary decompositions
    above H_1?   "YES" if some H_i with i >= 2 has torsion

and a verdict listing the linklessly embeddable graphs with torsion above H_1.  Real output
at k = 3 (from the stored data):

    graph               k  planar  edges torsion                    above H_1?
    K_5                 3   False     10 H_1: Z/2                   no
    apex(planar8)       3   False     26 H_1: Z/2   H_2: (Z/2)^4    YES
    planar9             3    True     21 none                       no
    K_{2,2,2}           3    True     12 none                       no
    ...
    Linklessly embeddable graphs WITH torsion above H_1: [('apex(planar8)', 3, False)]

USAGE
-----
    python3 experiments/linkless_controls.py [--k 3]

    --k LIST   (default "3") comma-separated particle numbers, e.g. --k 3,4.

Results are stored in data/linkless_controls.json under "NAME|K" with record
{"H": {degree: [free rank, [invariant factors]]}, "planar", "vertices", "edges", "note"}
(older records carry only "H").  The shipped data hold every control at k = 3 and K_5 at k = 4;
a run with --k 3 is instant, and --k 3,4 computes the remaining k = 4 cases, which takes a few
minutes and up to a few GB (the apex graphs are the largest; `gbt.lean` is used above 600,000
reduced cells).  Needs networkx.

Exit code: always 0.
"""

import argparse
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import from_networkx, reduced_cell_count, paths              # noqa: E402
from gbt import fast, lean                                           # noqa: E402
from gbt.report import has_torsion, homology_line, torsion_line        # noqa: E402

LEAN_THRESHOLD = 600_000
DATA_FILE = "linkless_controls.json"


def maximal_planar(n, seed):
    """A planar triangulation on n vertices, grown by stacking into a random face.

    Start from K_4 with its four triangular faces; each new vertex v is joined to the three
    corners of a face chosen by random.Random(seed), which is replaced by the three new faces.
    Every step keeps the graph a triangulation, so the result is maximal planar (3n - 6 edges).
    Returns a networkx graph.
    """
    import networkx as nx
    G = nx.complete_graph(4)
    rng = random.Random(seed)
    faces = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    for v in range(4, n):
        f = faces.pop(rng.randrange(len(faces)))
        for x in f:
            G.add_edge(v, x)
        faces += [(v, f[0], f[1]), (v, f[0], f[2]), (v, f[1], f[2])]
    return G


def controls():
    """[(name, Graph, is planar, note)] -- all linklessly embeddable."""
    import networkx as nx
    out = []

    def add(name, Gx, note):
        """Convert Gx to a gbt Graph and record (name, Graph, planarity of Gx, note)."""
        out.append((name, from_networkx(Gx), nx.check_planarity(Gx)[0], note))

    add("K_5", nx.complete_graph(5), "Kuratowski, linklessly embeddable")
    add("K_{3,3}", nx.complete_bipartite_graph(3, 3), "Kuratowski")
    g = nx.complete_bipartite_graph(3, 3)
    g.add_edge(0, 1)
    add("K_{3,3}+e", g, "denser non-planar control")
    for n, seed in [(7, 1), (8, 2)]:
        g = maximal_planar(n, seed)
        a = max(g.nodes()) + 1
        for v in list(g.nodes()):
            g.add_edge(a, v)
        add("apex(planar%d)" % n, g, "apex: linklessly embeddable, not planar")
    add("planar9", maximal_planar(9, 3), "maximal planar")
    add("K_{2,2,2}", nx.complete_multipartite_graph(2, 2, 2), "octahedron: PLANAR")
    g = nx.complete_graph(6)
    g.remove_edges_from([(0, 1), (2, 3), (4, 5)])
    add("K_6 - PM", g, "K_6 minus a perfect matching (= the octahedron; consistency check)")
    return out


def main():
    """Compute (or load) every control at every k, print the table and the verdict."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", default="3", help="comma-separated particle numbers")
    args = ap.parse_args()
    ks = [int(x) for x in args.k.split(",")]

    res = paths.load_json(DATA_FILE)
    rows = []
    for name, G, planar, note in controls():
        for k in ks:
            key = "%s|%d" % (name, k)
            sz = reduced_cell_count(G, k)
            if key not in res:
                mod = lean if sz > LEAN_THRESHOLD else fast
                t0 = time.time()
                H, raw, red = mod.compute(G, k, verbose=(sz > LEAN_THRESHOLD))
                print("%-16s n=%-2d m=%-2d planar=%-5s k=%d  %s cells  [%.0fs]"
                      % (name, G.N, G.M, planar, k, format(sum(raw.values()), ","),
                         time.time() - t0))
                print("                 %s" % homology_line(H))
                print("                 torsion: %s" % torsion_line(H), flush=True)
                res[key] = {"H": paths.record(H), "planar": planar,
                            "vertices": G.N, "edges": G.M, "note": note}
                paths.save_json(DATA_FILE, res)
            # metadata comes from the live graph, not the stored record: results cached
            # by an earlier version of this script may carry only the homology
            H = {int(d): tuple(v) for d, v in res[key]["H"].items()}
            rows.append((name, k, planar, G.M, H, note))

    print("\n" + "=" * 86)
    print("linklessly embeddable graphs: is H_2 torsion-free?")
    print("=" * 86)
    print("%-16s %4s %7s %6s %-26s %s"
          % ("graph", "k", "planar", "edges", "torsion", "above H_1?"))
    offenders = []
    for name, k, planar, m, H, note in rows:
        above = has_torsion(H, min_degree=2)
        if above:
            offenders.append((name, k, planar))
        print("%-16s %4d %7s %6d %-26s %s"
              % (name, k, planar, m, torsion_line(H), "YES" if above else "no"))
    print()
    if offenders:
        print("Linklessly embeddable graphs WITH torsion above H_1: %s" % offenders)
        print("Combined with petersen_family.py, this rules out any Ko--Park analogue")
        print("for intrinsic linking in degree two: the implication fails both ways.")
        if any(p for _, _, p in offenders):
            print("Note that at least one of them is PLANAR.")
    else:
        print("None of these controls has torsion above H_1 at the k tried.")
        print("Try --k 4: the octahedron's class only appears with four particles.")
    print("\nwrote data/%s" % DATA_FILE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
