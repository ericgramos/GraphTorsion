#!/usr/bin/env python3
"""
The other half of the Ko--Park question: linklessly embeddable graphs WITH H_2 torsion.

`petersen_family.py` shows that being intrinsically linked does not force torsion above H_1.
This script shows the converse fails too, which is the more interesting direction: graphs
that are linklessly embeddable -- several of them planar, and so as far from the Petersen
family as possible -- nevertheless carry Z/2 in H_2.

The controls:

    K_5, K_{3,3}        the Kuratowski graphs: non-planar, but linklessly embeddable.
                        They carry the Ko--Park Z/2 in H_1 and nothing above it.
    K_{3,3}+e           a slightly denser non-planar graph, same conclusion.
    apex(planar)        a planar graph plus one vertex joined to everything.  Apex graphs
                        are exactly the standard family of linklessly embeddable graphs
                        that are not planar, so they are the sharpest controls available.
    maximal planar      a planar triangulation, grown by stacking a vertex into a random
                        face, hence linklessly embeddable with as many edges as planarity
                        allows.
    octahedron          K_{2,2,2}: planar, 6 vertices, 12 edges -- and Z/2 in H_2 at k=4.
    K_6 - perfect matching
                        K_6 with three independent edges removed.  This graph IS the
                        octahedron K_{2,2,2}, so this row is a consistency check on the
                        graph construction (it must agree with the previous row), not a
                        new control.

The planar triangulations are generated from a fixed seed, so the run is reproducible.

    python3 experiments/linkless_controls.py [--k 3,4]
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
    """A planar triangulation on n vertices, grown by stacking into a random face."""
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
