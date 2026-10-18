#!/usr/bin/env python3
"""
Edge-minimal carriers of H_2 torsion at k = 4: the octahedron is not the only one.

THE QUESTION
------------
The octahedron is the carrier of H_2 torsion at k = 4 with the fewest vertices (six; nothing on
five vertices carries any torsion above H_1).  Is it the only MINIMAL carrier, in the sense that
every graph with H_2 torsion at k = 4 contains it (up to subdivision)?

THE ANSWER
----------
No.  A graph is called edge-minimal here if it carries torsion in H_2(B_4) and deleting any
single edge (keeping the graph connected) leaves H_2(B_4) torsion-free.  This script tests
edge-minimality and octahedron-containment for

    * the cube Q_3 (8 vertices, 12 edges);
    * every 7-vertex carrier with 12 edges (three graphs, one of them planar);
    * the 7-vertex, 13-edge graph  K_5 + a vertex joined to two vertices of the K_5 + a pendant
      edge  -- a carrier with a vertex of DEGREE ONE whose leaf cannot be removed: deleting the
      pendant edge kills the torsion.

All of them are edge-minimal and none contains the octahedron as a subgraph, and none has a
vertex of degree 2, so none is a subdivision of another.  The last example also shows that the
restriction to minimum degree >= 3 in the small-graph sweeps is not harmless: hanging edges can
create H_2 torsion.

The 7-vertex graphs are read from data/small_graphs_n7_k4_d1.json, which is produced by

    python3 experiments/sweep_small_graphs.py --vertices 7 --k 4 --mindeg 1

(all 853 connected 7-vertex graphs, a few minutes); this script recomputes everything it reports.

    python3 experiments/minimal_carriers.py

Requires networkx.  About a minute.
"""

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import cube, fast, from_networkx, paths, primary   # noqa: E402
from gbt.report import has_torsion                            # noqa: E402

K = 4
LEAF_EXAMPLE = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4),
                (3, 4), (3, 5), (4, 5), (5, 6)]


def h2_torsion(g):
    import networkx as nx
    H, _, _ = fast.compute(from_networkx(nx.convert_node_labels_to_integers(g)), K)
    return H.get(2, (0, []))[1], has_torsion(H, min_degree=2)


def edge_minimal(g):
    """(True if every single-edge deletion kills the H_2 torsion, number of deletions tested).

    Deleting a pendant edge also deletes the leaf it leaves isolated; deletions that
    disconnect the graph otherwise are skipped.
    """
    import networkx as nx
    tested = 0
    for e in list(g.edges()):
        h = g.copy()
        h.remove_edge(*e)
        h.remove_nodes_from([v for v in list(h) if h.degree(v) == 0])   # a deleted leaf
        if not nx.is_connected(h):
            continue
        tested += 1
        if h2_torsion(h)[1]:
            return False, tested
    return True, tested


def main():
    import networkx as nx
    from networkx.algorithms import isomorphism
    octa = nx.octahedral_graph()
    cands = [("octahedron", octa), ("cube Q_3", nx.Graph(cube().edges))]
    src = os.path.join(paths.DATA, "small_graphs_n7_k4_d1.json")
    if os.path.exists(src):
        with open(src) as fh:
            rows = json.load(fh)
        for r in rows:
            if r["torsion_above_H1"] and r["edges"] == 12:
                g = nx.Graph([tuple(e) for e in r["edgelist"]])
                cands.append(("7 vertices, 12 edges%s" % (" (planar)" if r["planar"] else ""),
                              g))
    else:
        print("(data/small_graphs_n7_k4_d1.json not found: skipping the 7-vertex 12-edge graphs)")
    cands.append(("K_5 + vertex + pendant edge", nx.Graph(LEAF_EXAMPLE)))

    hdr = "%-34s %3s %3s %-14s %-10s %-12s %s" % ("graph", "N", "M", "degrees", "H_2(B_4)",
                                                  "edge-minimal", "contains octahedron")
    print(hdr)
    print("-" * len(hdr))
    for name, g in cands:
        t, has = h2_torsion(g)
        mini, tested = edge_minimal(g) if has else (False, 0)
        contains = isomorphism.GraphMatcher(g, octa).subgraph_is_monomorphic()
        degs = "".join(str(d) for d in sorted((d for _, d in g.degree()), reverse=True))
        print("%-34s %3d %3d %-14s %-10s %-12s %s"
              % (name, g.number_of_nodes(), g.number_of_edges(), degs,
                 primary(t) if t else "free", "yes (%d)" % tested if mini else "no",
                 "yes" if contains else "no"), flush=True)
    print("-" * len(hdr))
    print("\nEvery listed graph other than the octahedron is an edge-minimal carrier that does")
    print("not contain the octahedron, so the octahedron is the carrier with the fewest vertices")
    print("but not the unique minimal one.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
