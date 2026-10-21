#!/usr/bin/env python3
"""
Edge-minimal carriers of H_2 torsion at k = 4: the octahedron is not the only one.

THE QUESTION
------------
The octahedron is the carrier of H_2 torsion at k = 4 with the fewest vertices (six; nothing on
five vertices carries any torsion above H_1).  Is it the only MINIMAL carrier, in the sense that
every graph with H_2 torsion at k = 4 contains it (up to subdivision)?  If it were, "contains an
octahedron" would be the degree-two analogue of Kuratowski's theorem at k = 4.

THE ANSWER
----------
No.  A graph is called edge-minimal here if it carries torsion above H_1 at k = 4 and every
single-edge deletion that keeps the graph connected leaves it torsion-free above H_1 (deleting
a pendant edge also deletes the leaf it isolates).  This script tests edge-minimality and
octahedron-containment for

    * the octahedron itself (for comparison);
    * every 7-vertex carrier with 12 edges (three graphs, one of them planar), read from the
      stored exhaustive sweep;
    * the 7-vertex, 13-edge graph  K_5 + a vertex joined to two vertices of the K_5 + a pendant
      edge  (LEAF_EXAMPLE: edges among {0..4}, then {3,5}, {4,5}, {5,6}) -- a carrier with a
      vertex of DEGREE ONE whose leaf cannot be removed: deleting the pendant edge kills the
      torsion.

All of them are edge-minimal and none except the octahedron contains the octahedron as a
subgraph.  Containment is tested as a subgraph monomorphism (networkx GraphMatcher).  None of the
graphs has a vertex of degree 2, and none can contain a subdivided octahedron either: that
needs six vertices of degree >= 4, which none of the 7-vertex graphs has.  The last example
also shows that the restriction to minimum degree >= 3 in some small-graph sweeps is not
harmless: hanging edges can create H_2 torsion.

WHAT IS PRINTED
---------------
One row per graph:

    graph                 a description
    N, M                  vertices and edges
    degrees               the degree sequence, largest first, written as a string of digits
    H_2(B_4)              the torsion of H_2 at k = 4, as a primary decomposition, or "free"
    edge-minimal          "yes (T)" with T the number of single-edge deletions tested, or "no"
    contains octahedron   "yes" / "no"

then a three-line conclusion.  Real output (13 seconds):

    graph                           N   M degrees  H_2(B_4) edge-minimal contains octahedron
    octahedron                      6  12 444444   Z/2      yes (12)     yes
    7 vertices, 12 edges            7  12 4443333  Z/2      yes (12)     no
    7 vertices, 12 edges (planar)   7  12 4443333  Z/2      yes (12)     no
    7 vertices, 12 edges            7  12 4443333  Z/2      yes (12)     no
    K_5 + vertex + pendant edge     7  13 5544431  Z/2      yes (13)     no

(columns narrowed here to fit; the script pads them more widely).

USAGE
-----
    python3 experiments/minimal_carriers.py

No options.  The 7-vertex graphs are read from data/small_graphs_n7_k4_d1.json, which is
produced by

    python3 experiments/sweep_small_graphs.py --vertices 7 --k 4 --mindeg 1

(all 853 connected 7-vertex graphs, a few minutes); if that file is missing those rows are
skipped with a note.  Everything reported is recomputed here with `gbt.fast`; nothing is
written.  Requires networkx.  About 15 seconds, well under 1 GB.  Exit code: always 0.
"""

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import fast, from_networkx, paths, primary   # noqa: E402
from gbt.report import has_torsion                            # noqa: E402

K = 4
LEAF_EXAMPLE = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4),
                (3, 4), (3, 5), (4, 5), (5, 6)]


def h2_torsion(g):
    """(invariant factors of tors H_2(B_4 g), True if some H_i with i >= 2 has torsion).

    g is a networkx graph with any node labels; it is relabelled to 0..N-1 for gbt.  Note that
    the second entry looks at every degree >= 2, not only H_2; edge-minimality uses it.
    """
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
    """Assemble the candidate graphs, test each, print the table and the conclusion."""
    import networkx as nx
    from networkx.algorithms import isomorphism
    octa = nx.octahedral_graph()
    cands = [("octahedron", octa)]
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
