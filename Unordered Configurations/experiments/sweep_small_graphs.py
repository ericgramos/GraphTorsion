#!/usr/bin/env python3
"""
Exhaustive sweeps over all small graphs: where does H_2 torsion first appear?

It enumerates every
connected graph on a given number of vertices with a given minimum degree -- using
networkx's graph atlas, which is complete up to seven vertices -- and computes
H_*(B_k) for each.  The questions it answers:

    * Does any graph on at most five vertices have torsion in H_i for i >= 2?  (No.  This
      is what rules out Kunneth products with an H_0 factor as a source of H_2 torsion, so
      products need at least 8 vertices: see decomposability_K8.py.)
    * Which is the smallest graph with H_2 torsion at k = 3?  (K_7, on 7 vertices and 21
      edges.)
    * At k = 4?  (The octahedron K_{2,2,2}, on 6 vertices and 12 edges -- which is PLANAR.
      That is the observation that kills any hope of a planarity criterion in degree 2.)
    * Is the property of having H_2 torsion closed under adding edges?  (No: see
      density_interpolation.py.  Hence it is not minor-closed, and there is no
      Kuratowski-type theorem above degree one.)

The minimum-degree restriction makes the sweep much cheaper.  A vertex of degree 2 can be
smoothed away without changing the configuration spaces (they are homeomorphism invariants),
so --mindeg 2 loses nothing but subdivisions.  A vertex of degree 1 CANNOT in general be
deleted -- a hanging edge changes H_*(B_k) -- so --mindeg 3 excludes the graphs with leaves;
run with --mindeg 1 to include them.

    python3 experiments/sweep_small_graphs.py [--vertices N] [--k K] [--mindeg D]

The default (--vertices 6 --k 3,4 --mindeg 3) is one of the stored sweeps; all of them, in increasing cost:

    --vertices 5 --k 2,3,4,5 --mindeg 2      the "no small carrier" check (seconds)
    --vertices 6 --k 4       --mindeg 3      19 graphs; finds the octahedron (seconds)
    --vertices 7 --k 4       --mindeg 3      150 graphs (minutes); stored in data/
    --vertices 7 --k 3       --mindeg 3      150 graphs; K_7 is the only carrier (seconds)
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import from_networkx, paths                    # noqa: E402
from gbt import fast                                   # noqa: E402
from gbt.report import has_torsion, homology_line, torsion_line   # noqa: E402


def atlas_graphs(nvert, mindeg):
    """All connected graphs on `nvert` vertices with minimum degree >= `mindeg`.

    networkx's atlas is complete for up to 7 vertices; beyond that this raises, rather
    than silently sweeping an incomplete list.
    """
    import networkx as nx
    from networkx.generators.atlas import graph_atlas_g
    if nvert > 7:
        raise ValueError("the graph atlas only goes to 7 vertices; "
                         "for more, generate the graphs with nauty/geng instead")
    out = []
    for Gx in graph_atlas_g():
        if Gx.number_of_nodes() != nvert or not nx.is_connected(Gx):
            continue
        if min(d for _, d in Gx.degree()) < mindeg:
            continue
        out.append(Gx)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vertices", type=int, default=6)
    ap.add_argument("--k", default="3,4", help="comma-separated particle numbers")
    ap.add_argument("--mindeg", type=int, default=3)
    args = ap.parse_args()
    ks = [int(x) for x in args.k.split(",")]

    import networkx as nx
    graphs = atlas_graphs(args.vertices, args.mindeg)
    print("%d connected graphs on %d vertices with min degree >= %d, at k in %s\n"
          % (len(graphs), args.vertices, args.mindeg, ks), flush=True)

    hdr = ("%4s %7s %7s %4s %-34s %s"
           % ("#", "edges", "planar", "k", "homology", "torsion"))
    print(hdr)
    print("-" * len(hdr))

    rows = []
    for i, Gx in enumerate(graphs):
        G = from_networkx(Gx)
        planar = nx.check_planarity(Gx)[0]
        degseq = sorted((d for _, d in Gx.degree()), reverse=True)
        for k in ks:
            t0 = time.time()
            H, _, _ = fast.compute(G, k)
            tl = torsion_line(H)
            interesting = has_torsion(H, min_degree=2)
            rows.append({"edges": Gx.number_of_edges(), "planar": planar, "k": k,
                         "degrees": degseq, "torsion": tl,
                         "torsion_above_H1": interesting,
                         "edgelist": sorted(tuple(e) for e in Gx.edges())})
            print("%4d %7d %7s %4d %-34s %s%s"
                  % (i, Gx.number_of_edges(), planar, k, homology_line(H), tl,
                     "   <== torsion above H_1" if interesting else ""), flush=True)
            del t0

    print("-" * len(hdr))

    hits = [r for r in rows if r["torsion_above_H1"]]
    print("\n%d of %d (graph, k) pairs have torsion in H_i for some i >= 2"
          % (len(hits), len(rows)))
    if hits:
        least = min(r["edges"] for r in hits)
        print("fewest edges among them: %d" % least)
        for r in sorted(hits, key=lambda r: (r["edges"], r["k"])):
            print("   m=%2d k=%d planar=%-5s degrees=%s  torsion %s"
                  % (r["edges"], r["k"], r["planar"], r["degrees"], r["torsion"]))
        if any(r["planar"] for r in hits):
            print("\nNOTE: at least one PLANAR graph carries torsion above H_1.")
            print("Planarity is therefore not an obstruction in degree 2;")
            print("see theta_class.py for where that class comes from.")
    else:
        print("no graph of this size carries torsion above H_1 at these k.")

    name = "small_graphs_n%d_k%s_d%d.json" % (args.vertices, args.k.replace(",", "-"),
                                              args.mindeg)
    paths.save_json(name, rows)
    print("\nwrote data/%s" % name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
