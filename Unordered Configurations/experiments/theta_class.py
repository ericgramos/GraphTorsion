#!/usr/bin/env python3
"""
Where planar H_2 torsion comes from: a FREE class on an embedded theta subgraph.

THE PUZZLE
----------
An and Knudsen (J. Topol. 15 (2022) 666-691) prove that for a planar graph, H_2 of the
configuration spaces is generated -- under embedding, disjoint union and edge stabilization --
by three atoms: the cycle, the star and the theta graph.  All three are torsion-free in H_2.
In fact only the theta has nonzero H_2 at all: configuration spaces of trees have free
fundamental group and so no H_2, and UConf_k of a circle is a circle.  And
H_2(B_k theta) = 0 for k <= 3 and Z for k = 4 (for the theta with three paths) -- which is
exactly the particle number at which planar torsion starts appearing (`platonic_solids.py`).

And yet the octahedron and the dodecahedron are both planar and both have Z/2 in H_2(B_4).  A
direct sum of torsion-free groups has no torsion, so the resolution has to be that the
generation statement is about images, not summands: a FREE class on an embedded atom can map
to a TORSION class.  This script verifies that directly, and it supports the conjecture that
torsion in these spaces is never carried by a generator of the Z[E]-module of homologies; it
is always created by relations.

(This also corrects a plausible-sounding but false inference: that an edge-minimal carrier of
torsion must be "indecomposable" in the sense of not coming from a subgraph.  Minimality says
every edge is needed to make the class torsion; it does not say the class is not the image of
something on a subgraph.)

THE TEST
--------
Let d_3 be the third differential of the UNREDUCED Swiatkowski complex of Gamma at weight k;
the torsion of coker d_3 is the torsion of H_2 (the "baseline").  Find a theta subgraph -- two
hubs u, w joined by three internally vertex-disjoint paths -- build ITS unreduced complex on the
same vertex set, and compute a Z-basis of its 2-cycles Z_2 (a kernel lattice).  The unreduced
complex is functorial for inclusions of subgraphs (no basepoint edges are involved), so each
2-cycle pushes forward to a 2-cycle of Gamma by renaming edges.  Adjoin all pushed cycles to
the columns of d_3 and recompute the torsion of the cokernel.

Why this decides the question: let z represent the Z/2, so 2z is a boundary.  If the cokernel
of [d_3 | pushed cycles] is torsion-free, then z lies in (boundaries + pushed cycles), i.e.
[z] = i_*[y] for a cycle y on the theta: the torsion class is the image of a class of the
theta, whose H_2 is free.  If the torsion survives, the theta does not reach it.  The test is
run for the first theta found, and then for all of them together.

Note that the theta subgraph need not have paths of length 2: configuration spaces are
invariant under subdivision, so every theta graph has the same H_*(B_k), and a graph of girth
at least four (no two vertices with three common neighbours) needs paths of length three or
more.  The thetas are found by a max-flow search (one per pair of vertices that admits three
disjoint paths), not by looking for common neighbours.

WHAT IS PRINTED
---------------
    <graph>: N=.., M=.., k=..
      weight-k complex {degree: cells}            (unreduced)
      baseline torsion of coker(d_3) = Z/2
      <t> theta subgraphs found
      a single theta subgraph, hubs u and w: its edges, its own H_*(B_k theta) (torsion-free),
        the rank of Z_2(theta), and "torsion after adjoining it: ..."
      all <t> thetas together: "torsion after adjoining all of them: ..."
      a one-line summary and a verdict paragraph.
Progress lines from the cokernel computation ("... cols -> peeled ..., residual ...") are
interleaved.  Excerpt (the default, about 20 seconds):

      baseline torsion of coker(d_3) = Z/2
      15 theta subgraphs found (one per vertex pair admitting three disjoint paths)
      a single theta subgraph, hubs 0 and 1:
        edges [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4)]
        H_*(B_4 theta) = H_1=Z^3  H_2=Z  H_3=0  H_4=0
        Z_2(theta) has rank 1051; 1051 nonzero pushed columns  [1s]
        torsion after adjoining it: 0
      baseline Z/2;  one theta -> 0;  all thetas -> 0
      A single embedded theta subgraph accounts for the whole torsion. ...

USAGE
-----
    python3 experiments/theta_class.py [--graph octahedron] [--k 4] [--max-thetas N]

    --graph NAME     any zero-argument graph constructor in gbt (default octahedron; e.g.
                     dodecahedron, prism, tetrahedron).
    --k K            number of particles (default 4).
    --max-thetas N   stop after N theta subgraphs (default: all of them).

COST, OUTPUT FILES AND EXIT CODE
--------------------------------
About 20 seconds for the octahedron; graphs with more vertices take much longer (the unreduced
complex is used throughout, and the kernel lattice of each theta is the slow part).  Writes no
data file.  Exits with 1 if the graph has torsion but no theta subgraph, and 0 otherwise --
including when the thetas do NOT account for the torsion, which is reported in the text.
Requires networkx.
"""

import argparse
import itertools
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import gbt                                                            # noqa: E402
from gbt import Graph, swiatkowski, primary, fmt                       # noqa: E402
from gbt import fast                                                  # noqa: E402
from gbt.linalg import coker_torsion, indexed_columns, kernel_basis    # noqa: E402


def theta_subgraphs(G, limit=None):
    """Theta subgraphs of G: a list of ((u, w), edge list), one per vertex pair u < w that is
    joined by at least three internally vertex-disjoint paths.

    Uses a max-flow computation (networkx's node_disjoint_paths) rather than looking for
    common neighbours, so it finds thetas whose paths are longer than 2 -- which is what
    graphs of girth at least four, such as the dodecahedron, need.  Only the three shortest
    of the (at most four) disjoint paths found for each pair are used, which keeps the
    subgraph small.  The edges are given in G's vertex labels as sorted pairs.  `limit`
    stops the search after that many thetas.
    """
    import networkx as nx
    Gx = nx.Graph(G.edges)
    out = []
    for u, w in itertools.combinations(sorted(Gx.nodes()), 2):
        try:
            paths = list(itertools.islice(nx.node_disjoint_paths(Gx, u, w), 4))
        except Exception:                   # no path at all (networkx raises): skip pair
            continue
        if len(paths) < 3:
            continue
        paths = sorted(paths, key=len)[:3]
        edges = []
        for pth in paths:
            edges += [(pth[i], pth[i + 1]) for i in range(len(pth) - 1)]
        edges = sorted({tuple(sorted(e)) for e in edges})
        # a genuine theta: three internally disjoint paths, so no vertex other than the
        # hubs has degree 3 in the subgraph
        out.append(((u, w), edges))
        if limit and len(out) >= limit:
            break
    return out


def push_columns(sub_edges, eG, G, index2, k):
    """Z-basis of Z_2(subgraph, weight k), pushed into C_2(G) as columns.

    Arguments: the subgraph's edges `sub_edges` (pairs in G's labels), the map
    eG {frozenset edge: index in G}, the graph G, `index2` {degree-2 generator of G's
    unreduced complex: row index} as returned by indexed_columns(C, 3), and k.

    The subgraph is built on G's own vertex labels (unused vertices are isolated), so the
    inclusion is the identity on vertices and only the edge indices have to be translated.
    Both complexes are UNREDUCED, which is what makes the inclusion a chain map.

    Returns (columns, rank of Z_2, sizes of the subgraph's complex); each column is a dict
    {row index of G's degree-2 generators: coefficient}, zero columns dropped.
    """
    TH = Graph(G.N, sub_edges)
    eT = {i: frozenset(e) for i, e in enumerate(TH.edges)}
    CT = swiatkowski(TH, k)
    cols2, _, gens2_sub = indexed_columns(CT, 2)
    Z2 = kernel_basis(cols2, len(CT.gens[1]))

    def push(gen):
        """Image in G's complex of a basis element (mono, st) of the subgraph's complex:
        parked edges and moving-edge states are renamed from subgraph to G edge indices."""
        mono, st = gen
        nst = [0] * G.N
        for v, x in enumerate(st):
            if x == 1:
                nst[v] = 1
            elif x >= 2:                    # moving along subgraph edge x-2
                nst[v] = 2 + eG[eT[x - 2]]
        return (tuple(sorted(eG[eT[e]] for e in mono)), tuple(nst))

    out = []
    for kv in Z2:
        col = {}
        for j, c in kv.items():
            r = index2[push(gens2_sub[j])]
            col[r] = col.get(r, 0) + c
            if col[r] == 0:
                del col[r]
        if col:
            out.append(col)
    return out, len(Z2), CT.sizes()


def main():
    """Run the test of the module docstring on --graph at --k particles and print the verdict.

    Returns 1 if the graph has torsion but no theta subgraph, otherwise 0."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", default="octahedron",
                   help="any zero-argument graph constructor in gbt, e.g. octahedron, "
                        "dodecahedron, prism, tetrahedron")
    ap.add_argument("--k", type=int, default=4)
    ap.add_argument("--max-thetas", type=int, default=None,
                   help="cap on how many theta subgraphs to push (default: all)")
    args = ap.parse_args()
    k = args.k

    G = getattr(gbt, args.graph)()
    eG = {frozenset(e): i for i, e in enumerate(G.edges)}
    print("%s: N=%d, M=%d, k=%d" % (args.graph, G.N, G.M, k), flush=True)

    C = swiatkowski(G, k)
    print("  weight-%d complex %s" % (k, C.sizes()), flush=True)
    # the rows of d_3 are the degree-2 generators: index2 is where pushed cycles land
    d3, index2, _ = indexed_columns(C, 3)
    baseline = coker_torsion(d3, "%s d_3" % args.graph)
    print("  baseline torsion of coker(d_3) = %s\n" % (primary(baseline) or "none"),
          flush=True)
    if not baseline:
        print("  no torsion here, so there is nothing to explain.")
        print("  (run with --graph octahedron or dodecahedron, and --k 4)")
        return 0

    thetas = theta_subgraphs(G, limit=args.max_thetas)
    print("  %d theta subgraphs found (one per vertex pair admitting three disjoint paths)"
          % len(thetas), flush=True)
    if not thetas:
        print("  none: no theta subgraph, so this mechanism cannot be the explanation.")
        return 1

    # --- one theta on its own ------------------------------------------------
    (u, w), edges = thetas[0]
    print("\n  a single theta subgraph, hubs %d and %d:" % (u, w))
    print("    edges %s" % (edges,))
    relabel = {v: i for i, v in enumerate(sorted({x for e in edges for x in e}))}
    Hth, _, _ = fast.compute(
        Graph(len(relabel), [(relabel[a], relabel[b]) for a, b in edges]), k)
    print("    H_*(B_%d theta) = %s"
          % (k, "  ".join("H_%d=%s" % (d, fmt(*Hth[d]))
                          for d in sorted(Hth) if d > 0)))
    print("    -- torsion-free, as An--Knudsen's atoms must be", flush=True)

    t0 = time.time()
    cols, nz, szs = push_columns(edges, eG, G, index2, k)
    print("    theta weight-%d complex %s" % (k, szs))
    print("    Z_2(theta) has rank %d; %d nonzero pushed columns  [%.0fs]"
          % (nz, len(cols), time.time() - t0), flush=True)
    # adjoin the pushed cycles as extra columns: torsion of C_2 / (B_2 + i(Z_2 theta))
    after_one = coker_torsion(d3 + cols, "+ i_*(Z_2 of one theta)")
    print("    torsion after adjoining it: %s" % (primary(after_one) or "0"), flush=True)

    # --- and all of them -----------------------------------------------------
    allcols = []
    for (_, edges2) in thetas:
        c, _, _ = push_columns(edges2, eG, G, index2, k)
        allcols.extend(c)
    print("\n  all %d thetas together: %d pushed columns" % (len(thetas), len(allcols)),
          flush=True)
    after_all = coker_torsion(d3 + allcols, "+ i_*(all thetas)")
    print("  torsion after adjoining all of them: %s" % (primary(after_all) or "0"),
          flush=True)

    print()
    print("  baseline %s;  one theta -> %s;  all thetas -> %s"
          % (primary(baseline), primary(after_one) or "0", primary(after_all) or "0"))
    if not after_one:
        print()
        print("  A single embedded theta subgraph accounts for the whole torsion.  Its own")
        print("  H_2 is free, so a FREE class on the theta maps to the TORSION class on")
        print("  %s: the torsion is created by the relations the embedding" % args.graph)
        print("  imposes, not carried by a new generator.")
    elif not after_all:
        print()
        print("  The theta subgraphs collectively account for the torsion, though no single")
        print("  one does on its own.")
    else:
        print()
        print("  The theta classes do NOT account for the torsion; %s survives."
              % primary(after_all))
    return 0


if __name__ == "__main__":
    sys.exit(main())
