#!/usr/bin/env python3
"""
look_at_the_complex.py -- see the Swiatkowski complex itself: its cells, its differential, and
what partial reduction does to it.

WHAT IT IS FOR
--------------
The fast pipelines turn every cell into an integer and never show it to you.  When you want to
understand a computation -- check a sign by hand, see which cells a cycle uses, or convince
yourself that the model is right -- you need the cells in readable form.  This script builds
the complex of a small graph with the reference pipeline (gbt.swiatkowski + gbt.ChainComplex),
prints every generator in words together with its boundary, checks d o d = 0, computes the
homology, and then repeats with the complex reduced at every vertex.

The default graph is the tripod K_{1,3} with two particles.  Its configuration space B_2 K_{1,3}
is homotopy equivalent to a circle: two particles at the centre of a "Y" can trade places by
using the three arms, and that exchange loop generates H_1 = Z.  You should see exactly that.

USAGE
-----
    python3 examples/look_at_the_complex.py [GRAPH] [K] [--max-cells N]

GRAPH and K are as in examples/compute.py (default: the tripod "0-1 0-2 0-3" and K = 2).
--max-cells N (default 60) limits how many generators are printed per degree; the homology is
computed from the whole complex regardless.

READING THE OUTPUT
------------------
A generator of the unreduced complex is printed as, for example,

    [deg 1]  parked: e2={0,3} | moving: 0 along e1={0,2}

meaning: one particle is parked in the interior of edge 2 (the edge {0,3}), and a particle at
vertex 0 is moving out along edge 1 (this is the generator h_{0,e1}).  "sitting: 3" means a
particle sits at vertex 3.  The degree is the number of moving particles; the weight (number of
particles) is parked + sitting + moving and always equals K.  The boundary is printed as a
signed sum of such generators: d h_{v,e} = e - v, i.e. "the moving particle either finishes
its move into edge e (parked on e) or stays at v (sitting)", with the Koszul sign (-1)^(i-1)
for the i-th moving vertex in increasing vertex order.

In the reduced complex the moving generators are the differences h~_{v,e} = h_{v,e} - h_{v,e0},
where e0 is the edge of least index at v, and sitting particles no longer occur at reduced
vertices.  They are printed as "moving: v along e (minus e0)".

Nothing is written to data/.  Takes a second for the default graph.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import swiatkowski                                                 # noqa: E402
from gbt.report import homology_line                                        # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compute import parse_graph                                             # noqa: E402


def edge_name(G, e):
    """'e3={1,4}': edge index and endpoints, the way edges are referred to in the output."""
    u, v = G.edges[e]
    return "e%d={%d,%d}" % (e, u, v)


def describe_cell(G, cell, reduced_at=()):
    """A generator (mono, st) of the Swiatkowski complex in words.

    `cell` is the pair (mono, st) used throughout gbt: `mono` is a sorted tuple of edge indices
    (one entry per parked particle, with repetition), and `st[v]` is 0 (v empty), 1 (a particle
    sits at v) or 2 + e (a particle at v moves along edge e).  `reduced_at` lists the vertices
    at which the complex was reduced; there a moving state means the difference generator
    h_{v,e} - h_{v,e0(v)}, and the text says so.
    """
    mono, st = cell
    parts = []
    if mono:
        parts.append("parked: " + ", ".join(edge_name(G, e) for e in mono))
    sitting = [v for v, s in enumerate(st) if s == 1]
    if sitting:
        parts.append("sitting: " + ", ".join(str(v) for v in sitting))
    moving = []
    for v, s in enumerate(st):
        if s >= 2:
            text = "%d along %s" % (v, edge_name(G, s - 2))
            if v in reduced_at:
                e0 = min(i for i, (a, b) in enumerate(G.edges) if v in (a, b))
                text += " (minus e%d)" % e0
            moving.append(text)
    if moving:
        parts.append("moving: " + ", ".join(moving))
    return " | ".join(parts) if parts else "(empty: no particles)"


def show(G, k, reduce_at, max_cells):
    """Build the complex reduced at `reduce_at`, print generators with boundaries, check
    d o d = 0 and print the homology.  Returns the homology dict."""
    C = swiatkowski(G, k, reduce_at=reduce_at)
    print("generators by degree: %s" % C.sizes())
    print("d o d = 0: %s" % C.check_d2())
    for d in sorted(C.gens):
        gens = sorted(C.gens[d], key=repr)          # a deterministic order for printing
        print("\n  degree %d (%d generators)" % (d, len(gens)))
        for c in gens[:max_cells]:
            print("    [deg %d]  %s" % (d, describe_cell(G, c, reduce_at)))
            if d > 0:
                bd = C.B[d].get(c, {})
                terms = ["%+d * (%s)" % (coef, describe_cell(G, f, reduce_at))
                         for f, coef in sorted(bd.items(), key=lambda kv: repr(kv[0]))]
                print("             d = " + ("\n                 ".join(terms) if terms else "0"))
        if len(gens) > max_cells:
            print("    ... %d more (raise --max-cells to see them)" % (len(gens) - max_cells))
    # homology() works on any complex; reducing first (discrete Morse theory) only makes it
    # faster.  reduce_fully() modifies C in place, which is fine since we are done printing.
    C.reduce_fully()
    H = C.homology()
    print("\nhomology: H_0=%s  %s" % ("Z^%d" % H[0][0] if H[0][0] > 1 else "Z",
                                      homology_line(H)))
    return H


def main():
    """Show the unreduced complex, then the complex reduced at every vertex, and compare."""
    ap = argparse.ArgumentParser()
    ap.add_argument("graph", nargs="?", default="0-1 0-2 0-3")
    ap.add_argument("k", nargs="?", type=int, default=2)
    ap.add_argument("--max-cells", type=int, default=60)
    args = ap.parse_args()
    G = parse_graph(args.graph)
    print("graph: %d vertices, edges %s;  k = %d"
          % (G.N, ", ".join(edge_name(G, e) for e in range(G.M)), args.k))

    print("\n=========== the UNREDUCED complex (reduce_at = ()) ===========")
    H1 = show(G, args.k, (), args.max_cells)

    # A vertex of degree 1 has no non-basepoint edge, so it cannot be reduced (see
    # gbt.swiatkowski); reduce only at the vertices of degree >= 2.
    R = tuple(v for v in range(G.N) if G.deg(v) >= 2)
    print("\n=========== REDUCED at the vertices %s ===========" % (R,))
    H2 = show(G, args.k, R, args.max_cells)

    print("\nsame homology: %s" % (H1 == H2))
    return 0 if H1 == H2 else 1


if __name__ == "__main__":
    sys.exit(main())
