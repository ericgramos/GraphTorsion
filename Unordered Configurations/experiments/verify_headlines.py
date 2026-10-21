#!/usr/bin/env python3
"""
Independent recomputation of the headline torsion groups.

WHY
---
The stored results in data/ were mostly produced by `gbt.fast` (Morse reduction of the fully
reduced complex, then Smith normal form of the whole residual).  This script recomputes the most
important torsion groups by two other routes, and checks them against expected values
written into this script (the values stored in data/ and quoted in the README):

  (1) `gbt.lean.torsion_in_degree`: builds only the degree slices C_d and C_{d+1}, peels unit
      pivots without a priority queue, and takes the Smith form of what is left.  It uses
      tors H_d = tors coker(d_{d+1}), so it never looks at d_d.  It shares the complex builder
      with `gbt.fast` but none of the reduction or bookkeeping code.
  (2) Abrams' discretised configuration space UD^k of the subdivided graph, a different chain
      complex altogether, built here from scratch.  It shares no code with the Swiatkowski
      builders beyond the Graph class.  Each edge is subdivided enough to meet the
      Abrams--Prue--Scrimshaw hypotheses (paths between essential vertices of length >= k-1,
      cycles of length >= k+1).

The subdivision factor of each edge is chosen by hand per case (the "edges xT" in the output)
and is NOT checked by the code.

WHAT IS CHECKED
---------------
Route 1 (23 jobs), all torsion of H_2 unless stated: K_7 to K_10 at k = 3; K_6 and K_7 at
k = 4; K_5 and K_{3,3} at k = 5 (torsion-free); K_{4,5} at k = 3; K_{5,5} at k = 4 (H_2 and
H_3); K_{5,6} at k = 4 (H_3); K_{4,5} plus three edges inside the large side at k = 3; the
octahedron at k = 3 and 4; the triangular and pentagonal prisms at k = 4; K_6 and the Petersen
graph of the Petersen family at k = 4, and the Petersen graph at k = 3; the apex graph of
`linkless_controls.py` at k = 3; the dodecahedron and icosahedron at k = 4.
Route 2 (6 jobs): K_5 and K_{3,3} at k = 2 (H_1), the octahedron at k = 3 and 4, K_{4,5}
and K_7 at k = 3.

WHAT IS PRINTED
---------------
One line per job,

    <route> <graph> k=<k>  tors H_<d> = <group>   ok | FAIL (expected ...)   <cells> cells  <time>

where cells is the number of cells actually built (the two slices for route 1, the two degrees
of Abrams' complex for route 2), then "ALL HEADLINE VALUES CONFIRMED by both routes." or
"<m> DISAGREEMENT(S): [...]".  Excerpt (`--skip-large`, about 40 seconds):

    == route 1: slice pipeline gbt.lean on the Swiatkowski complex ==
    lean    K_7                        k=3  tors H_2 = Z/2               ok        15,400 cells  0s
    ...
    == route 2: Abrams' discretised model UD^k of the subdivided graph ==
    Abrams  K_5 (edges x1)             k=2  tors H_1 = Z/2               ok            45 cells  0s
    ...
    ALL HEADLINE VALUES CONFIRMED by both routes.

USAGE, COST AND EXIT CODE
-------------------------
    python3 experiments/verify_headlines.py [--skip-large]

    --skip-large   drop the dodecahedron, the icosahedron and H_2 of K_{5,5} at k = 4.

About ten minutes in all, and a few GB for the icosahedron (1.3 million cells in two slices);
about 40 seconds with --skip-large.  Writes nothing to data/.  Exits with 1 if any value
disagrees with the expected one, else 0.
"""

import argparse
import itertools
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import (Graph, complete, complete_bip, dodecahedron, from_edges,         # noqa: E402
                 icosahedron, lean, octahedron, petersen_family, primary, prism_graph,
                 subdivide)
from gbt.linalg import peel_sweeps, sparse_smith                                 # noqa: E402
from gbt.sizes import reduced_cell_count_by_degree                               # noqa: E402


def maximal_planar_apex(n, seed):
    """apex(planar n) of linkless_controls.py: a stacked triangulation plus a cone vertex.

    Start from K_4, add vertices 4..n-1 each stacked into a face chosen with
    random.Random(seed), then join a new vertex n to every vertex.  Returns a gbt Graph on
    n+1 vertices (linklessly embeddable, not planar).  Must stay identical to the
    construction in linkless_controls.py, whose stored value it checks."""
    edges = {(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)}
    rng = random.Random(seed)
    faces = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    for v in range(4, n):
        f = faces.pop(rng.randrange(len(faces)))
        for x in f:
            edges.add((min(v, x), max(v, x)))
        faces += [(v, f[0], f[1]), (v, f[0], f[2]), (v, f[1], f[2])]
    apex = n
    edges |= {(v, apex) for v in range(n)}
    return from_edges(sorted(edges))


def k45_plus_three_in_B():
    """K_{4,5} on A = {0..3}, B = {4..8} plus the first three B--B edges, as in
    density_interpolation.py."""
    base = [(a, b) for a in range(4) for b in range(4, 9)]
    bb = list(itertools.combinations(range(4, 9), 2))[:3]
    return Graph(9, base + bb)


def abrams_torsion(G, k, d):
    """tors H_d(UD^k G) = tors coker(d_{d+1}) in Abrams' cubical complex, built here.

    G must already be subdivided enough for UD^k G to model B_k G (not checked).  A cell is
    a set of k "tokens" (vertices or edges of G) with pairwise disjoint closures; its degree
    is its number of edges.  Only the cells of degrees d and d+1 are enumerated.  The
    boundary of a cell replaces its j-th edge token (in token order) by either endpoint, with
    sign (-1)^j and + for the first endpoint u, - for the second.

    Returns (sorted invariant factors > 1 of the torsion, number of cells built)."""
    toks = [("v", u) for u in range(G.N)] + [("e", i) for i in range(G.M)]
    sup = [{t[1]} if t[0] == "v" else set(G.edges[t[1]]) for t in toks]
    cells = {d: [], d + 1: []}

    def rec(start, chosen, used, ne):
        """Enumerate token sets in increasing token order, pruning above d+1 edges; `used` is
        the set of vertices already covered, `ne` the number of edges chosen."""
        if ne > d + 1:
            return
        if len(chosen) == k:
            if ne in cells:
                cells[ne].append(tuple(chosen))
            return
        for i in range(start, len(toks)):
            if sup[i] & used:
                continue
            rec(i + 1, chosen + [i], used | sup[i], ne + (toks[i][0] == "e"))

    rec(0, [], set(), 0)
    idx = {c: i for i, c in enumerate(cells[d])}
    cols = []
    for c in cells[d + 1]:
        col = {}
        es = [i for i in c if toks[i][0] == "e"]
        for j, ei in enumerate(es):
            u, v = G.edges[toks[ei][1]]
            for w in (u, v):             # vertex token index == vertex number
                # the face is again a sorted token tuple, so it can be looked up in idx
                face = tuple(sorted((set(c) - {ei}) | {w}))
                s = (-1) ** j * (1 if w == u else -1)
                r = idx[face]
                col[r] = col.get(r, 0) + s
        cols.append({r: x for r, x in col.items() if x})
    res, _ = peel_sweeps(cols)
    _, divs = sparse_smith(res, verbose=False)
    return sorted(x for x in divs if x not in (0, 1)), len(cells[d]) + len(cells[d + 1])


def main():
    """Run both routes and compare with the expected values; returns 1 on any disagreement."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-large", action="store_true")
    args = ap.parse_args()
    pf = dict(petersen_family())

    # name, graph, k, degree, expected torsion, large?
    lean_jobs = [
        ("K_7", complete(7), 3, 2, [2], False),
        ("K_8", complete(8), 3, 2, [2] * 21, False),
        ("K_9", complete(9), 3, 2, [2] * 28, False),
        ("K_10", complete(10), 3, 2, [2] * 36, False),
        ("K_6", complete(6), 4, 2, [2], False),
        ("K_7", complete(7), 4, 2, [2, 2], False),
        ("K_5", complete(5), 5, 2, [], False),
        ("K_{3,3}", complete_bip(3, 3), 5, 2, [], False),
        ("K_{4,5}", complete_bip(4, 5), 3, 2, [2, 2, 2, 6], False),
        ("K_{5,5}", complete_bip(5, 5), 4, 2, [2] * 17, True),
        ("K_{5,5}", complete_bip(5, 5), 4, 3, [2], False),
        ("K_{5,6}", complete_bip(5, 6), 4, 3, [2, 2, 2, 2, 6], False),
        ("K_{4,5} + 3 edges in B", k45_plus_three_in_B(), 3, 2, [2] * 12 + [6], False),
        ("octahedron", octahedron(), 3, 2, [], False),
        ("octahedron", octahedron(), 4, 2, [2], False),
        ("triangular prism", prism_graph(3), 4, 2, [], False),
        ("pentagonal prism", prism_graph(5), 4, 2, [2], False),
        ("K_6", pf["K_6"], 4, 2, [2], False),
        ("Petersen", pf["Petersen"], 4, 2, [2], False),
        ("Petersen", pf["Petersen"], 3, 2, [], False),
        ("apex(planar 8)", maximal_planar_apex(8, 2), 3, 2, [2, 2, 2, 2], False),
        ("dodecahedron", dodecahedron(), 4, 2, [2], True),
        ("icosahedron", icosahedron(), 4, 2, [2], True),
    ]
    # name, graph, edge subdivision, k, degree, expected torsion
    abrams_jobs = [
        ("K_5", complete(5), 1, 2, 1, [2]),
        ("K_{3,3}", complete_bip(3, 3), 1, 2, 1, [2]),
        ("octahedron", octahedron(), 2, 3, 2, []),
        ("octahedron", octahedron(), 3, 4, 2, [2]),
        ("K_{4,5}", complete_bip(4, 5), 2, 3, 2, [2, 2, 2, 6]),
        ("K_7", complete(7), 2, 3, 2, [2]),
    ]

    bad = []

    def show(route, name, k, d, got, exp, cells, dt):
        """Print one result line and record a disagreement in `bad`."""
        ok = sorted(got) == sorted(exp)
        if not ok:
            bad.append((route, name, k, d))
        print("%-7s %-26s k=%d  tors H_%d = %-26s %-4s %11s cells %5.0fs"
              % (route, name, k, d, primary(got) if got else "0",
                 "ok" if ok else "FAIL (expected %s)" % (primary(exp) if exp else "0"),
                 format(cells, ","), dt), flush=True)

    print("== route 1: slice pipeline gbt.lean on the Swiatkowski complex ==")
    for name, G, k, d, exp, large in lean_jobs:
        if large and args.skip_large:
            continue
        sl = reduced_cell_count_by_degree(G, k)
        t0 = time.time()
        got = lean.torsion_in_degree(G, k, d)
        show("lean", name, k, d, got, exp, sl.get(d, 0) + sl.get(d + 1, 0), time.time() - t0)

    print("\n== route 2: Abrams' discretised model UD^k of the subdivided graph ==")
    for name, G, t, k, d, exp in abrams_jobs:
        t0 = time.time()
        got, cells = abrams_torsion(subdivide(G, t), k, d)
        show("Abrams", "%s (edges x%d)" % (name, t), k, d, got, exp, cells, time.time() - t0)

    print()
    if bad:
        print("%d DISAGREEMENT(S): %s" % (len(bad), bad))
        return 1
    print("ALL HEADLINE VALUES CONFIRMED by both routes.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
