#!/usr/bin/env python3
"""
Is the H_2 torsion supported on a proper subgraph?  Single-edge deletions and greedy descent.

Two questions about how tightly the torsion is tied to the whole graph.

  (1) EDGE MINIMALITY.  Delete one edge at a time and recompute.  If every single-edge
      deletion kills the torsion, the carrier is called edge-minimal.  (This is a statement
      about single deletions only.  Having torsion is not known to be monotone under edge
      deletion, so it does not by itself rule out a smaller subgraph, two or more edges down,
      that carries torsion again; the exhaustive sweeps of `sweep_small_graphs.py` are what
      settle such questions.)  This holds for K_7 at k=3 (all 21 deletions) and for the octahedron
      at k=4 (all 12).  It does NOT hold for K_6 at k=4: every one of the 15 deletions leaves
      a graph that still has Z/2 in H_2, which is the first hint that the K_6 class really
      lives on a smaller subgraph -- the greedy descent in (2) then finds it.

  (2) GREEDY DESCENT.  Starting from K_6 at k=4, repeatedly drop any edge that can be dropped
      while keeping torsion in H_2, until none can.  This locates a LOCALLY minimal carrier
      inside a given graph, and is how the octahedron was first found.  Greedy descent is
      order-dependent: a single pass in the natural edge order gets stuck at 13 edges, while
      shuffling the order reaches the 12-edge 4-regular graph -- the octahedron -- so the
      script makes several randomised restarts and reports the best. The authority on
      minimality is not this search but the exhaustive six-vertex sweep in
      `sweep_small_graphs.py`, which shows 12 edges is the true minimum and the octahedron
      the unique graph achieving it.

  (3) STABILIZATION IN k.  By An--Drummond-Cole--Knudsen the groups (+)_k H_i(B_k Gamma)
      form a module over Z[E], where multiplying by an edge e adds a particle on e (edge
      stabilization).  Edge stabilization H_i(B_k) -> H_i(B_{k+1}) is injective
      (An--Drummond-Cole--Knudsen, "Subdivisional spaces and graph braid groups",
      Doc. Math. 24 (2019), Proposition 5.21), so torsion present at k persists at every
      larger k.  The script shows K_6 passing from
      torsion-free H_2 at k=2,3 to Z/2 at k=4: the class is born at k=4.  (It only PRINTS
      the groups; it does not test persistence beyond k = 4.)

A warning about interpretation, since this is an easy mistake to make.  Edge
minimality does NOT mean the class is indecomposable, nor that it fails to come from a
subgraph: a homomorphism can send a free class to a torsion class, and that is exactly what
happens here -- see `theta_class.py`.  Minimality says every edge is needed for the class to be
torsion, not that the class is new.

WHAT IS PRINTED
---------------
For a deletion test (cases OCT-4, K6-4, K7-3):

    === NAME: H_d = <group> ===
       NAME - (u, v):  H_d = <group after deleting the first edge>      (a sample)
       T of M single-edge deletions tested (S skipped: would leave a leaf)
       every deletion kills the H_d torsion => EDGE-MINIMAL
    or
       torsion SURVIVES deleting: [((u, v), [invariant factors]), ...]
       => not edge-minimal

A deletion that would leave a vertex of degree < 2 is skipped (none of the three graphs tested
has such an edge).  Real output for the octahedron (3 seconds):

    === octahedron: H_2 = Z^95 + Z/2 ===
       octahedron - (0, 2):  H_2 = Z^65
       12 of 12 single-edge deletions tested (0 skipped: would leave a leaf)
       every deletion kills the H_2 torsion => EDGE-MINIMAL

For the descent: one line per restart, "restart r: stops at m edges, degrees [...]", then the
best graph found (edges, degree sequence, torsion) and, if it is 4-regular on 6 vertices with
12 edges, the remark that it is the octahedron.  For the stabilization case: H_2(B_k K_6) for
k = 2, 3, 4, e.g.

    === K_6: torsion in H_2 as k grows (edge stabilization) ===
       k=2:  H_2 = Z^9
       k=3:  H_2 = Z^129
       k=4:  H_2 = Z^264 + Z/2

USAGE
-----
    python3 experiments/edge_minimality.py [--case K7-3|K6-4|OCT-4|descent|stab|all]

    --case   (default all) which part to run:
             OCT-4    deletions from the octahedron at k = 4            (seconds)
             K6-4     deletions from K_6 at k = 4                        (about 15 seconds)
             stab     H_2(B_k K_6) for k = 2, 3, 4                       (a few seconds)
             descent  six greedy descents from K_6 at k = 4              (about 40 seconds)
             K7-3     deletions from K_7 at k = 3, 21 recomputations     (about 7 seconds)
             all      all of the above, in the order OCT-4, K6-4, stab, descent, K7-3

Every computation uses `gbt.fast`; memory stays well under 1 GB.  Nothing is written to
data/.  Exit code: always 0 (the findings are printed, not asserted).
"""

import argparse
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import Graph, complete, octahedron, fmt                 # noqa: E402
from gbt import fast                                            # noqa: E402
from gbt.report import torsion_line                             # noqa: E402


def torsion_in_degree(G, k, d):
    """(free rank, invariant factors > 1) of H_d(B_k G), computed in full with `gbt.fast`.

    Degrees above the top of the complex are absent from the homology dict; they are reported
    as (0, []).  The callers make sure G has no vertex of degree < 2 (fast needs degree >= 1).
    """
    H, _, _ = fast.compute(G, k)
    return H.get(d, (0, []))


def deletion_test(name, G, k, d):
    """Does every single-edge deletion kill the torsion in degree d?

    name    label used in the printout;  G  the graph;  k  particles;  d  homological degree.
    Deletions that would leave a vertex of degree < 2 are skipped and counted.
    Returns None if G itself has no torsion in degree d, True if every tested deletion kills
    it (edge-minimal), False if some deletion keeps it (the survivors are printed).
    """
    free, tors = torsion_in_degree(G, k, d)
    print("=== %s: H_%d = %s ===" % (name, d, fmt(free, tors)), flush=True)
    if not tors:
        print("   no torsion here; nothing to test\n")
        return None
    survivors = []
    skipped = 0
    for i in range(G.M):
        sub = G.delete_edge(i)
        if min(sub.deg(v) for v in range(sub.N)) < 2:
            # a degree-0 or degree-1 vertex changes the homotopy type for a trivial
            # reason (fast.compute also needs every degree >= 1), so these are not
            # informative deletions
            skipped += 1
            continue
        f2, t2 = torsion_in_degree(sub, k, d)
        if t2:
            survivors.append((G.edges[i], t2))
        if i == 0:
            print("   %s - %s:  H_%d = %s" % (name, str(G.edges[i]), d, fmt(f2, t2)),
                  flush=True)
    tested = G.M - skipped
    print("   %d of %d single-edge deletions tested (%d skipped: would leave a leaf)"
          % (tested, G.M, skipped))
    if survivors:
        print("   torsion SURVIVES deleting: %s" % survivors)
        print("   => not edge-minimal\n", flush=True)
        return False
    print("   every deletion kills the H_%d torsion => EDGE-MINIMAL\n" % d, flush=True)
    return True


def _descend(n, k, d, order):
    """One greedy descent: repeatedly drop the first droppable edge, in the given order.

    n  number of vertices;  k, d  particles and degree;  order  the starting edge list.
    An edge is droppable if the remaining graph still uses all n vertices, has minimum degree
    >= 2, and still has torsion in H_d.  After each drop the scan restarts from the beginning.
    Returns the edge list at which no further edge can be dropped (a local minimum).
    """
    cur = list(order)
    while True:
        dropped = False
        for i in range(len(cur)):
            cand = [e for j, e in enumerate(cur) if j != i]
            if len({x for e in cand for x in e}) < n:
                continue
            sub = Graph(n, cand)
            if min(sub.deg(v) for v in range(n)) < 2:
                continue
            _, t2 = torsion_in_degree(sub, k, d)
            if t2:
                cur = cand
                dropped = True
                break
        if not dropped:
            return cur


def greedy_descent(n, k, d, restarts=6, seed=0):
    """Randomised greedy descent from K_n: drop edges while keeping torsion in degree d.

    Restart 0 uses the natural edge order of K_n, restarts 1.. a shuffled order drawn from
    random.Random(seed), so the run is reproducible.  Prints each restart's stopping size and
    the best (smallest) graph found with its torsion.  Returns None.
    """
    print("=== greedy descent from K_%d at k=%d (%d randomised restarts) ==="
          % (n, k, restarts), flush=True)
    G = complete(n)
    rng = random.Random(seed)
    t0 = time.time()
    best = None
    for r in range(restarts):
        order = list(G.edges)
        if r:                                  # pass 0 uses the natural edge order
            rng.shuffle(order)
        cur = _descend(n, k, d, order)
        sub = Graph(n, cur)
        degs = sorted((sub.deg(v) for v in range(n)), reverse=True)
        print("   restart %d: stops at %d edges, degrees %s" % (r, len(cur), degs),
              flush=True)
        if best is None or len(cur) < len(best):
            best = cur
    sub = Graph(n, best)
    H, _, _ = fast.compute(sub, k)
    degs = sorted((sub.deg(v) for v in range(n)), reverse=True)
    print("   BEST: %d edges, degree sequence %s   [%.0fs]"
          % (len(best), degs, time.time() - t0))
    print("   edges: %s" % sorted(best))
    print("   torsion: %s" % torsion_line(H))
    if len(best) == 12 and degs == [4] * 6:
        print("   (4-regular on 6 vertices with 12 edges: this is the octahedron)")
    print("   NOTE: greedy descent only finds LOCAL minima -- different restarts stop at")
    print("   different sizes.  sweep_small_graphs.py --vertices 6 --k 4 settles the true")
    print("   minimum exhaustively.")
    print(flush=True)


def stabilization(n, ks, d):
    """Print H_d(B_k K_m) for each k in ks, where n is a label of the form "K_m".

    Only prints the groups; the closing sentence about non-decreasing torsion is the theorem
    the reader should see illustrated, not something this function checks.
    """
    print("=== %s: torsion in H_%d as k grows (edge stabilization) ===" % (n, d))
    G = complete(int(n.split("_")[1]))
    for k in ks:
        H, _, _ = fast.compute(G, k)
        print("   k=%d:  H_%d = %s" % (k, d, fmt(*H.get(d, (0, [])))), flush=True)
    print("   torsion is non-decreasing in k, as edge stabilization requires\n")


def main():
    """Run the case(s) selected by --case; always return 0."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="all",
                    choices=("K7-3", "K6-4", "OCT-4", "descent", "stab", "all"))
    args = ap.parse_args()
    c = args.case

    if c in ("OCT-4", "all"):
        deletion_test("octahedron", octahedron(), 4, 2)
    if c in ("K6-4", "all"):
        deletion_test("K_6", complete(6), 4, 2)
    if c in ("stab", "all"):
        stabilization("K_6", (2, 3, 4), 2)
    if c in ("descent", "all"):
        greedy_descent(6, 4, 2)
    if c in ("K7-3", "all"):
        deletion_test("K_7", complete(7), 3, 2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
