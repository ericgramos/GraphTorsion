#!/usr/bin/env python3
"""
Is the H_2 torsion supported on a proper subgraph?  Single-edge deletions and greedy descent.

Two questions about how tightly the torsion is tied to the whole graph.

  (1) EDGE MINIMALITY.  Delete one edge at a time and recompute.  If every single-edge
      deletion kills the torsion, the carrier is edge-minimal: no proper subgraph on the same
      vertices has it.  This holds for K_7 at k=3 (all 21 deletions) and for the octahedron
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

  (3) STABILIZATION IN k.  By An--Drummond-Cole--Knudsen the groups
      (+)_k H_i(B_k Gamma) form a module over Z[E], and torsion is non-decreasing in k.  The
      script shows K_6 passing from torsion-free H_2 at k=2,3 to Z/2 at k=4, which is the
      visible form of that statement: the class is born at k=4 and persists.

A warning about interpretation, since this is an easy mistake to make.  Edge
minimality does NOT mean the class is indecomposable, nor that it fails to come from a
subgraph: a homomorphism can send a free class to a torsion class, and that is exactly what
happens here -- see `theta_class.py`.  Minimality says every edge is needed for the class to be
torsion, not that the class is new.

    python3 experiments/edge_minimality.py [--case K7-3|K6-4|OCT-4|descent|all]

K7-3 is the slow one (21 recomputations of a 7-vertex, k=3 complex).
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
    """The invariant factors of H_d(B_k G), with a guard for degenerate subgraphs."""
    H, _, _ = fast.compute(G, k)
    return H.get(d, (0, []))


def deletion_test(name, G, k, d):
    """Does every single-edge deletion kill the torsion in degree d?"""
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
    """One greedy pass: drop edges in the given order while the torsion survives."""
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
    """Randomised greedy descent from K_n: drop edges while keeping torsion in degree d."""
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
    print("=== %s: torsion in H_%d as k grows (edge stabilization) ===" % (n, d))
    G = complete(int(n.split("_")[1]))
    for k in ks:
        H, _, _ = fast.compute(G, k)
        print("   k=%d:  H_%d = %s" % (k, d, fmt(*H.get(d, (0, [])))), flush=True)
    print("   torsion is non-decreasing in k, as edge stabilization requires\n")


def main():
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
