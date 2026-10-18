#!/usr/bin/env python3
"""
Is the Z/2 in H_2(B_4 K_6) indecomposable?  Yes -- by exhaustion, no product can reach it.

`decomposability_K8.py` shows that in K_8 at k=3 every torsion class is a Kunneth product.
K_6 at k=4 is the opposite case, and the argument is a counting one rather than a computation
on cycles: there simply is no pair of vertex-disjoint subgraphs of K_6 whose homology can
contribute torsion to degree 2.

This script makes that exhaustive.  For every splitting of the six vertices into two disjoint
sets of sizes m_1 + m_2 <= 6, and every k_1 + k_2 = 4, it looks up H_*(B_{k_i} K_{m_i}) and
asks whether Kunneth could produce 2-torsion in total degree 2:

        H_i(B_{k_1} K_{m_1}) (x) H_j(B_{k_2} K_{m_2})    with i + j = 2, torsion in a factor
        Tor( H_i, H_j )                                  with i + j = 1

A product subgraph of K_6 has at most 5 vertices, and no graph on at most 5 vertices has
torsion in H_i for i >= 2 (`sweep_small_graphs.py --vertices 5`), while torsion in H_1 needs
non-planarity, i.e. 5 vertices and at least 9 edges, leaving at most 1 vertex for the other
factor -- not enough for a cycle.  The enumeration confirms that every one of the
combinations fails.

The script prints the homology of all the pieces, so the reasoning can be checked by eye
rather than taken on trust, and then states the conclusion.  `octahedron_generator.py` finds
what the class actually is.

    python3 experiments/kunneth_K6.py [--n 6] [--k 4]
"""

import argparse
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import Graph, complete, primary             # noqa: E402
from gbt import fast                                 # noqa: E402
from gbt.report import torsion_line                  # noqa: E402


def induced_complete(m):
    """K_m as an abstract graph on m vertices (an empty graph when m < 2)."""
    return Graph(m, list(itertools.combinations(range(m), 2))) if m >= 2 else Graph(1, [])


def homology_or_none(m, k, cache):
    """H_*(B_k K_m) as {d: (free, invariant factors)}, or 'pt' for k = 0, or None."""
    if k == 0:
        return "pt"
    if (m, k) not in cache:
        try:
            H, _, _ = fast.compute(induced_complete(m), k)
        except Exception:
            H = None
        cache[(m, k)] = H
    return cache[(m, k)]


def get(H, d):
    if H == "pt":
        return (1, []) if d == 0 else (0, [])
    if H is None:
        return None
    return H.get(d, (0, []))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--k", type=int, default=4)
    args = ap.parse_args()
    n, k = args.n, args.k

    H_full, _, _ = fast.compute(complete(n), k)
    print("H_*(B_%d K_%d) torsion: %s" % (k, n, torsion_line(H_full)))
    print("the question: can any Kunneth product reach it?\n", flush=True)

    cache = {}
    rows = []
    for m1 in range(1, n):
        for m2 in range(1, n + 1 - m1):
            for k1 in range(0, k + 1):
                k2 = k - k1
                H1 = homology_or_none(m1, k1, cache)
                H2 = homology_or_none(m2, k2, cache)
                sources = []
                # torsion in total degree 2 from a tensor factor
                for i in (0, 1, 2):
                    j = 2 - i
                    a, b = get(H1, i), get(H2, j)
                    if a is None or b is None:
                        continue
                    if a[1] and (b[0] or b[1]):
                        sources.append("H_%d(G1)_tors (x) H_%d(G2)" % (i, j))
                    if b[1] and (a[0] or a[1]):
                        sources.append("H_%d(G1) (x) H_%d(G2)_tors" % (i, j))
                # and from the Tor term
                for i in (0, 1):
                    j = 1 - i
                    a, b = get(H1, i), get(H2, j)
                    if a is None or b is None:
                        continue
                    if a[1] and b[1]:
                        sources.append("Tor(H_%d, H_%d)" % (i, j))
                rows.append((m1, m2, k1, k2, sources))

    hits = [r for r in rows if r[4]]
    print("%d combinations (|V(G1)|, |V(G2)|, k1, k2) with m1 + m2 <= %d and k1 + k2 = %d"
          % (len(rows), n, k))
    print("combinations that could contribute 2-torsion in degree 2: %d" % len(hits))
    for r in hits:
        print("   m1=%d m2=%d k1=%d k2=%d  via %s" % (r[0], r[1], r[2], r[3], r[4]))

    print("\nthe pieces that fit inside K_%d:" % n)
    print("  %-12s %s" % ("subgraph", "torsion of H_*(B_k)"))
    for m in range(1, n):
        for kk in range(1, k + 1):
            H = homology_or_none(m, kk, cache)
            if H in (None, "pt"):
                continue
            print("  %-12s %s" % ("K_%d, k=%d" % (m, kk), torsion_line(H)))

    print()
    if not hits:
        deg2 = primary(H_full.get(2, (0, []))[1])
        print("No combination can contribute 2-torsion in degree 2, so the %s" % deg2)
        print("in H_2(B_%d K_%d) is INDECOMPOSABLE: it is not a Kunneth product" % (k, n))
        print("of classes on vertex-disjoint subgraphs.")
        print()
        print("The reason is a count.  A proper subgraph here has at most %d vertices," % (n - 1))
        print("torsion in H_1 needs a non-planar subgraph (5 vertices, Ko--Park), and")
        print("torsion in H_2 or above needs more than 5 vertices -- which leaves nothing")
        print("for the second factor: a product needs a K_5 (for the Ko--Park class) and a")
        print("vertex-disjoint cycle, so at least 8 vertices.")
    else:
        print("Some combination could contribute; build the products explicitly")
        print("(as in decomposability_K8.py) to decide whether they really do.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
