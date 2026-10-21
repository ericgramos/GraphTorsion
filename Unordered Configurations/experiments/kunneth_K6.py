#!/usr/bin/env python3
"""
Is the Z/2 in H_2(B_4 K_6) indecomposable?  Yes -- by exhaustion, no product can reach it.

THE QUESTION
------------
`decomposability_K8.py` finds that in K_8 at k = 3 the torsion of H_2 disappears once the
Kunneth product classes are adjoined.  K_6 at k = 4 is the opposite case, and the argument is a
counting one rather than a computation on cycles: there simply is no pair of vertex-disjoint
subgraphs of K_6 whose homology can contribute torsion to degree 2.  So the Z/2 of
H_2(B_4 K_6) is not a Kunneth product, and something else must explain it
(`octahedron_generator.py` shows it is the image of the octahedron's class).

THE METHOD
----------
A product class comes from vertex-disjoint subgraphs G_1, G_2 on m_1 and m_2 vertices with
k_1 + k_2 = k particles.  For every (m_1, m_2) with m_1 + m_2 <= n and every k_1 + k_2 = k, the
script takes the LARGEST possible pieces, the complete graphs K_{m_1} and K_{m_2}, computes
their homology with `gbt.fast`, and asks whether Kunneth could produce torsion in total degree 2:

        H_i(B_{k_1} K_{m_1}) (x) H_j(B_{k_2} K_{m_2})    with i + j = 2, torsion in a factor
                                                         and the other factor nonzero
        Tor( H_i, H_j )                                  with i + j = 1, torsion in both

A factor with 0 particles is a point (H_0 = Z).  The reason this is enough is the vertex count
below, not the choice of complete pieces (see PITFALLS).

THE COUNT
---------
A product subgraph of K_6 has at most 5 vertices, and no graph on at most 5 vertices has
torsion in H_i for i >= 2 (`sweep_small_graphs.py`, all graphs on 4 and 5 vertices).  Torsion
in H_1 needs a non-planar graph (Ko--Park), hence a K_5 or K_{3,3} minor, hence at least 5
vertices; that leaves at most 1 vertex for the other factor -- not enough for a cycle, so the
other factor has no H_1, and H_1 (x) H_1 and Tor(H_1, H_0) give nothing either.  The
enumeration confirms that all 75 combinations fail.

WHAT IS PRINTED
---------------
    H_*(B_k K_n) torsion: ...              the group being explained
    C combinations (|V(G1)|, |V(G2)|, k1, k2) with m1 + m2 <= n and k1 + k2 = k
    combinations that could contribute 2-torsion in degree 2: H
       m1=.. m2=.. k1=.. k2=..  via [...]  (one line per possible contribution, if any)

then the torsion of every piece K_m with m < n and 1 <= k' <= k (so the reasoning can be
checked by eye) and a verdict: "INDECOMPOSABLE" plus the counting argument, or, if some
combination could contribute, a note that the products must be built explicitly.  Real output
for the defaults (2 seconds):

    H_*(B_4 K_6) torsion: H_1: Z/2   H_2: Z/2
    75 combinations (|V(G1)|, |V(G2)|, k1, k2) with m1 + m2 <= 6 and k1 + k2 = 4
    combinations that could contribute 2-torsion in degree 2: 0
      K_5, k=2     H_1: Z/2
    ...
    No combination can contribute 2-torsion in degree 2, so the Z/2
    in H_2(B_4 K_6) is INDECOMPOSABLE: ...

USAGE
-----
    python3 experiments/kunneth_K6.py [--n 6] [--k 4]

    --n N   (default 6) the complete graph K_N whose degree-2 torsion is examined
    --k K   (default 4) the number of particles

Nothing is written to data/.  Exit code: always 0.

PITFALLS
--------
* Only complete pieces are computed.  Torsion is not monotone under taking subgraphs (K_9 has
  no Z/3 although it contains K_{4,5}), so a smaller piece could in principle carry torsion
  that K_m does not.  The conclusion rests on the vertex count, which covers every subgraph;
  the computed table is a check of that count, not a substitute for it.
* A piece on one vertex has no edges, which `gbt.fast` cannot handle; `homology_or_none` then
  returns None and that combination contributes nothing.  This is correct (a point carries
  no torsion), but such combinations are counted among the 75 without being computed.
* Hits are only possible contributions; a nonzero count would not prove decomposability.
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
    """H_*(B_k K_m) as {d: (free, invariant factors)}, or 'pt' for k = 0, or None.

    None means the computation failed, which happens exactly for m = 1 (a graph with no edges,
    which `gbt.fast` does not accept).  Results are memoised in `cache`, keyed by (m, k).
    """
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
    """(free, invariant factors) of degree d of a piece, or None if its homology is unknown.

    H is a homology dict, the marker "pt" for a factor with no particles (a point: Z in degree
    0, nothing else), or None.  Missing degrees are (0, []).
    """
    if H == "pt":
        return (1, []) if d == 0 else (0, [])
    if H is None:
        return None
    return H.get(d, (0, []))


def main():
    """Enumerate the splittings, list the possible Kunneth contributions, print the verdict."""
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
