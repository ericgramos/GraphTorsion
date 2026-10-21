#!/usr/bin/env python3
"""
Does the implementation compute what it claims to?  Five independent checks.

    (1) d o d = 0 on the Swiatkowski complex, reduced and unreduced, on a spread of graphs.
    (2) The closed-form cell counts of `gbt.sizes` agree with the complexes actually built.
    (3) Partial reduction does not change homology: reducing at no vertices, at one side
        of a bipartite graph, and at every vertex all give the same answer.
    (4) The Swiatkowski model agrees with Abrams' discretised model on cases small enough
        for the latter, and both reproduce the values known from the literature:
        B_2(K_5) and B_2(K_{3,3}) are the closed non-orientable surfaces N_7 and N_5
        (Abrams, PhD thesis 2000), so H_1 = Z^6 + Z/2 and Z^4 + Z/2.
    (5) "Torsion is a cokernel": for every degree d,
            torsion of coker(d_{d+1})  ==  torsion of H_d,
        while the free ranks differ by exactly rank(d_d).  This is the lemma that lets one
        read the torsion of H_{n-1} off a single differential without ever computing
        ker(d_{n-1}) -- see the discussion below the check.

Two smaller groups ride along: (3b) `gbt.fast` and `gbt.lean` give identical homology, and
under (4) the path P_4 and the cycle C_4 have the homology of a point and of a circle
(unordered configurations on an interval form a contractible space, and on a circle a space
homotopy equivalent to a circle), so nothing exotic may appear there.

Nothing here depends on any of the torsion results; it is purely a test of the machinery.
Run it first if anything looks surprising.

WHAT EACH CHECK COMPARES
------------------------
    (1)  `ChainComplex.check_d2()` on the unreduced and fully reduced complexes of K_5 (k=2),
         K_{3,3} (k=3), the octahedron (k=3) and theta(3,2) (k=4), and on K_{4,5} (k=3)
         reduced at the B side only.
    (2)  the number of cells actually built (`fast.build` for the reduced complex,
         `swiatkowski` for the unreduced one) against `reduced_cell_count` / `cell_count`,
         for K_6 (k=4), K_{4,5} (k=3), the octahedron (k=4) and the triangular prism (k=3).
    (3)  the full homology of K_{3,3}, K_{3,4} (k=2) and K_{4,5} (k=3) computed three ways --
         reduce_at none, B only, all -- with the reference pipeline, and against `gbt.fast`.
    (3b) `fast` against `lean` on K_{3,3} (k=3), K_6 (k=4), K_{4,5} (k=3).
    (4)  Abrams' model on the subdivided graph against the Swiatkowski complex: K_5 at k=2
         (no subdivision needed) and the octahedron at k=3 (edges subdivided twice), plus
         d o d = 0 for Abrams' complex; then the literature values for B_2 K_5 and B_2 K_{3,3},
         and the path P_4 and cycle C_4 at k = 2, 3.
    (5)  (a) on the fully reduced complexes of K_5 and K_{3,3} at k = 2, BEFORE any Morse
         reduction, and (b) on the Morse-reduced complexes of K_{4,5} (k=3), K_6 (k=4) and the
         octahedron (k=4): for every degree d, the non-unit invariant factors of d_{d+1}
         (computed with `sparse_smith`) equal the torsion of H_d, and
         (#d-cells - rank d_{d+1}) - (free rank of H_d) = rank d_d.

WHAT IS PRINTED
---------------
A section header per group, then one line per check:

    [ok  ] NAME: VALUE                      the check passed; VALUE is what was computed
    [FAIL] NAME: VALUE   (expected WANT)    the check failed

and finally "ALL CHECKS PASSED" (67 checks) or "N CHECK(S) FAILED: [names]".  Excerpt:

    == (2) closed-form cell counts ==
      [ok  ] reduced count, K_6, k=4: 71220
      [ok  ] unreduced count, K_6, k=4: 176580
    == (3) the amount of reduction does not matter ==
      [ok  ] K_{4,5}, k=3: reduce_at none == B only: {0: (1, ()), 1: (12, (2,)), ...}

USAGE
-----
    python3 experiments/model_sanity.py          (about 15 seconds here; well under 1 GB)

No options; writes nothing.  Exit code 1 if any check fails, 0 otherwise.  The checks run at
import time (the script has no main function), so do not import this file from another script.
"""

import os
import sys
import time

# make the repository root importable no matter where this is run from
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import (Graph, abrams, complete, complete_bip, octahedron, prism,        # noqa: E402
                 theta, swiatkowski, subdivide, cell_count, reduced_cell_count)
from gbt import fast, lean                                                        # noqa: E402
from gbt.linalg import indexed_columns, sparse_smith                              # noqa: E402

FAILURES = []


def check(name, got, want):
    """Print one "[ok  ]" or "[FAIL]" line comparing got with want; record failures.

    Comparison is plain ==, so both sides must have the same shape (dicts of tuples, ints...).
    Returns True if they agree.
    """
    ok = (got == want)
    print("  [%s] %s: %s%s" % ("ok  " if ok else "FAIL", name, got,
                               "" if ok else "   (expected %s)" % (want,)))
    if not ok:
        FAILURES.append(name)
    return ok


def hom(G, k, pipeline="fast"):
    """{degree: (free rank, tuple of invariant factors)} via the named pipeline."""
    H, _, _ = (fast if pipeline == "fast" else lean).compute(G, k)
    return {d: (H[d][0], tuple(H[d][1])) for d in H}


# ---------------------------------------------------------------------------
print("== (1) d o d = 0 ==")
for name, G, k in [("K_5, k=2", complete(5), 2),
                   ("K_{3,3}, k=3", complete_bip(3, 3), 3),
                   ("octahedron, k=3", octahedron(), 3),
                   ("theta(3,2), k=4", theta(3, 2), 4)]:
    for tag, R in [("unreduced", ()),
                   ("fully reduced", tuple(range(G.N)))]:
        C = swiatkowski(G, k, reduce_at=R)
        check("%s [%s]" % (name, tag), C.check_d2(), True)
# partial reduction at one side of a bipartite graph
Gb = complete_bip(4, 5)
C = swiatkowski(Gb, 3, reduce_at=tuple(range(4, 9)))
check("K_{4,5}, k=3 [reduced at B only]", C.check_d2(), True)

# ---------------------------------------------------------------------------
print("\n== (2) closed-form cell counts ==")
for name, G, k in [("K_6, k=4", complete(6), 4),
                   ("K_{4,5}, k=3", complete_bip(4, 5), 3),
                   ("octahedron, k=4", octahedron(), 4),
                   ("prism, k=3", prism(), 3)]:
    built_red = sum(fast.build(G, k)[0].values())
    check("reduced count, %s" % name, built_red, reduced_cell_count(G, k))
    built_raw = sum(swiatkowski(G, k).sizes().values())
    check("unreduced count, %s" % name, built_raw, cell_count(G, k))

# ---------------------------------------------------------------------------
print("\n== (3) the amount of reduction does not matter ==")
for name, G, k, halves in [("K_{3,3}, k=2", complete_bip(3, 3), 2, (3, 3)),
                           ("K_{3,4}, k=2", complete_bip(3, 4), 2, (3, 4)),
                           ("K_{4,5}, k=3", complete_bip(4, 5), 3, (4, 5))]:
    a, b = halves
    answers = {}
    for tag, R in [("none", ()),
                   ("B only", tuple(range(a, a + b))),
                   ("all", tuple(range(a + b)))]:
        C = swiatkowski(G, k, reduce_at=R)
        C.reduce_fully()
        answers[tag] = {d: (v[0], tuple(v[1])) for d, v in C.homology().items()}
    check("%s: reduce_at none == B only" % name, answers["none"], answers["B only"])
    check("%s: reduce_at none == all" % name, answers["none"], answers["all"])
    check("%s: reduce_at all == fast pipeline" % name, answers["all"], hom(G, k))

print("\n== (3b) the two fast pipelines agree ==")
for name, G, k in [("K_{3,3}, k=3", complete_bip(3, 3), 3),
                   ("K_6, k=4", complete(6), 4),
                   ("K_{4,5}, k=3", complete_bip(4, 5), 3)]:
    check("fast == lean, %s" % name, hom(G, k, "fast"), hom(G, k, "lean"))

# ---------------------------------------------------------------------------
print("\n== (4) Abrams' model, and values known from the literature ==")
# Abrams requires sufficient subdivision; t is chosen by hand per case and the
# resulting cubical complex is huge, so only the two smallest cases are feasible.
for name, G, k, t in [("K_5, k=2", complete(5), 2, 1),
                      ("octahedron, k=3", octahedron(), 3, 2)]:
    t0 = time.time()
    S = swiatkowski(G, k, reduce_at=tuple(range(G.N)))
    S.reduce_fully()
    A = abrams(subdivide(G, t), k)
    check("Abrams d o d = 0, %s" % name, A.check_d2(), True)
    A.reduce_fully()
    check("Abrams == Swiatkowski, %s (subdivided %dx, %.0fs)" % (name, t, time.time() - t0),
          S.homology(), A.homology())

print("\n  B_2(K_5) is the surface N_7 and B_2(K_{3,3}) is N_5 (Abrams 2000):")
check("H_1(B_2 K_5)", hom(complete(5), 2)[1], (6, (2,)))
check("H_2(B_2 K_5)", hom(complete(5), 2)[2], (0, ()))
check("H_1(B_2 K_{3,3})", hom(complete_bip(3, 3), 2)[1], (4, (2,)))
check("H_2(B_2 K_{3,3})", hom(complete_bip(3, 3), 2)[2], (0, ()))
print("  (the Z/2 is the Ko--Park class; H_2 vanishes because these are the")
print("   NON-orientable closed surfaces, whose top integral homology is Z/2 in H_1)")

print("\n  trees and cycles: configuration spaces are homotopy equivalent to a point")
print("  or a circle, so nothing exotic may appear --")
PATH = Graph(4, [(0, 1), (1, 2), (2, 3)])
CYCLE = Graph(4, [(0, 1), (1, 2), (2, 3), (0, 3)])
for k in (2, 3):
    H = hom(PATH, k)
    check("path P_4, k=%d" % k, (H[0], H.get(1, (0, ()))), ((1, ()), (0, ())))
    H = hom(CYCLE, k)
    check("cycle C_4, k=%d" % k, (H[0], H.get(1, (0, ()))), ((1, ()), (1, ())))

# ---------------------------------------------------------------------------
print("\n== (5) torsion is a cokernel ==")
print("  coker(d_{d+1}) = C_d / im d_{d+1} forgets to intersect with ker d_d, so it is")
print("  bigger than H_d -- but only by the free group im d_d.  Hence the torsion agrees")
print("  and the free ranks differ by exactly rank d_d.\n")

print("  (a) on the fully reduced complex, before any discrete Morse reduction:")
for name, G, k in [("K_5, k=2", complete(5), 2),
                   ("K_{3,3}, k=2", complete_bip(3, 3), 2)]:
    C = swiatkowski(G, k, reduce_at=tuple(range(G.N)))
    H = C.homology()
    for d in sorted(C.gens):
        if d == 0:
            continue
        cols_up, _, _ = indexed_columns(C, d + 1)
        cols_here, _, _ = indexed_columns(C, d)
        rank_up, divs_up = sparse_smith(cols_up, verbose=False) if cols_up else (0, [])
        rank_here, _ = sparse_smith(cols_here, verbose=False) if cols_here else (0, [])
        check("%s, degree %d: tors coker d_%d == tors H_%d" % (name, d, d + 1, d),
              tuple(sorted(x for x in divs_up if x not in (0, 1))),
              tuple(sorted(H[d][1])))
        check("%s, degree %d: free ranks differ by rank d_%d = %d"
              % (name, d, d, rank_here),
              (len(C.gens[d]) - rank_up) - H[d][0], rank_here)

print("\n  (b) on Morse-reduced complexes, which reaches the interesting cases:")
for name, G, k in [("K_{4,5}, k=3", complete_bip(4, 5), 3),
                   ("K_6, k=4", complete(6), 4),
                   ("octahedron, k=4", octahedron(), 4)]:
    sizes, B = fast.morse(*fast.build(G, k))
    H = fast.homology(sizes, B)
    for d in sorted(sizes):
        if d == 0 or not sizes[d]:
            continue
        up = B.get(d + 1) or []
        here = B.get(d) or []
        rank_up, divs_up = sparse_smith(up, verbose=False) if up else (0, [])
        rank_here, _ = sparse_smith(here, verbose=False) if here else (0, [])
        check("%s, degree %d: tors coker d_%d == tors H_%d" % (name, d, d + 1, d),
              tuple(sorted(x for x in divs_up if x not in (0, 1))),
              tuple(sorted(H[d][1])))
        check("%s, degree %d: free ranks differ by rank d_%d = %d"
              % (name, d, d, rank_here),
              (sizes[d] - rank_up) - H[d][0], rank_here)
print()
print("  This is why one can read the torsion of H_{n-1} off the single matrix d_n:")
print("  ker d_{n-1} changes the free rank of H_{n-1} and nothing else.")

print()
if FAILURES:
    print("%d CHECK(S) FAILED: %s" % (len(FAILURES), FAILURES))
    sys.exit(1)
print("ALL CHECKS PASSED")
