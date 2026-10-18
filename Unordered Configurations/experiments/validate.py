#!/usr/bin/env python3
"""
Regression suite: every headline number that is cheap enough to recompute.

Run this after changing anything in `gbt/`.  It exits non-zero on the first failure, and each
line names the claim it is checking, so a failure points at the statement it would invalidate.

    python3 experiments/validate.py            (about five minutes)

The heavy counterparts -- the million-cell cases -- are in `validate_slow.py`.

What is checked, in order:

    1.  Values known independently from the literature (Abrams' surfaces, Ko--Park).
    2.  The top-degree vanishing for a = n+1, b <= 2n-1.
    3.  The exact torsion of H_{n-1}(B_n K_{n+1,b}), which is the
        cokernel of the Wilson inclusion matrix, including the Z/3 of K_{4,5}.
    4.  The dimension formula dim H_n = max(0, C(b,n) - C(b,n-1)).
    5.  Complete graphs: 2-torsion only, with multiplicity rank H_1 for n >= 8 at k = 3.
    6.  Planar graphs with torsion: the octahedron and the cube.
    7.  The three An--Knudsen planar atoms are torsion-free.
    8.  Internal consistency: d^2 = 0; the closed-form cell counts; reduce_at makes no
        difference; `fast` and `lean` agree; Abrams' model agrees with Swiatkowski.
"""

import os
import sys
from math import comb

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import (Graph, complete, complete_bip, cube, octahedron, prism,   # noqa: E402
                 theta, swiatkowski, subdivide, abrams,
                 cell_count, reduced_cell_count, primary)
from gbt import fast, lean                                                # noqa: E402

FAILURES = []


def check(name, got, want):
    ok = (got == want)
    print("  [%s] %s: %s%s" % ("ok  " if ok else "FAIL", name, got,
                               "" if ok else "   (expected %s)" % (want,)))
    if not ok:
        FAILURES.append(name)


def hom(G, k, pipeline="fast"):
    H, _, _ = (fast if pipeline == "fast" else lean).compute(G, k)
    return {d: (H[d][0], tuple(H[d][1])) for d in H}


STAR = Graph(5, [(0, 1), (0, 2), (0, 3), (0, 4)])
CYCLE4 = Graph(4, [(0, 1), (1, 2), (2, 3), (0, 3)])

# ---------------------------------------------------------------------------
print("== 1. values known from the literature ==")
# Abrams: B_2(K_5) = N_7 and B_2(K_{3,3}) = N_5, closed non-orientable surfaces
check("H_1(B_2 K_5)", hom(complete(5), 2)[1], (6, (2,)))
check("H_2(B_2 K_5)", hom(complete(5), 2)[2], (0, ()))
check("H_1(B_2 K_{3,3})", hom(complete_bip(3, 3), 2)[1], (4, (2,)))
# Ko--Park: H_1 has torsion exactly for non-planar graphs, and then a single Z/2
check("H_1(B_2 K_4) torsion-free (planar)", hom(complete(4), 2)[1][1], ())
check("H_1(B_3 K_5) = Z/2 (non-planar)", hom(complete(5), 3)[1][1], (2,))
check("H_1(B_4 octahedron) torsion-free (planar)", hom(octahedron(), 4)[1][1], ())

print("\n== 2. top-degree vanishing: H_n(B_n K_{n+1,b}) = 0 for b <= 2n-1 ==")
for (a, b, n) in [(3, 3, 2), (4, 4, 3), (4, 5, 3)]:
    check("dim H_%d(B_%d K_%d,%d)  [b <= 2n-1]" % (n, n, a, b),
          hom(complete_bip(a, b), n)[n][0], 0)
print("  and a <= n always vanishes:")
for (a, b, n) in [(2, 5, 2), (3, 7, 3), (2, 8, 4)]:
    check("dim H_%d(B_%d K_%d,%d)  [a <= n]" % (n, n, a, b),
          hom(complete_bip(a, b), n).get(n, (0, ()))[0], 0)

print("\n== 3. the torsion is coker of the Wilson matrix ==")
check("K_{3,3} n=2, H_1", hom(complete_bip(3, 3), 2)[1], (4, (2,)))
check("K_{4,4} n=3, H_2", hom(complete_bip(4, 4), 3)[2], (80, (2,)))
check("K_{4,5} n=3, H_2  [the first odd torsion]",
      hom(complete_bip(4, 5), 3)[2], (155, (2, 2, 2, 6)))
check("K_{4,5} n=3, H_2 primary", primary(hom(complete_bip(4, 5), 3)[2][1]),
      "(Z/2)^4 + Z/3")
check("K_{4,6} n=3, H_2", hom(complete_bip(4, 6), 3)[2], (259, (2, 2, 2, 2, 6)))
check("K_{4,7} n=3, H_2 primary", primary(hom(complete_bip(4, 7), 3)[2][1]),
      "(Z/2)^6 + Z/3")

print("\n== 4. dim H_n = max(0, C(b,n) - C(b,n-1)) for a = n+1 ==")
for (b, n) in [(3, 2), (4, 2), (5, 2), (6, 2), (7, 2), (4, 3), (5, 3), (6, 3), (7, 3)]:
    want = max(0, comb(b, n) - comb(b, n - 1))
    check("dim H_%d(B_%d K_%d,%d)" % (n, n, n + 1, b),
          hom(complete_bip(n + 1, b), n)[n][0], want)

print("\n== 5. complete graphs: 2-torsion only ==")
check("K_7 n=3, H_2", hom(complete(7), 3)[2], (350, (2,)))
check("K_8 n=3, H_2", hom(complete(8), 3)[2], (734, tuple([2] * 21)))
check("K_9 n=3, H_2", hom(complete(9), 3)[2], (1401, tuple([2] * 28)))
check("K_6 n=4, H_2", hom(complete(6), 4)[2], (264, (2,)))
print("  multiplicity vs rank H_1 = C(n-1,2), for n >= 8 at k = 3:")
for n in (8, 9):
    H = hom(complete(n), 3)
    check("K_%d: m((Z/2)) == C(%d,2)" % (n, n - 1),
          (sum(1 for x in H[2][1] if x % 2 == 0), H[1][0]),
          (comb(n - 1, 2), comb(n - 1, 2)))

print("\n== 6. planar graphs can carry 2-torsion ==")
check("octahedron n=4, H_2", hom(octahedron(), 4)[2], (95, (2,)))
check("cube Q_3 n=4, H_2", hom(cube(), 4)[2], (29, (2,)))
check("octahedron n=3, H_2 (torsion-free)", hom(octahedron(), 3)[2], (47, ()))
check("prism n=4, H_2 (control, torsion-free)", hom(prism(), 4)[2][1], ())

print("\n== 7. the three An--Knudsen planar atoms are torsion-free ==")
for name, G in [("cycle C_4", CYCLE4), ("star K_{1,4}", STAR), ("theta", theta(3, 2))]:
    for k in (3, 4):
        H = hom(G, k)
        check("%s, k=%d: no torsion anywhere" % (name, k),
              any(H[d][1] for d in H), False)

print("\n== 8. internal consistency ==")
for name, G, k in [("K_5 k=2", complete(5), 2), ("K_{3,3} k=3", complete_bip(3, 3), 3),
                   ("octahedron k=3", octahedron(), 3)]:
    for tag, R in [("unreduced", ()), ("reduced", tuple(range(G.N)))]:
        check("d^2 = 0, %s [%s]" % (name, tag),
              swiatkowski(G, k, reduce_at=R).check_d2(), True)
for name, G, k in [("K_6 k=4", complete(6), 4), ("K_{4,5} k=3", complete_bip(4, 5), 3)]:
    check("reduced cell count, %s" % name,
          sum(fast.build(G, k)[0].values()), reduced_cell_count(G, k))
    check("unreduced cell count, %s" % name,
          sum(swiatkowski(G, k).sizes().values()), cell_count(G, k))
G = complete_bip(4, 5)
ans = {}
for tag, R in [("none", ()), ("B only", tuple(range(4, 9))), ("all", tuple(range(9)))]:
    C = swiatkowski(G, 3, reduce_at=R)
    C.reduce_fully()
    ans[tag] = {d: (v[0], tuple(v[1])) for d, v in C.homology().items()}
check("K_{4,5} k=3: reduce_at none == B only", ans["none"], ans["B only"])
check("K_{4,5} k=3: reduce_at none == all", ans["none"], ans["all"])
check("K_{4,5} k=3: reduce_at all == fast", ans["all"], hom(G, 3))
for name, G, k in [("K_{3,3} k=3", complete_bip(3, 3), 3), ("K_6 k=4", complete(6), 4),
                   ("K_{4,5} k=3", complete_bip(4, 5), 3)]:
    check("fast == lean, %s" % name, hom(G, k, "fast"), hom(G, k, "lean"))
for name, G, k, t in [("K_5 k=2", complete(5), 2, 1),
                      ("octahedron k=3", octahedron(), 3, 2)]:
    S = swiatkowski(G, k, reduce_at=tuple(range(G.N)))
    S.reduce_fully()
    A = abrams(subdivide(G, t), k)
    A.reduce_fully()
    check("Abrams == Swiatkowski, %s" % name, S.homology(), A.homology())

print()
if FAILURES:
    print("%d CHECK(S) FAILED: %s" % (len(FAILURES), FAILURES))
    sys.exit(1)
print("ALL CHECKS PASSED")
