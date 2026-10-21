#!/usr/bin/env python3
"""
An explicit top-degree chain X with  dX = (N - n + 1) y,  y primitive.

SETTING
-------
Gamma = K_{n+1,N} with parts A = {a_0 < ... < a_n} of size n+1 and B of size N, and n
particles.  A *matching configuration* is a top-degree basis element of the Swiatkowski
complex

        h_{v_1, t(v_1)} (x) ... (x) h_{v_n, t(v_n)},        v_1 < ... < v_n,

whose n edges {v, t(v)} are pairwise vertex disjoint.  Write S for the set of sources and
split it as S_A = S n A, S_B = S n B.  The A-vertices touched by the matching are
S_A u t(S_B); there are n of them, so since |A| = n+1 exactly one A-vertex is left over --
the "hole".  Then

    sigma(S,t) := [ S_A increasing | t(S_B) ordered by its B-SOURCE | the hole ]

lists every element of A exactly once, so it is a permutation of A, and we set

    X := sum over matching configurations of  (-1)^{|S_A|} sgn sigma(S,t) . (S,t).

CLAIM
-----
    (P1)  the "park" half of dX vanishes identically, term by term;
    (P2)  the "sit" half is (N - n + 1) y with every coefficient of y equal to +-1.

Hence dX = (N-n+1) y.  At N = 2n-1 this reads dX = n y, which together with the vanishing
theorem H_n(B_n K_{n+1,2n-1}) = 0 forces [y] to have order exactly n in H_{n-1}: Z/n
appears, giving Z/2 for K_{3,3} (the Ko--Park class), Z/3 for K_{4,5} and Z/4 for K_{5,7}.

WHY IT IS TRUE
--------------
X is a sum of products of flip chains  z(e) = h_{w,e} - h_{a,e},  one per edge of a matching,
antisymmetrised over the matchings with a given B-side; expanding the
products gives exactly the closed formula above.  Since  d z(e) = a - w  has no parked term,
(P1) is automatic; the terms with a particle sitting in B cancel in pairs because exchanging a
partner with the hole is a transposition; and the terms with a particle sitting in A do not
depend on which of the N-n+1 free B-vertices was used, which gives (P2).  The script
`wilson_block.py` checks the product form against the complex itself.

WHAT THIS SCRIPT DOES
---------------------
Builds X straight from the closed formula (no complex is built) and applies the differential
term by term.  The differential of a top-degree element (S, t) is

    d( h_{v_1} (x) ... (x) h_{v_n} ) = sum_i (-1)^(i-1) h_{v_1} (x) ... (e_i - v_i) ... h_{v_n},

so each source v_i contributes a "park" term (the particle parks on its edge e_i) and a "sit"
term (it sits at v_i), with the Koszul sign (-1)^(i-1) for the i-th source in increasing order.
The terms of dX are collected in a dictionary keyed by the resulting cell, and the script
reports |supp X|, |supp dX|, whether the park half vanished, the gcd g of dX, whether dX / g is
primitive, and whether g = N - n + 1.

It also runs the COMPETING reading of sigma, in which t(S_B) is listed in increasing order of
its A-value rather than of its B-source, so that the convention is pinned down by computation
rather than by recollection.  Under the wrong reading the park half does not vanish (18 terms
on K_{3,3}), so dX is not a multiple of a sit-only chain, and on K_{4,5} the gcd of dX even
collapses to 1.  The competing reading is tried on the cases with n <= 3.

WHAT IS PRINTED
---------------
For each case and convention a block

      K_{3,3}, n=2, sigma ordered by-source:
         |supp X| = 72,  |supp dX| = 36  (park 0, sit 36)
         coefficients of dX: [2],  gcd = 2
         park half vanishes identically  <-- (P1) holds
         dX = 2 . y   with y primitive: True,  y coefficients [1]
         predicted N-n+1 = 2   -> MATCH
         N = 2n-1, so this gives  Z/2  in H_1(B_2 K_{3,3})

("coefficients of dX" is the set of absolute values occurring; "y coefficients" likewise for
y = dX / gcd, which should be [1]).  The last line only appears when N = 2n-1.  At the end:

      All 9 cases satisfy (P1) and (P2) under the by-source convention.
      The by-value convention fails in every case tried:
      the ordering in the definition of sigma is forced.

USAGE
-----
    python3 experiments/torsion_class_X.py [n,N ...]

With no arguments it runs nine cases, (n, N) = (2,3), (2,4), (2,5), (2,6), (3,5), (3,6),
(3,7), (4,7), (4,8); the three with N = 2n-1 give Z/n.  Each argument is a pair n,N, e.g.

    python3 experiments/torsion_class_X.py 2,3 4,7

COST, OUTPUT FILES AND EXIT CODE
--------------------------------
The default run takes a few seconds (K_{5,8} with n = 4, with 134,400 matching
configurations, is the slow case) and little memory.  Writes no data file.  Exits with 1 if the
by-source convention fails in some case, or if the by-value convention does NOT fail in every
case tried; otherwise 0.

PITFALL
-------
The check "y primitive" is gcd(y) = 1; that every coefficient of y is +-1 is shown in the "y
coefficients" line but not enforced by the exit code.
"""

import itertools
import os
import sys
from collections import defaultdict
from math import gcd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def sgn_of_sequence(seq):
    """Sign (+1 or -1) of the permutation presented by `seq`, a sequence of distinct
    comparable items: (-1)^(number of inversions)."""
    s = 1
    a = list(seq)
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if a[i] > a[j]:
                s = -s
    return s


def matching_configurations(nA, nB, n):
    """All matching configurations of K_{nA,nB} with n particles.

    Vertices: A = 0..nA-1 and B = nA..nA+nB-1.  Yields (S, t) with S the sorted tuple of
    source vertices and t the dict source -> target.  Each of the n matching edges may be
    marked at either endpoint, which is the factor 2^n in the count.
    """
    A = list(range(nA))
    B = list(range(nA, nA + nB))
    for A0 in itertools.combinations(A, n):                 # A-vertices used
        for B0 in itertools.combinations(B, n):             # B-vertices used
            for perm in itertools.permutations(range(n)):   # the matching itself
                pairs = [(A0[i], B0[perm[i]]) for i in range(n)]
                for marks in itertools.product((0, 1), repeat=n):
                    # marks[i] = 0: the A-endpoint carries the mark; 1: the B-endpoint
                    t = {}
                    for (a, b), mk in zip(pairs, marks):
                        if mk == 0:
                            t[a] = b
                        else:
                            t[b] = a
                    yield tuple(sorted(t)), t


def coefficient(S, t, nA, convention):
    """The coefficient (-1)^{|S_A|} sgn sigma(S,t) of the matching configuration (S, t) in X.

    `convention` is "by-source" (the correct one: the targets of the B-sources are listed in
    increasing order of their B-source) or "by-value" (listed in increasing order of the
    A-vertex itself).  nA = n+1 is the size of A, so exactly one A-vertex (the hole) is
    unused and goes last."""
    SA = sorted(v for v in S if v < nA)
    SB = sorted(v for v in S if v >= nA)
    if convention == "by-source":
        tail = [t[v] for v in SB]                 # ordered by the B-source  (correct)
    else:
        tail = sorted(t[v] for v in SB)           # ordered by the A-value   (wrong)
    used = set(SA) | set(tail)
    hole = [v for v in range(nA) if v not in used]
    assert len(hole) == 1, "A must have exactly one spare vertex"
    seq = SA + tail + hole
    assert sorted(seq) == list(range(nA))
    return (-1) ** len(SA) * sgn_of_sequence(seq)


def boundary(S, t):
    """The Swiatkowski differential of the top-degree basis element (S,t).

        d(h_1 (x) ... (x) h_n) = sum_i (-1)^{i-1} h_1 (x)...(x) (e_i - v_i) (x)...(x) h_n

    A parked term is recorded as ('park', edge, rest, rt) and a sitting term as
    ('sit', v, rest, rt), where `rest` is the tuple of the other sources and `rt` the sorted
    pairs (source, target) of the other moving particles; these tuples identify the cell.
    Returns a list of (term, sign).
    """
    out = []
    order = list(S)
    for i, v in enumerate(order):
        sgn = (-1) ** i                 # Koszul sign of the i-th source (counting from 0)
        rest = tuple(x for x in order if x != v)
        rt = tuple(sorted((x, t[x]) for x in rest))
        edge = (min(v, t[v]), max(v, t[v]))
        out.append((("park", edge, rest, rt), sgn))
        out.append((("sit", v, rest, rt), -sgn))
    return out


def dX(n, N, convention):
    """Apply the differential to X for K_{n+1,N} with n particles.

    Returns (dict {boundary term: nonzero coefficient}, number of matching configurations,
    i.e. |supp X|).  Terms are those produced by boundary()."""
    nA = n + 1
    acc = defaultdict(int)
    count = 0
    for S, t in matching_configurations(nA, N, n):
        c = coefficient(S, t, nA, convention)
        count += 1
        for term, s in boundary(S, t):
            acc[term] += c * s
    return {k: v for k, v in acc.items() if v}, count


def report(n, N, convention):
    """Print the block described in the module docstring for one (n, N) and one convention.

    Returns True if (P1) and (P2) hold: no park terms, and dX = (N-n+1) y with y primitive."""
    d, nterms = dX(n, N, convention)
    park = {k: v for k, v in d.items() if k[0] == "park"}
    sit = {k: v for k, v in d.items() if k[0] == "sit"}
    g = 0
    for v in d.values():
        g = gcd(g, v)

    print("  K_{%d,%d}, n=%d, sigma ordered %s:" % (n + 1, N, n, convention))
    print("     |supp X| = %d,  |supp dX| = %d  (park %d, sit %d)"
          % (nterms, len(d), len(park), len(sit)))
    print("     coefficients of dX: %s,  gcd = %d"
          % (sorted({abs(v) for v in d.values()}), g))
    ok = True
    if park:
        print("     PARK HALF NONZERO on %d terms   <-- (P1) FAILS" % len(park))
        ok = False
    else:
        print("     park half vanishes identically  <-- (P1) holds")
    if g:
        prim = {k: v // g for k, v in d.items()}
        pg = 0
        for v in prim.values():
            pg = gcd(pg, v)
        print("     dX = %d . y   with y primitive: %s,  y coefficients %s"
              % (g, pg == 1, sorted({abs(v) for v in prim.values()})))
        print("     predicted N-n+1 = %d   -> %s"
              % (N - n + 1, "MATCH" if g == N - n + 1 else "MISMATCH"))
        if g != N - n + 1 or pg != 1:
            ok = False
    else:
        ok = False
    if N == 2 * n - 1 and ok:
        print("     N = 2n-1, so this gives  Z/%d  in H_%d(B_%d K_{%d,%d})"
              % (n, n - 1, n, n + 1, N))
    print(flush=True)
    return ok


def main():
    """Run the cases given on the command line (or the nine defaults) under both conventions.

    Returns 1 if the by-source convention fails anywhere or the by-value convention does not
    fail everywhere, else 0."""
    args = sys.argv[1:]
    if args:
        cases = [tuple(int(x) for x in a.split(",")) for a in args]
    else:
        # nine cases; the three with N = 2n-1 give Z/n
        cases = [(2, 3), (2, 4), (2, 5), (2, 6),
                 (3, 5), (3, 6), (3, 7),
                 (4, 7), (4, 8)]

    print("===== sigma's middle block ordered BY B-SOURCE (the correct convention) =====\n")
    good = all(report(n, N, "by-source") for (n, N) in cases)

    print("===== and now BY A-VALUE, the reading that does not work =====\n")
    bad_cases = [(n, N) for (n, N) in cases if n <= 3]       # enough to make the point
    broken = [(n, N) for (n, N) in bad_cases
              if not report(n, N, "by-value")]

    print("-" * 74)
    if not good:
        print("FAILURE: the by-source convention did not give dX = (N-n+1) y somewhere.")
        return 1
    print("All %d cases satisfy (P1) and (P2) under the by-source convention." % len(cases))
    if len(broken) == len(bad_cases):
        print("The by-value convention fails in every case tried:")
        print("the ordering in the definition of sigma is forced.")
    else:
        print("NOTE: the by-value convention did NOT fail everywhere; recheck the sign convention.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
