#!/usr/bin/env python3
"""
The torsion of H_{n-1}(B_n K_{n+1,b}) is the cokernel of an inclusion matrix.

THE STATEMENT
-------------
Let W = W_{t,k}(m) be the unsigned inclusion matrix of t-subsets against k-subsets of an
m-set: rows indexed by the t-subsets, columns by the k-subsets, and

        W[S][T] = 1  if S is a subset of T,   0 otherwise.

Take the unreduced Swiatkowski complex of K_{n+1,b} with n particles and b >= n (reducing at
the vertices of the b-side changes nothing essential; `wilson_block.py` checks the argument
step by step).  In suitable bases the top differential d_n is

        d_n  ~  [ I  0 ]
                [ 0  W_{n-1,n}(b) ]
                [ 0  0 ]

an identity block (the "flip" relations, which cancel) plus one genuinely interesting block,
which is the inclusion matrix.  Two consequences:

        H_n(B_n K_{n+1,b})             =  ker W_{n-1,n}(b)
        torsion H_{n-1}(B_n K_{n+1,b}) =  torsion coker W_{n-1,n}(b).

The first gives the top-degree vanishing for b <= 2n-1, via Gottlieb's rank theorem.  The
second is the sharp statement: it says exactly which torsion appears, and it is what produces
odd and non-squarefree torsion -- Z/3 in H_2(B_3 K_{4,5}) and Z/4 in H_3(B_4 K_{5,7}) --
rather than merely proving that some exists.

WHAT THIS SCRIPT CHECKS
-----------------------
    (1) Wilson's diagonal form of W_{t,k}(m): for t <= k and t + k <= m the matrix is
        equivalent to a diagonal one with entries C(k-i, t-i) repeated
        C(m,i) - C(m,i-1) times, i = 0..t  (R. M. Wilson, European J. Combin. 11 (1990)
        609-615).  Checked against an exact Smith normal form computed here.

        Note that Wilson's list is a DIAGONAL form, not a divisibility chain: for
        W_{3,5}(9) it reads 10, 6^8, 3^27, 1^48, and 6 does not divide 10.  The invariant
        that the two forms share is the multiset of elementary divisors -- the primary
        decomposition -- so that is what is compared.

    (2) Gottlieb's rank theorem: W_{t,k}(m) has full row rank C(m,t) when t <= k and
        t + k <= m (D. H. Gottlieb, Proc. Amer. Math. Soc. 17 (1966) 1233-1237).  Hence
        the numerical equivalence behind the vanishing,
        b <= 2n-1  <=>  C(b,n-1) >= C(b,n)  <=>  ker W_{n-1,n}(b) = 0.

    (3) The homological claim: coker W_{n-1,n}(b) computed from the matrix, against the
        torsion of H_{n-1}(B_n K_{n+1,b}) computed from the Swiatkowski complex, for every
        (n, b) small enough to build.

Check (3) is the one that could fail, and it is a real test of the cokernel statement: the two
sides are computed by completely different routes -- one a C(b,n-1) x C(b,n) combinatorial
matrix, the other a Morse reduction of a complex with tens of thousands of cells.

WHICH CASES
-----------
    (1), (2)  every W_{t,k}(m) with 2 <= m <= 12, 1 <= t < k, t + k <= m and k - t <= 2 (so
              that the exact Smith form stays affordable), and at most --matrix-cap entries
              C(m,t) * C(m,k).
    vanishing for 2 <= n <= 7 and n <= b <= 3n with C(b,n-1) * C(b,n) <= 800,000: the rank of
              W_{n-1,n}(b), exactly by Smith form below --matrix-cap entries and otherwise
              over F_1000003 (marked *; safe, since rank mod p never exceeds the rational
              rank and every diagonal entry of W_{n-1,n}(b) is at most n < p).
    (3)       2 <= n <= 4 and n+1 <= b <= 11, with at most --cap reduced cells and at most
              --matrix-cap matrix entries: 15 cases at the defaults, K_{3,3} to K_{4,9}.

WHAT IS PRINTED
---------------
Three tables.  (1)+(2): columns  W_{t,k}(m) | rows | cols | rank | full? | elementary divisors
(as "q^multiplicity") | vs Wilson ("matches", or "MISMATCH, predicted ..."; "RANK NOT FULL" is
appended if Gottlieb fails).  The vanishing table: columns  n | b | C(b,n-1) | C(b,n) |
rank W | dim ker | b <= 2n-1?.  (3): columns  K_{n+1,b} | n | cells | coker W torsion |
H_{n-1} torsion | ok/MISMATCH.  Then the number of cases of (3), and the verdict.  Excerpt
(the defaults, about 7 seconds):

    W_{t,k}(m)        rows    cols    rank full?  elementary divisors        vs Wilson
    W_{3,4}(7)          35      35      35 yes    2^14 3^6 4                 matches
    ...
    K_{4,5}           3     19,976 (Z/2)^4 + Z/3          (Z/2)^4 + Z/3          ok
    15 case(s) of the cokernel statement checked
    ALL CHECKS PASSED: Wilson's diagonal form, Gottlieb's rank theorem, and the
    cokernel statement hold on every case computed.

USAGE
-----
    python3 experiments/wilson_cokernel.py [--cap 150000] [--matrix-cap 20000]

    --cap N           skip a homological check (3) whose fully reduced complex has more than
                      N cells (default 150,000).
    --matrix-cap N    skip inclusion matrices with more than N entries (rows x columns) in the
                      exact checks (default 20,000).

COST, OUTPUT FILES AND EXIT CODE
--------------------------------
About 7 seconds at the defaults, little memory.  Writes nothing to data/.  Exits with 1 if any
check fails (a Wilson mismatch, a rank deficiency, a failure of the equivalence, or a cokernel
mismatch), otherwise 0.
"""

import argparse
import itertools
import os
import sys
from collections import Counter
from math import comb

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import complete_bip, primary, reduced_cell_count   # noqa: E402
from gbt import fast                                       # noqa: E402
from gbt.linalg import rank_mod_p, sparse_smith              # noqa: E402


def inclusion_columns(t, k, m):
    """W_{t,k}(m) by columns: (list of {row index: 1}, number of rows).

    Rows are the t-subsets of {0,...,m-1} in lexicographic order, columns the k-subsets (also
    lexicographic); column T has a 1 in the row of every t-subset of T.
    """
    rows = list(itertools.combinations(range(m), t))
    ri = {s: i for i, s in enumerate(rows)}
    cols = [{ri[S]: 1 for S in itertools.combinations(T, t)}
            for T in itertools.combinations(range(m), k)]
    return cols, len(rows)


def wilson_diagonal(t, k, m):
    """Wilson's predicted diagonal entries of W_{t,k}(m), as a Counter.

    Valid for t <= k and t + k <= m.  The i-th block contributes C(k-i, t-i) with
    multiplicity C(m,i) - C(m,i-1), where C(m,-1) = 0.
    """
    out = Counter()
    for i in range(t + 1):
        mult = comb(m, i) - (comb(m, i - 1) if i >= 1 else 0)
        if mult:
            out[comb(k - i, t - i)] += mult
    return out


def elementary_divisors(divs):
    """Multiset of prime powers occurring in a list of diagonal entries.

    This is the invariant shared by every diagonal form of a matrix, so it is what
    Wilson's list and the Smith normal form can be compared on.
    """
    out = Counter()
    for d in divs:
        if d in (0, 1):
            continue
        x = d
        p = 2
        while x > 1:
            # split off the full power of the least prime p dividing x
            while p * p <= x and x % p:
                p += 1
            if p * p > x:
                p = x
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            out[p ** e] += 1
            p = 2
    return out


def show(counter):
    """Format a Counter of prime powers as "2^14 3^6 4" (multiplicity 1 is not written)."""
    return " ".join("%d^%d" % (q, m) if m > 1 else "%d" % q
                    for q, m in sorted(counter.items())) or "-"


def main():
    """Run checks (1), (2), the vanishing equivalence and (3); returns 1 on any failure."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap", type=int, default=150_000,
                    help="skip the homological check above this many reduced cells")
    ap.add_argument("--matrix-cap", type=int, default=20_000,
                    help="skip inclusion matrices with more than this many entries")
    args = ap.parse_args()
    failures = []

    # ------------------------------------------------------------ (1) and (2)
    print("== (1) Wilson's diagonal form  and  (2) Gottlieb's rank theorem ==")
    print("   exact SNF of W_{t,k}(m), compared with the closed form on")
    print("   elementary divisors (see the docstring for why not on the raw list)\n")
    print("   restricted to t < k and t + k <= m, which is the regime of the theorem")
    print("   (W_{k,k} is the identity, and t + k > m is outside Wilson's hypothesis)\n")
    hdr = ("%-14s %7s %7s %7s %-6s %-26s %s"
           % ("W_{t,k}(m)", "rows", "cols", "rank", "full?", "elementary divisors",
              "vs Wilson"))
    print(hdr)
    print("-" * len(hdr))
    for m in range(2, 13):
        for k in range(2, m + 1):
            for t in range(1, k):
                # k - t <= 2 keeps the exact SNF affordable: the elementary divisors of
                # W_{t,k} grow like C(k,t), and so does the entry growth in the residual.
                if t + k > m or k - t > 2:
                    continue
                if comb(m, t) * comb(m, k) > args.matrix_cap:
                    continue
                cols, nr = inclusion_columns(t, k, m)
                rank, divs = sparse_smith(cols, verbose=False)
                got = elementary_divisors(divs)
                full = (rank == nr)

                pred = elementary_divisors(
                    [d for d, c in wilson_diagonal(t, k, m).items() for _ in range(c)])
                if pred == got:
                    verdict = "matches"
                else:
                    verdict = "MISMATCH, predicted %s" % show(pred)
                    failures.append("Wilson form W_{%d,%d}(%d)" % (t, k, m))
                if not full:
                    verdict += "; RANK NOT FULL"
                    failures.append("Gottlieb rank W_{%d,%d}(%d)" % (t, k, m))

                print("%-14s %7d %7d %7d %-6s %-26s %s"
                      % ("W_{%d,%d}(%d)" % (t, k, m), nr, len(cols), rank,
                         "yes" if full else "no", show(got), verdict), flush=True)
    print("-" * len(hdr))

    # ----------------------------------------- the equivalence behind the vanishing
    print("\n== the numerical fact behind the top-degree vanishing ==")
    print("   ker W_{n-1,n}(b) = 0   <=>   C(b,n-1) >= C(b,n)   <=>   b <= 2n-1\n")
    print("   rank marked * is computed over F_1000003 rather than by an exact SNF, which")
    print("   is safe here: rank mod p never exceeds the rational rank, and the diagonal")
    print("   entries of W_{n-1,n}(b) are all <= n (Wilson), so p cannot divide one.\n")
    print("   %3s %4s %10s %10s %9s %9s %12s"
          % ("n", "b", "C(b,n-1)", "C(b,n)", "rank W", "dim ker", "b <= 2n-1?"))
    for n in range(2, 8):
        for b in range(n, 3 * n + 1):
            if comb(b, n - 1) * comb(b, n) > 800_000:
                continue
            cols, nr = inclusion_columns(n - 1, n, b)
            if comb(b, n - 1) * comb(b, n) <= args.matrix_cap:
                rank, _ = sparse_smith(cols, verbose=False)
                mark = ""
            else:
                rank = rank_mod_p(cols, 1000003)
                mark = "*"
            kerdim = len(cols) - rank
            agree = ((kerdim == 0) == (b <= 2 * n - 1))
            if not agree:
                failures.append("equivalence at (n,b) = (%d,%d)" % (n, b))
            print("   %3d %4d %10d %10d %8d%s %9d %12s%s"
                  % (n, b, comb(b, n - 1), comb(b, n), rank, mark or " ", kerdim,
                     "yes" if b <= 2 * n - 1 else "no", "" if agree else "   MISMATCH"),
                  flush=True)

    # ------------------------------------------------------------------- (3)
    print("\n== (3) coker W_{n-1,n}(b)  vs  torsion H_{n-1}(B_n K_{n+1,b}) ==")
    print("   left: the inclusion matrix.  right: the Swiatkowski complex.\n")
    hdr = ("%-14s %4s %10s %-22s %-22s %s"
           % ("K_{n+1,b}", "n", "cells", "coker W torsion", "H_{n-1} torsion", ""))
    print(hdr)
    print("-" * len(hdr))
    tested = 0
    for n in range(2, 5):
        for b in range(n + 1, 12):
            G = complete_bip(n + 1, b)
            sz = reduced_cell_count(G, n)
            if sz > args.cap or comb(b, n - 1) * comb(b, n) > args.matrix_cap:
                continue
            cols, _ = inclusion_columns(n - 1, n, b)
            _, divs = sparse_smith(cols, verbose=False)

            H, _, _ = fast.compute(G, n)
            hom_tors = H.get(n - 1, (0, []))[1]

            # compare primary decompositions: equal exactly when the groups are isomorphic
            same = (primary(divs) == primary(hom_tors))
            tested += 1
            if not same:
                failures.append("cokernel statement at K_{%d,%d}, n=%d" % (n + 1, b, n))
            print("%-14s %4d %10s %-22s %-22s %s"
                  % ("K_{%d,%d}" % (n + 1, b), n, format(sz, ","),
                     primary(divs), primary(hom_tors),
                     "ok" if same else "MISMATCH"), flush=True)
    print("-" * len(hdr))
    print("%d case(s) of the cokernel statement checked" % tested)

    print()
    if failures:
        print("%d FAILURE(S): %s" % (len(failures), failures))
        return 1
    print("ALL CHECKS PASSED: Wilson's diagonal form, Gottlieb's rank theorem, and the")
    print("cokernel statement hold on every case computed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
