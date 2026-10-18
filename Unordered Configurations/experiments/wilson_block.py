#!/usr/bin/env python3
"""
The top differential of K_{n+1,b}, checked step by step: why W sits next to an identity block.

THE ARGUMENT BEING CHECKED
--------------------------
Gamma = K_{n+1,b}, n particles, the UNREDUCED Swiatkowski complex.  The rows of the top
differential d_n come in three kinds:

    parked      n-1 moving particles and one parked on an edge
    B-sitting   n-1 moving particles and one sitting at a vertex of B
    A-sitting   n-1 moving particles and one sitting at a vertex of A

Write d' for the parked and B-sitting rows and d'' for the A-sitting rows.

    (1) Flip chains.  For every n-subset B_0 of B,
            d X_{B_0} = (-1)^{n-1} sum_{T subset B_0} Y_T,
        with X_{B_0}, Y_T the explicit chains of gbt.flipchains.  In particular d' X_{B_0} = 0.
    (2) X_{B_0} is the flip-chain class X restricted to B_0: its coefficients are
        (-1)^{|S_A|} sgn sigma(S,t) (see gbt.flipchains.closed_formula_X).
    (3) The identity block.  d' has every invariant factor equal to 1, and
        its kernel has rank C(b,n) -- so it is spanned by the X_{B_0}.
    (4) The Y_T have pairwise disjoint supports and all their entries are +-1, so they extend
        to a basis of the A-sitting cells.
    (5) The conclusion.  The non-unit invariant factors of d_n are those of W_{n-1,n}(b).
    (6) The completion F of the row basis can be taken to be standard basis vectors --
        unit pivots in the rows R' clear every non-reference column.  (Not needed for the
        argument; this is a computational observation.)

Every item is computed from the complex itself, not from the formulas being tested.

    python3 experiments/wilson_block.py [n,b ...]

Default cases: (n,b) = (2,3), (2,4), (2,5), (3,4), (3,5), (3,6) -- about a minute in total.
Exits non-zero if anything fails.
"""

import itertools
import os
import sys
from math import comb

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import complete_bip, swiatkowski                                 # noqa: E402
from gbt.flipchains import X_B0, Y_T, add_into, closed_formula_X          # noqa: E402
from gbt.linalg import sparse_smith                                       # noqa: E402


def row_kind(r, a):
    """'P' (parked), 'VA' (sitting in A) or 'VB' (sitting in B) for an (n-1)-cell r."""
    mono, st = r
    if mono:
        return "P"
    v = st.index(1)
    return "VA" if v < a else "VB"


def nonunit(cols):
    """Rank and non-unit invariant factors of a matrix given by {row: coeff} columns."""
    rows = sorted({r for c in cols for r in c}, key=repr)
    ri = {r: i for i, r in enumerate(rows)}
    rk, divs = sparse_smith([{ri[r]: x for r, x in c.items()} for c in cols], verbose=False)
    return rk, sorted(d for d in divs if d not in (0, 1))


def unit_pivot(cols, allowed):
    """Greedy elimination with +-1 pivots in the allowed rows; returns #columns not pivoted.

    Each step picks a column with a +-1 entry in an allowed row r, and subtracts multiples of
    it from the other columns to clear row r.  Only column operations among the given columns
    are used, so if every column gets pivoted the pivot rows form a square submatrix of
    determinant +-1.
    """
    cols = {i: dict(c) for i, c in enumerate(cols)}
    byrow = {}
    for i, c in cols.items():
        for r in c:
            byrow.setdefault(r, set()).add(i)
    progress = True
    while progress:
        progress = False
        for r in list(byrow):
            if r not in allowed:
                continue
            cand = [i for i in byrow[r] if abs(cols[i].get(r, 0)) == 1]
            if not cand:
                continue
            p = min(cand, key=lambda i: len(cols[i]))
            pc = cols.pop(p)
            for rr in pc:
                byrow[rr].discard(p)
            for i in list(byrow[r]):
                f = cols[i][r] * pc[r]
                for rr, x in pc.items():
                    y = cols[i].get(rr, 0) - f * x
                    if y:
                        if rr not in cols[i]:
                            byrow.setdefault(rr, set()).add(i)
                        cols[i][rr] = y
                    elif rr in cols[i]:
                        del cols[i][rr]
                        byrow[rr].discard(i)
            del byrow[r]
            progress = True
    return len(cols)


def wilson(n, b):
    """Columns of W_{n-1,n}(b): n-subsets -> the (n-1)-subsets they contain."""
    return [{T: 1 for T in itertools.combinations(S, n - 1)}
            for S in itertools.combinations(range(b), n)]


def check(n, b):
    a = n + 1
    G = complete_bip(a, b)
    D = swiatkowski(G, n).B[n]                 # columns: n-cells -> {(n-1)-cell: coeff}
    print("K_{%d,%d}, n=%d: %d n-cells" % (a, b, n, len(D)), flush=True)
    ok = True

    # (1) and (2): flip chains
    sigma_formula = closed_formula_X(n, b)
    good1 = good2 = True
    for B0 in itertools.combinations(range(a, a + b), n):
        x = X_B0(n, b, B0)
        good2 &= all(sigma_formula[c] == v for c, v in x.items())
        bd = {}
        for c, v in x.items():
            for r, y in D[c].items():
                bd[r] = bd.get(r, 0) + v * y
        bd = {r: v for r, v in bd.items() if v}
        pred = {}
        for T in itertools.combinations(B0, n - 1):
            add_into(pred, Y_T(n, b, T), (-1) ** (n - 1))
        good1 &= (bd == pred)
    print("   (1) dX_{B_0} = (-1)^(n-1) sum_T Y_T for every B_0:          %s" % good1)
    print("   (2) X_{B_0} has the coefficients (-1)^|S_A| sgn sigma:       %s" % good2)
    ok &= good1 and good2

    # (3) the identity block
    dprime = [{r: x for r, x in col.items() if row_kind(r, a) != "VA"} for col in D.values()]
    rk, nu = nonunit(dprime)
    good3 = (not nu) and (len(D) - rk == comb(b, n))
    print("   (3) d' = parked + B-sitting rows: rank %d, non-unit factors %s, kernel rank %d"
          " = C(b,n) = %d:  %s" % (rk, nu or "none", len(D) - rk, comb(b, n), good3))
    ok &= good3

    # (4) the Y_T
    supports = {}
    entries = set()
    for T in itertools.combinations(range(a, a + b), n - 1):
        y = Y_T(n, b, T)
        assert all(row_kind(r, a) == "VA" for r in y)
        supports[T] = set(y)
        entries |= {abs(v) for v in y.values()}
    allsupp = [r for s in supports.values() for r in s]
    good4 = (len(allsupp) == len(set(allsupp))) and entries == {1}
    print("   (4) the %d chains Y_T: disjoint supports, entries +-1:          %s"
          % (len(supports), good4))
    ok &= good4

    # (6) the completion F of the row basis can be taken standard, for the
    #     reference configurations c_0(B_0) = h_{a_1,{a_1,w_1}} ... h_{a_n,{a_n,w_n}}
    refs = set()
    for B0 in itertools.combinations(range(a, a + b), n):
        st = [0] * G.N
        for i, w in enumerate(B0):
            st[i + 1] = 2 + G.edges.index((i + 1, w))
        refs.add(((), tuple(st)))
    rprime = {r for col in D.values() for r in col if row_kind(r, a) != "VA"}
    left = unit_pivot([{r: x for r, x in D[c].items() if r in rprime}
                       for c in D if c not in refs], rprime)
    good6 = (left == 0)
    print("   (6) unit pivots in the rows R' clear all %d non-reference columns:  %s"
          % (len(D) - len(refs), good6))
    ok &= good6

    # (5) the conclusion
    _, nu_d = nonunit(list(D.values()))
    _, nu_w = nonunit(wilson(n, b))
    good5 = (nu_d == nu_w)
    print("   (5) non-unit factors of d_n %s = those of W_{%d,%d}(%d) %s:  %s"
          % (nu_d or "none", n - 1, n, b, nu_w or "none", good5))
    ok &= good5
    print(flush=True)
    return ok


def main():
    args = sys.argv[1:]
    cases = ([tuple(int(x) for x in s.split(",")) for s in args] if args else
             [(2, 3), (2, 4), (2, 5), (3, 4), (3, 5), (3, 6)])
    results = [check(n, b) for n, b in cases]
    if all(results):
        print("ALL CHECKS PASSED: in every case d_n ~ W (+) I (+) 0, for the reasons above.")
        return 0
    print("SOME CHECK FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
