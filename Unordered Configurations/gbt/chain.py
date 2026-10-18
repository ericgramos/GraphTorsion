"""
Chain complexes of free abelian groups: d^2 test, discrete Morse reduction, homology.

This is the reference pipeline.  It keeps generators as opaque hashable labels (for the
Swiatkowski complex, the (mono, st) pairs), which makes it easy to inspect and to debug
but costs a great deal of Python object overhead.  Use it when you want to see what the
generators are; use `gbt.fast` or `gbt.lean` when you want answers for large complexes.

THE REDUCTION
-------------
`reduce_fully` is algebraic discrete Morse theory (Forman).  If <d c, r> = +-1 then the
pair (c, r) can be cancelled by a unimodular change of basis, so the invariant factors of
the reduced complex are those of the original with some 1's removed; homology is
unchanged.  Pivots are taken in order of increasing Markowitz fill estimate
(|column| - 1) * (|row| - 1), which keeps the residual matrices sparse.

Cancelling the pair (c in degree d, r in degree d-1) requires THREE updates, and the
third is the one that is easy to forget:

    * delete the column c from B[d];
    * clear the row r from the other columns of B[d] (the Gaussian step);
    * delete the column r from B[d-1].

Omitting the last leaves a complex whose homology is wrong -- it was the cause of a real
bug during this project, which is why `validate.py` cross-checks two independent
pipelines against each other and against Abrams' model.

WHERE TORSION COMES FROM
------------------------
H_d = ker d_d / im d_{d+1}.  So the free rank in degree d is

        #generators(d) - rank d_d - rank d_{d+1}

and the torsion of H_d is read off the invariant factors of d_{d+1} -- the differential
coming IN, not the one going out.  This is the other mistake that is easy to make, and
it shifts all reported torsion by one degree.
"""

import heapq
from collections import Counter, defaultdict

from .smith import smith


class ChainComplex:
    """A finitely generated chain complex of free abelian groups.

    Attributes
    ----------
    gens : dict
        gens[d] is the set of generators in degree d (any hashable labels).
    B : dict
        B[d] is the matrix of d_d, stored by columns:
        {column generator: {row generator: coefficient}}.
    """

    def __init__(self, gens, B):
        self.gens = {d: set(g) for d, g in gens.items()}
        self.B = {d: {c: dict(col) for c, col in m.items()} for d, m in B.items()}

    def sizes(self):
        """{degree: number of generators}."""
        return {d: len(g) for d, g in sorted(self.gens.items())}

    def check_d2(self):
        """True iff d o d = 0.  A cheap but effective test of the sign conventions."""
        for d in self.B:
            if d - 1 not in self.B:
                continue
            for c, col in self.B[d].items():
                acc = defaultdict(int)
                for r, v in col.items():
                    for r2, v2 in self.B[d - 1].get(r, {}).items():
                        acc[r2] += v * v2
                if any(acc.values()):
                    return False
        return True

    def reduce_fully(self):
        """Cancel every pair of generators joined by a +-1 entry, top degree downwards.

        Modifies the complex in place and returns it, so calls can be chained.
        """
        for d in sorted(self.B, reverse=True):
            mat = self.B[d]
            rows = defaultdict(dict)                       # transpose index
            for c, col in mat.items():
                for r, v in col.items():
                    rows[r][c] = v
            heap = []
            for c, col in mat.items():
                for r, v in col.items():
                    if abs(v) == 1:
                        heapq.heappush(heap, ((len(col) - 1) * (len(rows[r]) - 1), r, c))
            lower = self.B.get(d - 1, {})
            upper = self.B.get(d + 1, {})
            while heap:
                f, r, c = heapq.heappop(heap)
                # the heap is lazy: skip entries invalidated by an earlier pivot
                if (c not in mat or r not in rows or r not in mat[c]
                        or abs(mat[c][r]) != 1):
                    continue
                cur = (len(mat[c]) - 1) * (len(rows[r]) - 1)
                if cur != f:
                    heapq.heappush(heap, (cur, r, c))      # stale cost; re-insert
                    continue
                pv = mat[c][r]
                pcol = dict(mat[c])
                prow = dict(rows[r])
                # clear row r from every other column, using column c
                for c2, v2 in prow.items():
                    if c2 == c:
                        continue
                    fac = -v2 * pv                         # pv = +-1, so 1/pv = pv
                    for r2, v1 in pcol.items():
                        if r2 == r:
                            continue
                        nv = mat[c2].get(r2, 0) + fac * v1
                        if nv:
                            mat[c2][r2] = nv
                            rows[r2][c2] = nv
                        else:
                            mat[c2].pop(r2, None)
                            rows.get(r2, {}).pop(c2, None)
                    mat[c2].pop(r, None)
                for r2 in pcol:
                    rows.get(r2, {}).pop(c, None)
                del mat[c]
                del rows[r]
                # the three updates: ... and the cancelled generators leave the
                # neighbouring differentials too
                for cc, col in list(upper.items()):
                    col.pop(c, None)
                lower.pop(r, None)
                self.gens[d].discard(c)
                self.gens[d - 1].discard(r)
                # fill-in may have created new unit entries
                for c2 in prow:
                    if c2 in mat:
                        for r2, v in mat[c2].items():
                            if abs(v) == 1:
                                heapq.heappush(
                                    heap,
                                    ((len(mat[c2]) - 1) * (len(rows[r2]) - 1), r2, c2))
            self.B[d] = mat
        # scrub dangling references to generators removed while another degree was active
        for d, mat in self.B.items():
            for c in list(mat):
                if c not in self.gens.get(d, ()):
                    del mat[c]
                    continue
                mat[c] = {r: v for r, v in mat[c].items()
                          if r in self.gens.get(d - 1, ())}
        return self

    def dense(self, d):
        """The matrix of d_d as (rows, nrows, ncols), with generators sorted by repr."""
        rows = sorted(self.gens.get(d - 1, ()), key=repr)
        cols = sorted(self.gens.get(d, ()), key=repr)
        ri = {r: i for i, r in enumerate(rows)}
        Mx = [[0] * len(cols) for _ in rows]
        for j, c in enumerate(cols):
            for r, v in self.B.get(d, {}).get(c, {}).items():
                Mx[ri[r]][j] = v
        return Mx, len(rows), len(cols)

    def homology(self, degrees=None):
        """{degree: (free rank, sorted list of non-unit invariant factors)}."""
        degs = sorted(self.gens) if degrees is None else degrees
        cache = {}

        def rk(d):
            if d not in cache:
                if d not in self.B or not self.gens.get(d):
                    cache[d] = (0, [])
                else:
                    Mx, nr, nc = self.dense(d)
                    cache[d] = smith(Mx, nr, nc)
            return cache[d]

        H = {}
        for d in degs:
            n = len(self.gens.get(d, ()))
            r_d, _ = rk(d)
            r_u, div_u = rk(d + 1)                 # torsion of H_d lives in d_{d+1}
            H[d] = (n - r_d - r_u,
                    sorted(x for x in div_u if x not in (0, 1)))
        return H


# =========================================================================
# Presentation
# =========================================================================

def fmt(free, tors):
    """Pretty-print a finitely generated abelian group from (rank, invariant factors).

    >>> fmt(350, [2])
    'Z^350 + Z/2'
    """
    parts = []
    if free:
        parts.append("Z^%d" % free if free > 1 else "Z")
    for t, m in sorted(Counter(tors).items()):
        parts.append("(Z/%d)^%d" % (t, m) if m > 1 else "Z/%d" % t)
    return " + ".join(parts) if parts else "0"


def primary(tors):
    """Invariant factors -> primary decomposition, as a string.

    The invariant factors of a group record it as a divisibility chain, which hides the
    prime structure: Z/6 is one factor but two primary summands.  Since the question
    of interest here is usually *which primes* occur, this is usually the form to print.

    >>> primary([2, 2, 6])
    '(Z/2)^3 + Z/3'
    """
    prim = Counter()
    for t in tors:
        x = t
        p = 2
        while x > 1:
            while p * p <= x and x % p:
                p += 1
            if p * p > x:
                p = x
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            prim[p ** e] += 1
            p = 2
    return " + ".join("(Z/%d)^%d" % (q, m) if m > 1 else "Z/%d" % q
                      for q, m in sorted(prim.items())) or "0"


def odd_part(x):
    """The largest odd divisor of x.  Used to test for odd torsion."""
    while x and x % 2 == 0:
        x //= 2
    return x
