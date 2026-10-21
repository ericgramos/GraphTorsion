"""
Chain complexes of free abelian groups: d^2 test, discrete Morse reduction, homology.

WHAT THIS MODULE PROVIDES
-------------------------
    ChainComplex   a finitely generated chain complex of free abelian groups, with
                   .sizes(), .check_d2(), .reduce_fully(), .dense(d), .homology()
    fmt            print a group from (free rank, invariant factors): 'Z^155 + (Z/2)^3 + Z/6'
    primary        invariant factors -> primary decomposition: [2,2,2,6] -> '(Z/2)^4 + Z/3'
    odd_part       the largest odd divisor of an integer

This is the reference pipeline.  It keeps generators as opaque hashable labels (for the
Swiatkowski complex, the (mono, st) pairs of `gbt.swiatkowski`), which makes it easy to
inspect and to debug but costs a great deal of Python object overhead.  Use it when you
want to see what the generators are; use `gbt.fast` or `gbt.lean` when you want answers
for large complexes.  Practical limit: a few thousand generators in all.

STORAGE CONVENTION
------------------
The matrix of d_d : C_d -> C_{d-1} is stored BY COLUMNS: B[d] maps each degree-d generator
c to its boundary, a dict {degree-(d-1) generator r: coefficient}.  Zero coefficients are
never stored.  gens[d] is the set of degree-d generators.

EXAMPLE
-------
    >>> from gbt import swiatkowski, complete
    >>> C = swiatkowski(complete(5), 2)          # unreduced complex of B_2 K_5
    >>> C.check_d2()
    True
    >>> C.reduce_fully().homology()
    {0: (1, []), 1: (6, [2]), 2: (0, [])}

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
bug during development, which is why `experiments/validate.py` cross-checks independent
pipelines against each other and against Abrams' model.

There is a fourth, silent, consequence: the cancelled generator c must also disappear as
a ROW of d_{d+1}.  Nothing else in d_{d+1} changes.  The reason is that the Gaussian step
replaces each other column c2 of B[d] by c2' = c2 + fac.c, a change of basis in C_d that
alters only the c-coordinate of a chain: sum a_{c2} c2 + a_c c = sum a_{c2} c2' +
(a_c - sum fac.a_{c2}) c.  So the rows of d_{d+1} indexed by surviving generators keep
their entries, and the c-row is simply dropped (it is the coordinate along the cancelled
pair).

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
        {column generator: {row generator: coefficient}}.  A degree d with no entry in
        B is treated as having d_d = 0.

    The constructor copies both arguments, so the caller's dictionaries are not shared.
    Degrees need not start at 0 and need not be contiguous.
    """

    def __init__(self, gens, B):
        """Build the complex from gens = {d: iterable of generators} and
        B = {d: {column generator: {row generator: coefficient}}}; both are copied."""
        self.gens = {d: set(g) for d, g in gens.items()}
        self.B = {d: {c: dict(col) for c, col in m.items()} for d, m in B.items()}

    def sizes(self):
        """{degree: number of generators}, in increasing order of degree."""
        return {d: len(g) for d, g in sorted(self.gens.items())}

    def check_d2(self):
        """True iff d o d = 0.  A cheap but effective test of the sign conventions.

        For every column c of d_d it composes with d_{d-1} entry by entry and checks that
        every coefficient of d_{d-1}(d_d c) vanishes.  Cost: linear in the number of
        nonzero entries times the typical column length.  A wrong Koszul sign almost
        always shows up here, which is why every builder is checked with it in
        `experiments/model_sanity.py`.
        """
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

        Algebraic discrete Morse reduction (Forman, Adv. Math. 134 (1998) 90-145).  For
        each degree d, from the top down, a pair (c, r) with <d c, r> = +-1 is cancelled by
        a unimodular change of basis; this does not change the homology, and it changes the
        Smith form of each differential only by removing a 1.  Pivots are taken in order of
        increasing Markowitz cost (|column c| - 1) * (|row r| - 1), the number of entries
        that the elimination can create, which keeps the matrices sparse.

        Modifies the complex in place (gens and B) and returns it, so calls can be chained:
        C.reduce_fully().homology().  Afterwards no differential has an entry +-1 (every unit
        entry, including those created by fill-in, is queued and eventually used), so what
        is left is the part of the complex that genuinely needs a Smith normal form.
        """
        for d in sorted(self.B, reverse=True):
            mat = self.B[d]
            # rows[r][c] = coefficient: the same matrix indexed by rows, kept in sync with
            # `mat` so that "all columns with an entry in row r" is a dictionary lookup
            rows = defaultdict(dict)                       # transpose index
            for c, col in mat.items():
                for r, v in col.items():
                    rows[r][c] = v
            # heap entries (Markowitz cost, row, column) for every unit entry; costs go stale
            # as the matrix changes and are re-checked when popped
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
                # update 2 of 3: clear row r from every other column c2, using column c
                # (c2 <- c2 - (v2/pv) c; the effect on d_{d+1} is explained in the module
                # docstring: only the c-row there is affected, and it is dropped below)
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
                # update 1 of 3: delete column c (and its row-index entries) from d_d
                for r2 in pcol:
                    rows.get(r2, {}).pop(c, None)
                del mat[c]
                del rows[r]
                # the cancelled generators leave the neighbouring differentials too:
                # c is no longer a row of d_{d+1}, and (update 3 of 3, the easy one to
                # forget) r is no longer a column of d_{d-1}
                for cc, col in list(upper.items()):
                    col.pop(c, None)
                lower.pop(r, None)
                self.gens[d].discard(c)
                self.gens[d - 1].discard(r)
                # fill-in may have created new unit entries in the columns just modified;
                # queue them (duplicates are harmless: stale entries are skipped above)
                for c2 in prow:
                    if c2 in mat:
                        for r2, v in mat[c2].items():
                            if abs(v) == 1:
                                heapq.heappush(
                                    heap,
                                    ((len(mat[c2]) - 1) * (len(rows[r2]) - 1), r2, c2))
            self.B[d] = mat
        # scrub dangling references to generators removed while another degree was active
        # (a row of d_d may still name a (d-1)-generator cancelled while processing d-1)
        for d, mat in self.B.items():
            for c in list(mat):
                if c not in self.gens.get(d, ()):
                    del mat[c]
                    continue
                mat[c] = {r: v for r, v in mat[c].items()
                          if r in self.gens.get(d - 1, ())}
        return self

    def dense(self, d):
        """The matrix of d_d as a dense list of rows: returns (Mx, nrows, ncols).

        Mx[i][j] is the coefficient of the i-th degree-(d-1) generator in the boundary of the
        j-th degree-d generator, both lists sorted by repr() so the result is deterministic.
        Memory is nrows * ncols Python ints, so call this only after `reduce_fully`.
        """
        rows = sorted(self.gens.get(d - 1, ()), key=repr)
        cols = sorted(self.gens.get(d, ()), key=repr)
        ri = {r: i for i, r in enumerate(rows)}
        Mx = [[0] * len(cols) for _ in rows]
        for j, c in enumerate(cols):
            for r, v in self.B.get(d, {}).get(c, {}).items():
                Mx[ri[r]][j] = v
        return Mx, len(rows), len(cols)

    def homology(self, degrees=None):
        """{degree: (free rank, sorted list of non-unit invariant factors)}.

        `degrees` restricts the computation to the listed degrees (default: all).  For each
        degree d it computes the Smith form of d_d and d_{d+1} densely (with `gbt.smith`),
        and returns

            free rank of H_d   =  #generators(d) - rank d_d - rank d_{d+1}
            torsion of H_d     =  the invariant factors of d_{d+1} other than 0 and 1.

        Example: {1: (6, [2])} means H_1 = Z^6 + Z/2.  Dense Smith forms are cubic in the
        matrix size, so on anything but a toy complex call `reduce_fully()` first.
        """
        degs = sorted(self.gens) if degrees is None else degrees
        cache = {}

        def rk(d):
            """(rank, diagonal of the Smith form) of d_d, computed once and cached."""
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

    `free` is the rank of the free part and `tors` a list of invariant factors (> 1, in
    any order; repeats are collected into powers).  The zero group prints as '0'.  The
    factors are printed as given, not split into primes: use `primary` for that.

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

    `tors` is a list of positive integers (any order, 1's are ignored); each is factored
    by trial division into prime powers p^e, and each p^e is one summand Z/p^e.  The
    result lists the summands in increasing order, collecting repeats into powers; an
    empty list gives '0'.  Only the torsion is printed (no free part).

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
    """The largest odd divisor of x (and 0 for x = 0).  Used to test for odd torsion:
    an invariant factor q has an odd prime divisor exactly when odd_part(q) > 1."""
    while x and x % 2 == 0:
        x //= 2
    return x
