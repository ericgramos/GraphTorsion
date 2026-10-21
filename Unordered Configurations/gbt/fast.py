"""
Integer-indexed pipeline for the FULLY REDUCED Swiatkowski complex: ~100x faster.

Same mathematics as `gbt.swiatkowski` + `gbt.chain`, but every generator is an integer
index into a per-degree list rather than a nested tuple, which removes almost all of the
Python object overhead.  K_6 at k = 4 (71,220 cells) takes about two seconds here against
four minutes through the reference pipeline.

Only the fully reduced complex is built (reduce_at = all vertices), because that is what
the sweeps need.  In the fully reduced complex every occupied vertex is moving, hence

        homological degree  ==  number of occupied vertices,

which is why the degree slices can be enumerated independently of each other -- `gbt.lean`
exploits that further to keep memory down.  Vertex states are 0 (empty) or 2+e for the
difference generator h_{v,e} - h_{v,e0(v)} with e != e0(v); see `gbt.swiatkowski` for the
conventions.

Every vertex must have degree >= 1.  A graph with an isolated vertex has to go through
`gbt.swiatkowski` with reduce_at=() instead.  (Leaves are fine: a degree-1 vertex simply
never occurs as an occupied vertex.)

    compute(G, k)  ->  (H, raw_sizes, reduced_sizes)
                       H[d] = (free rank, sorted list of non-unit invariant factors)
                       raw_sizes[d], reduced_sizes[d] = number of degree-d cells before and
                       after the discrete Morse reduction

THE THREE STAGES
----------------
    build     enumerate the cells degree by degree and write every boundary column;
    morse     cancel unit entries (algebraic discrete Morse theory), top degree down;
    homology  dense Smith normal form of what is left (usually a few hundred cells).

EXAMPLE
-------
    >>> from gbt import fast, complete_bip, primary
    >>> H, raw, red = fast.compute(complete_bip(4, 5), 3)
    >>> H[2][0], primary(H[2][1])
    (155, '(Z/2)^4 + Z/3')
    >>> sum(raw.values()), sum(red.values())
    (19976, 178)

PITFALLS
--------
* Memory is the limit: about 6*10^5 cells is comfortable on 7 GB.  Check
  `gbt.reduced_cell_count(G, k)` first and switch to `gbt.lean` above that.
* The generator lists and column indices are internal; nothing here can tell you what a
  surviving generator "is".  For that use the reference pipeline `gbt.swiatkowski`.
"""

import heapq
import itertools
import time
from collections import defaultdict

from .smith import smith
from .swiatkowski import edge_monomials


def build(G, k):
    """The fully reduced weight-k complex, integer-indexed.

    G is a `Graph` with no isolated vertex and k the number of particles.  Returns
    (sizes, B) where sizes[d] is the number of degree-d generators and B[d] is a list of
    columns, each a dict {row index in degree d-1: coefficient}, in the same order as the
    degree-d generator list.  The generator list itself is not returned; internally a
    degree-d generator is (edge monomial, vertex-state tuple) exactly as in
    `gbt.swiatkowski`, with every occupied vertex in a state 2 + e, e != e0(v).

    Cost: time and memory proportional to the number of cells times the degree (each
    boundary has at most 2d terms).
    """
    N, M = G.N, G.M
    # basepoint edge of every vertex: the least incident edge index (same choice as
    # gbt.swiatkowski, so the two builders produce literally the same complex)
    e0 = {v: min(G.inc[v]) for v in range(N)}
    mon = edge_monomials(M, k)

    # cells[d] = the degree-d basis as a list of (edge monomial, vertex state tuple).
    # The outer loop variable j IS the degree, as explained in the module docstring.
    cells = defaultdict(list)
    for j in range(0, min(k, N) + 1):
        for S in itertools.combinations(range(N), j):
            ch = [[2 + e for e in G.inc[v] if e != e0[v]] for v in S]
            if any(not c for c in ch):
                continue                              # a leaf has no difference generator
            for combo in itertools.product(*ch):
                st = [0] * N
                for v, x in zip(S, combo):
                    st[v] = x
                st = tuple(st)
                for m in mon[k - j]:
                    cells[j].append((m, st))

    # idx[d][cell] = position of the cell in cells[d]: how a face is turned into a row index
    idx = {d: {c: i for i, c in enumerate(cs)} for d, cs in cells.items()}
    B = {}
    for d in sorted(cells):
        if d == 0:
            continue
        low = idx[d - 1]
        cols = []
        for (m, st) in cells[d]:
            col = {}
            sgn = 1
            for v in range(N):
                if st[v] < 2:
                    continue
                e = st[v] - 2
                s1 = list(st)
                s1[v] = 0
                s1 = tuple(s1)
                # d h~_{v,e} = (park on e) - (park on the basepoint edge e0(v));
                # there is no sitting term, because every vertex is reduced
                for ee, c in ((e, sgn), (e0[v], -sgn)):
                    # the face: v emptied, one more particle parked on ee (re-sorted so
                    # the monomial is in canonical form and the dictionary lookup hits)
                    r = low[(tuple(sorted(m + (ee,))), s1)]
                    col[r] = col.get(r, 0) + c
                    if col[r] == 0:
                        del col[r]
                sgn = -sgn                            # Koszul sign of the next slot
            cols.append(col)
        B[d] = cols
    return {d: len(c) for d, c in cells.items()}, B


def morse(sizes, B):
    """Unit-pivot (algebraic discrete Morse) elimination, top degree downwards.

    Input: (sizes, B) as returned by `build`.  Cancelling a pair (c in degree d, r in degree
    d-1) joined by a +-1 entry removes the column c from B[d], clears the row r from the
    other columns of B[d], and removes the column r from B[d-1]; the row c of B[d+1] is
    dropped at the end, when the survivors are relabelled (see `gbt.chain` for why
    nothing else in B[d+1] changes).  Pivots are chosen by least Markowitz cost
    (|column| - 1)(|row| - 1) from a lazy heap.

    Returns (new sizes, new B) with the surviving generators relabelled to 0..n-1 in each
    degree, in increasing order of their old index.  Homology and non-unit invariant
    factors are unchanged.  The input B is not modified (it is copied first).
    """
    # columns become dicts keyed by their original index, so deletions keep indices stable;
    # zero columns are dropped immediately (they can never hold a pivot)
    B = {d: {i: dict(v) for i, v in enumerate(cols) if v} for d, cols in B.items()}
    alive = {d: set(range(n)) for d, n in sizes.items()}   # generators not yet cancelled
    for d in sorted(B, reverse=True):
        cols = B[d]
        # row index (transpose) of d_d, kept in sync with `cols`
        rows = defaultdict(dict)
        for c, v in cols.items():
            for r, x in v.items():
                rows[r][c] = x
        # every unit entry is a candidate pivot, keyed by its Markowitz cost
        heap = [((len(v) - 1) * (len(rows[r]) - 1), r, c)
                for c, v in cols.items() for r, x in v.items() if x in (1, -1)]
        heapq.heapify(heap)
        lower = B.get(d - 1)
        while heap:
            f, r, c = heapq.heappop(heap)
            # lazy heap: the entry may be gone, or no longer a unit, after earlier pivots
            v = cols.get(c)
            if v is None or v.get(r) not in (1, -1):
                continue
            rw = rows.get(r)
            if not rw:
                continue
            cur = (len(v) - 1) * (len(rw) - 1)
            if cur != f:
                heapq.heappush(heap, (cur, r, c))    # cost changed: requeue at true cost
                continue
            pv = v[r]
            pcol = v
            prow = dict(rw)
            # clear row r from every other column c2: c2 <- c2 - (v2 / pv) c, and since
            # pv = +-1 the factor -v2/pv is the integer -v2*pv
            for c2, v2 in prow.items():
                if c2 == c:
                    continue
                col2 = cols[c2]
                fac = -v2 * pv
                for r2, v1 in pcol.items():
                    if r2 == r:
                        continue
                    nv = col2.get(r2, 0) + fac * v1
                    if nv:
                        col2[r2] = nv
                        rows[r2][c2] = nv
                    else:
                        col2.pop(r2, None)
                        rows[r2].pop(c2, None)
                col2.pop(r, None)
            # delete the pivot column c from d_d, and the now-empty row r
            for r2 in pcol:
                rows[r2].pop(c, None)
            del cols[c]
            del rows[r]
            alive[d].discard(c)
            alive[d - 1].discard(r)
            if lower is not None:
                lower.pop(r, None)        # the cancelled (d-1)-generator must go too
            # queue unit entries created by fill-in in the columns just modified
            for c2 in prow:
                v3 = cols.get(c2)
                if v3:
                    for r2, x in v3.items():
                        if x in (1, -1):
                            heapq.heappush(
                                heap, ((len(v3) - 1) * (len(rows[r2]) - 1), r2, c2))

    # renumber the survivors 0..n-1 per degree; filtering each column's rows through
    # relab[d-1] is what drops the rows of generators cancelled in degree d-1 (the c-rows)
    relab = {d: {g: i for i, g in enumerate(sorted(a))} for d, a in alive.items()}
    NB = {}
    for d in sorted(B):
        cols = []
        for g in sorted(alive[d]):
            col = B[d].get(g, {})
            cols.append({relab[d - 1][r]: x for r, x in col.items()
                         if r in relab.get(d - 1, {})})
        NB[d] = cols
    return {d: len(a) for d, a in alive.items()}, NB


def homology(sizes, B):
    """{degree: (free rank, sorted list of non-unit invariant factors)}.

    Input: (sizes, B) in the format of `build` / `morse`, normally after `morse`.  Each
    differential is written out densely and handed to `gbt.smith.smith`, so this is only
    sensible on a reduced complex.

    As in `gbt.chain`: the free rank in degree d is #generators - rank d_d - rank d_{d+1},
    and the torsion of H_d is read off the invariant factors of d_{d+1}.
    """
    rk, dv = {}, {}
    for d in sorted(sizes):
        cols = B.get(d)
        if not cols or not sizes.get(d - 1):
            rk[d], dv[d] = 0, []
            continue
        nr, nc = sizes[d - 1], sizes[d]
        Mx = [[0] * nc for _ in range(nr)]
        for j, col in enumerate(cols):
            for r, x in col.items():
                Mx[r][j] = x
        rk[d], dv[d] = smith(Mx, nr, nc)
    top = max(sizes) + 1
    rk[top], dv[top] = 0, []
    return {d: (sizes[d] - rk.get(d, 0) - rk.get(d + 1, 0),
                [x for x in dv.get(d + 1, []) if x not in (0, 1)])
            for d in sorted(sizes)}


def compute(G, k, verbose=False):
    """Build, Morse-reduce, Smith-normalise: the integral homology of B_k G.

    G is a `Graph` with no isolated vertex, k >= 1 the number of particles.  With
    verbose=True one progress line (cells before/after reduction, seconds per stage) is
    printed.  Returns (H, raw sizes, reduced sizes):

        H[d]          (free rank of H_d, sorted list of invariant factors > 1 of tors H_d)
        raw[d]        number of degree-d cells of the fully reduced Swiatkowski complex
        reduced[d]    number of degree-d cells left after the discrete Morse reduction

    for d = 0, ..., min(k, N).  Example: H[2] == (155, [2, 2, 2, 6]) for K_{4,5} at k = 3.
    """
    t0 = time.time()
    sz, B = build(G, k)
    t1 = time.time()
    sz2, B2 = morse(sz, B)
    t2 = time.time()
    H = homology(sz2, B2)
    t3 = time.time()
    if verbose:
        print("   raw %d -> %d  [gen %.0fs morse %.0fs snf %.0fs]"
              % (sum(sz.values()), sum(sz2.values()), t1 - t0, t2 - t1, t3 - t2),
              flush=True)
    return H, sz, sz2
