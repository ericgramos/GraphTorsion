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
`gbt.swiatkowski` with reduce_at=() instead.

    compute(G, k)  ->  (H, raw_sizes, reduced_sizes)
                       H[d] = (free rank, non-unit invariant factors)
"""

import heapq
import itertools
import time
from collections import defaultdict

from .smith import smith
from .swiatkowski import edge_monomials


def build(G, k):
    """The fully reduced weight-k complex, integer-indexed.

    Returns (sizes, B) where sizes[d] is the number of degree-d generators and B[d] is a
    list of columns, each a dict {row index in degree d-1: coefficient}, in the same order
    as the degree-d generator list.
    """
    N, M = G.N, G.M
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

    Cancelling a pair (c in degree d, r in degree d-1) removes the column c from B[d],
    clears the row r from the other columns of B[d], and removes the column r from
    B[d-1].  Returns (new sizes, new B) with generators relabelled to 0..n-1 in each
    degree.
    """
    B = {d: {i: dict(v) for i, v in enumerate(cols) if v} for d, cols in B.items()}
    alive = {d: set(range(n)) for d, n in sizes.items()}
    for d in sorted(B, reverse=True):
        cols = B[d]
        rows = defaultdict(dict)
        for c, v in cols.items():
            for r, x in v.items():
                rows[r][c] = x
        heap = [((len(v) - 1) * (len(rows[r]) - 1), r, c)
                for c, v in cols.items() for r, x in v.items() if x in (1, -1)]
        heapq.heapify(heap)
        lower = B.get(d - 1)
        while heap:
            f, r, c = heapq.heappop(heap)
            v = cols.get(c)
            if v is None or v.get(r) not in (1, -1):
                continue
            rw = rows.get(r)
            if not rw:
                continue
            cur = (len(v) - 1) * (len(rw) - 1)
            if cur != f:
                heapq.heappush(heap, (cur, r, c))
                continue
            pv = v[r]
            pcol = v
            prow = dict(rw)
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
            for r2 in pcol:
                rows[r2].pop(c, None)
            del cols[c]
            del rows[r]
            alive[d].discard(c)
            alive[d - 1].discard(r)
            if lower is not None:
                lower.pop(r, None)        # the cancelled (d-1)-generator must go too
            for c2 in prow:
                v3 = cols.get(c2)
                if v3:
                    for r2, x in v3.items():
                        if x in (1, -1):
                            heapq.heappush(
                                heap, ((len(v3) - 1) * (len(rows[r2]) - 1), r2, c2))

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
    """{degree: (free rank, non-unit invariant factors)}.

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
    """Build, Morse-reduce, Smith-normalise.  Returns (H, raw sizes, reduced sizes)."""
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
