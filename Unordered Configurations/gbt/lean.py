"""
Memory-lean pipeline: the same answers as `gbt.fast`, in roughly a third of the memory.

Two changes make the difference.

  1.  Cells are packed into single Python ints (base 2^11 digits) instead of nested
      tuples, so the index dictionaries are much smaller.

  2.  The boundary matrices are built one degree at a time, from the top down.  This is
      possible because in the fully reduced complex the homological degree equals the
      number of occupied vertices, so each degree slice can be enumerated on its own
      without the others being in memory.  Only two index dictionaries -- degree d and
      degree d-1 -- are ever live at once.

Use this for complexes of roughly 1-4 million cells, where `gbt.fast` runs out of memory.
Above that the Python peel itself becomes the bottleneck (see the README caveats).

    compute(G, k)                  ->  (H, raw_sizes, reduced_sizes)   as in gbt.fast
    torsion_in_degree(G, k, d)     ->  invariant factors > 1 of tors H_d(B_k G)

CELL KEYS
---------
A degree-d cell of the fully reduced weight-k complex has d moving slots (v, e), one per
occupied vertex v with a non-basepoint edge e, and k - d parked particles.  It is stored as
one integer written in base BASE:

    digit i       (i = 0..d-1)       v_i * M + e_i + 1   for the slots in increasing v
    digit d + j   (j = 0..k-d-1)     e'_j + 1            for the parked edges, sorted

The "+ 1" keeps every digit nonzero, so keys of different shapes never collide.  The digit
order is the canonical order (vertices increasing, monomial sorted), which is also the
order in which the Koszul signs are assigned: the i-th moving slot gets (-1)^i.

SINGLE-DEGREE TORSION
---------------------
Because tors H_d = tors coker d_{d+1} (see `torsion_in_degree`), the torsion of one degree
needs only the slices of degrees d and d+1.  This is how the largest computations in the
repository were done (for example K_{5,7} and K_{5,8} at k = 4, degree 3).

EXAMPLE
-------
    >>> from gbt import lean, complete_bip, octahedron, primary
    >>> primary(lean.torsion_in_degree(complete_bip(4, 5), 3, 2))
    '(Z/2)^4 + Z/3'
    >>> H, raw, red = lean.compute(octahedron(), 4)
    >>> H[2]
    (95, [2])
"""

import heapq
import itertools
import time
from collections import defaultdict

from .smith import smith
from .swiatkowski import edge_monomials

#: Base for the packed-integer cell keys.  Must exceed N*M for the graph in hand; 2^11
#: covers every graph attempted here (the largest is K_13, with 13*78 = 1014 < 2048).
BASE = 1 << 11


def _check_base(G):
    """Raise ValueError unless every digit of a cell key of G fits below BASE.

    The largest digit is N*M (a moving slot (N-1, M-1) gives (N-1)*M + M), so BASE must
    exceed N*M; the test is slightly more conservative than that.
    """
    if G.N * G.M + G.M + 1 >= BASE:
        raise ValueError("BASE too small for a graph with N=%d, M=%d; raise gbt.lean.BASE"
                         % (G.N, G.M))


def slice_cells(G, k, d, mon):
    """The packed integer keys of the degree-d cells, as a list.

    A cell has d occupied (moving) vertices and k-d parked particles.  The key packs the
    d moving slots as (v*M + e + 1) in digits 0..d-1, then the parked edge indices as
    (e + 1) in digits d..k-1.  Moving slots are written in increasing vertex order and
    parked edges in increasing index order, so the key is canonical.

    `mon` is `edge_monomials(G.M, k)`, passed in so that it is built only once.  The order
    of the list (vertex subsets in lexicographic order, then edge choices, then
    monomials) is the column order used by `boundary` and the row order used for the
    degree below.  Leaves (degree-1 vertices) have no non-basepoint edge and never occur.
    """
    N, M = G.N, G.M
    e0 = [min(G.inc[v]) for v in range(N)]
    out = []
    for S in itertools.combinations(range(N), d):
        ch = [[e for e in G.inc[v] if e != e0[v]] for v in S]
        if any(not c for c in ch):
            continue
        for combo in itertools.product(*ch):
            head = 0
            for i, (v, e) in enumerate(zip(S, combo)):
                head += (v * M + e + 1) * BASE ** i
            for m in mon[k - d]:
                key = head
                for j, ei in enumerate(m):
                    key += (ei + 1) * BASE ** (d + j)
                out.append(key)
    return out


def decode(key, d, kd, M):
    """Packed key -> (list of (vertex, edge) moving slots, tuple of parked edges).

    d is the number of moving slots, kd = k - d the number of parked particles, M the
    number of edges.  The slots come out in increasing vertex order and the parked edges
    sorted, i.e. in canonical order.
    """
    mov = []
    for i in range(d):
        c = (key // BASE ** i) % BASE - 1
        mov.append((c // M, c % M))
    mon = tuple(sorted((key // BASE ** (d + j)) % BASE - 1 for j in range(kd)))
    return mov, mon


def encode(mov, mon, M):
    """(moving slots, parked edges) -> packed key.  Inverse of `decode`.

    `mov` is an iterable of (vertex, edge) pairs and `mon` of edge indices, in any order:
    both are sorted here, so the result is the canonical key that `slice_cells` produces.
    """
    key = 0
    for i, (v, e) in enumerate(sorted(mov)):
        key += (v * M + e + 1) * BASE ** i
    for j, ei in enumerate(sorted(mon)):
        key += (ei + 1) * BASE ** (len(mov) + j)
    return key


def boundary(G, k, d, mon, idx_low):
    """Columns of d_d as dicts {row index in degree d-1: coefficient}, in slice order.

    `idx_low` maps the degree-(d-1) packed keys to their indices (built from
    `slice_cells(G, k, d - 1, mon)`).  For a cell with moving slots (v_1, e_1), ...,
    (v_d, e_d) in increasing vertex order, the boundary is

        sum_i (-1)^(i-1) [ (v_i emptied, park on e_i) - (v_i emptied, park on e0(v_i)) ]

    exactly as in `gbt.fast.build`.
    """
    M = G.M
    e0 = [min(G.inc[v]) for v in range(G.N)]
    cols = []
    for key in slice_cells(G, k, d, mon):
        mov, mo = decode(key, d, k - d, M)
        col = {}
        sgn = 1                                       # Koszul sign of the current slot
        for (v, e) in mov:                            # slots in increasing vertex order
            rest = [p for p in mov if p[0] != v]      # the cell with v emptied
            for ee, c in ((e, sgn), (e0[v], -sgn)):
                r = idx_low[encode(rest, mo + (ee,), M)]
                col[r] = col.get(r, 0) + c
                if col[r] == 0:
                    del col[r]
            sgn = -sgn                                # every slot is moving: flip each time
        cols.append(col)
    return cols


def reduce_top_down(G, k, verbose=False, stop_degree=1):
    """Build and Morse-reduce one degree at a time.  Returns (raw sizes, sizes, B).

    raw[d] is the number of degree-d cells before reduction, sizes[d] after, and B[d] the
    reduced d_d as a list of columns {row index: coefficient}, generators relabelled to
    0..n-1 per degree (the format of `gbt.fast.morse`).  Only the slice of degree d and the
    index dictionary of degree d-1 are built at each step, which is the memory saving.

    `stop_degree` is the lowest differential to process: the peel runs over
    d = top, top-1, ..., stop_degree and the lower ones are left alone.  The default 1
    reduces everything, which is what a full homology computation needs.

    Stopping early is useful because of the identity "torsion is a cokernel": the torsion of
    H_d depends only on d_{d+1}, so stopping at d+1 is enough to get it, and the peels at
    degrees <= d -- often the most expensive ones, since the low slices are the biggest --
    can be skipped entirely.  `torsion_in_degree` does exactly that.  The returned `sizes`
    and `B` are then only meaningful in degrees >= stop_degree - 1.
    """
    _check_base(G)
    N = G.N
    top = min(k, N)
    mon = edge_monomials(G.M, k)
    stop = max(1, min(stop_degree, top))
    # only the slices at or above stop_degree - 1 are ever needed
    slices = {d: slice_cells(G, k, d, mon) for d in range(stop - 1, top + 1)}
    raw = {d: len(s) for d, s in slices.items()}
    # dead[d] = indices (in slice order) of degree-d cells already cancelled
    dead = {d: set() for d in range(stop - 1, top + 1)}
    kept = {}                                        # kept[d] = the reduced columns of d_d

    for d in range(top, stop - 1, -1):
        idx_low = {key: i for i, key in enumerate(slices[d - 1])}
        cols = boundary(G, k, d, mon, idx_low)
        del idx_low                                  # the big dictionary goes immediately
        # a degree-d cell cancelled as a ROW while processing d+1 must not be a column
        # here: this is the "delete column r from d_{d-1}" update of the Morse step
        cols = {c: v for c, v in enumerate(cols) if v and c not in dead[d]}
        # rows of degree d-1 are only killed while processing degree d (below), so at this
        # point dead[d-1] is still empty and this loop is a safeguard that removes nothing
        for c in cols:
            for r in list(cols[c]):
                if r in dead[d - 1]:
                    del cols[c][r]
        rows = defaultdict(dict)
        for c, v in cols.items():
            for r, x in v.items():
                rows[r][c] = x
        # the same lazy-heap Markowitz peel as gbt.fast.morse; see the comments there
        heap = [((len(v) - 1) * (len(rows[r]) - 1), r, c)
                for c, v in cols.items() for r, x in v.items() if x in (1, -1)]
        heapq.heapify(heap)
        while heap:
            f, r, c = heapq.heappop(heap)
            v = cols.get(c)
            if v is None or v.get(r) not in (1, -1):
                continue                             # stale: entry gone or no longer a unit
            rw = rows.get(r)
            if not rw:
                continue
            cur = (len(v) - 1) * (len(rw) - 1)
            if cur != f:
                heapq.heappush(heap, (cur, r, c))    # stale cost: requeue at the true cost
                continue
            pv = v[r]
            prow = dict(rw)
            # clear row r from every other column, using the pivot column c
            for c2, v2 in prow.items():
                if c2 == c:
                    continue
                col2 = cols[c2]
                fac = -v2 * pv
                for r2, v1 in v.items():
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
            for r2 in v:
                rows[r2].pop(c, None)
            del cols[c]
            del rows[r]
            dead[d].add(c)                           # c leaves C_d (and later the rows of d_{d+1})
            dead[d - 1].add(r)                       # r leaves C_{d-1} (and the columns of d_{d-1})
            for c2 in prow:
                v3 = cols.get(c2)
                if v3:
                    for r2, x in v3.items():
                        if x in (1, -1):
                            heapq.heappush(
                                heap, ((len(v3) - 1) * (len(rows[r2]) - 1), r2, c2))
        kept[d] = cols
        del rows
        if verbose:
            print("      degree %d: %9s -> %6s alive"
                  % (d, format(raw[d], ","), format(raw[d] - len(dead[d]), ",")),
                  flush=True)

    sizes = {d: raw[d] - len(dead[d]) for d in raw}
    # renumber survivors; filtering rows through relab[d-1] drops the rows of the cells
    # cancelled in degree d-1, the last of the Morse updates
    relab = {d: {g: i for i, g in enumerate(sorted(set(range(raw[d])) - dead[d]))}
             for d in raw}
    B = {}
    for d in kept:
        cols = []
        for g in sorted(relab[d]):
            col = kept[d].get(g, {})
            cols.append({relab[d - 1][r]: x for r, x in col.items() if r in relab[d - 1]})
        B[d] = cols
    return raw, sizes, B


def homology(sizes, B):
    """{degree: (free rank, sorted list of non-unit invariant factors)}.

    Identical to `gbt.fast.homology` (dense Smith form of each reduced differential); kept
    here so that the module has no import dependence on `gbt.fast`.
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


def torsion_in_degree(G, k, d, verbose=False):
    """The invariant factors of the torsion of H_d(B_k G), as cheaply as possible.

    Two facts make this far cheaper than a full homology computation.

    First, "torsion is a cokernel": tors H_d = tors coker d_{d+1}.  (The sequence
    0 -> H_d -> coker d_{d+1} -> im d_d -> 0 is exact and splits, because im d_d is a
    subgroup of a free group.)  So only ONE
    differential is needed -- the degree-(d+1) slice and the degree-d slice of the complex --
    and the differentials below it never enter.  Nothing has to be reduced first either: the
    identity holds in the raw complex.

    Second, the peel is driven by `gbt.linalg.peel_sweeps` rather than by a Markowitz heap.
    The heap holds one entry per nonzero and grows on every fill-in event, and at this size
    it is a larger allocation than the matrix itself; sweeping instead keeps the memory at
    the size of the matrix.  The icosahedron at k = 4 (1.3 million cells) exhausts 7 GB
    through the heap and finishes in about three minutes this way.

    G is a `Graph` with no isolated vertex, k the number of particles, d >= 0 the degree.
    Returns the list of invariant factors > 1 of tors H_d(B_k G), in divisibility order
    (empty if H_d is torsion-free); `gbt.primary` turns it into a primary decomposition.
    The free rank is NOT computed.  Memory is roughly that of the two slices, whose sizes
    `gbt.sizes.reduced_cell_count_by_degree(G, k)` gives in advance.
    """
    from .linalg import peel_sweeps, sparse_smith
    _check_base(G)
    t0 = time.time()
    mon = edge_monomials(G.M, k)
    low = {key: i for i, key in enumerate(slice_cells(G, k, d, mon))}
    nrows = len(low)
    cols = boundary(G, k, d + 1, mon, low)
    del low
    if verbose:
        print("      raw d_%d: %s columns over %s rows [%.0fs]"
              % (d + 1, format(len(cols), ","), format(nrows, ","), time.time() - t0),
              flush=True)
    if not cols:
        return []
    res, killed = peel_sweeps(cols)
    del cols
    if verbose:
        print("      peeled %s pivots, residual %s columns [%.0fs]"
              % (format(killed, ","), format(len(res), ","), time.time() - t0), flush=True)
    if not res:
        return []
    _, divs = sparse_smith(res, verbose=False)
    return [x for x in divs if x not in (0, 1)]


def compute(G, k, verbose=False):
    """Build, Morse-reduce, Smith-normalise.  Returns (H, raw sizes, reduced sizes).

    Same arguments and the same return format as `gbt.fast.compute`, and the same answers;
    only the memory profile differs.  Raises ValueError if the graph is too large for
    BASE (see `_check_base`).
    """
    t0 = time.time()
    raw, sz, B = reduce_top_down(G, k, verbose)
    t1 = time.time()
    H = homology(sz, B)
    if verbose:
        print("      raw %s -> %s [morse %.0fs snf %.0fs]"
              % (format(sum(raw.values()), ","), format(sum(sz.values()), ","),
                 t1 - t0, time.time() - t1), flush=True)
    return H, raw, sz
