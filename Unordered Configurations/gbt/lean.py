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

    compute(G, k)  ->  (H, raw_sizes, reduced_sizes)
"""

import heapq
import itertools
import time
from collections import defaultdict

from .smith import smith
from .swiatkowski import edge_monomials

#: Base for the packed-integer cell keys.  Must exceed N*M for the graph in hand; 2^11
#: covers every case in data/ (the largest is K_13, with 13*78 = 1014 < 2048).
BASE = 1 << 11


def _check_base(G):
    if G.N * G.M + G.M + 1 >= BASE:
        raise ValueError("BASE too small for a graph with N=%d, M=%d; raise gbt.lean.BASE"
                         % (G.N, G.M))


def slice_cells(G, k, d, mon):
    """The packed integer keys of the degree-d cells.

    A cell has d occupied (moving) vertices and k-d parked particles.  The key packs the
    d moving slots as (v*M + e + 1) in digits 0..d-1, then the parked edge indices as
    (e + 1) in digits d..k-1.  Moving slots are written in increasing vertex order and
    parked edges in increasing index order, so the key is canonical.
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
    """Packed key -> (list of (vertex, edge) moving slots, tuple of parked edges)."""
    mov = []
    for i in range(d):
        c = (key // BASE ** i) % BASE - 1
        mov.append((c // M, c % M))
    mon = tuple(sorted((key // BASE ** (d + j)) % BASE - 1 for j in range(kd)))
    return mov, mon


def encode(mov, mon, M):
    """(moving slots, parked edges) -> packed key.  Inverse of `decode`."""
    key = 0
    for i, (v, e) in enumerate(sorted(mov)):
        key += (v * M + e + 1) * BASE ** i
    for j, ei in enumerate(sorted(mon)):
        key += (ei + 1) * BASE ** (len(mov) + j)
    return key


def boundary(G, k, d, mon, idx_low):
    """Columns of d_d as dicts {row index in degree d-1: coefficient}, in slice order.

    `idx_low` maps the degree-(d-1) packed keys to their indices.
    """
    M = G.M
    e0 = [min(G.inc[v]) for v in range(G.N)]
    cols = []
    for key in slice_cells(G, k, d, mon):
        mov, mo = decode(key, d, k - d, M)
        col = {}
        sgn = 1
        for (v, e) in mov:
            rest = [p for p in mov if p[0] != v]
            for ee, c in ((e, sgn), (e0[v], -sgn)):
                r = idx_low[encode(rest, mo + (ee,), M)]
                col[r] = col.get(r, 0) + c
                if col[r] == 0:
                    del col[r]
            sgn = -sgn
        cols.append(col)
    return cols


def reduce_top_down(G, k, verbose=False, stop_degree=1):
    """Build and Morse-reduce one degree at a time.  Returns (raw sizes, sizes, B).

    `stop_degree` is the lowest differential to process: the peel runs over
    d = top, top-1, ..., stop_degree and the lower ones are left alone.  The default 1
    reduces everything, which is what a full homology computation needs.

    Stopping early is useful because of Lemma "torsion is a cokernel": the torsion of
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
    dead = {d: set() for d in range(stop - 1, top + 1)}
    kept = {}

    for d in range(top, stop - 1, -1):
        idx_low = {key: i for i, key in enumerate(slices[d - 1])}
        cols = boundary(G, k, d, mon, idx_low)
        del idx_low                                  # the big dictionary goes immediately
        cols = {c: v for c, v in enumerate(cols) if v and c not in dead[d]}
        for c in cols:                               # drop rows killed at degree d+1
            for r in list(cols[c]):
                if r in dead[d - 1]:
                    del cols[c][r]
        rows = defaultdict(dict)
        for c, v in cols.items():
            for r, x in v.items():
                rows[r][c] = x
        heap = [((len(v) - 1) * (len(rows[r]) - 1), r, c)
                for c, v in cols.items() for r, x in v.items() if x in (1, -1)]
        heapq.heapify(heap)
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
            prow = dict(rw)
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
            dead[d].add(c)
            dead[d - 1].add(r)
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
    """{degree: (free rank, non-unit invariant factors)}; see `gbt.fast.homology`."""
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

    First, Lemma "torsion is a cokernel": tors H_d = tors coker d_{d+1}.  So only ONE
    differential is needed -- the degree-(d+1) slice and the degree-d slice of the complex --
    and the differentials below it never enter.  Nothing has to be reduced first either: the
    identity holds in the raw complex.

    Second, the peel is driven by `gbt.linalg.peel_sweeps` rather than by a Markowitz heap.
    The heap holds one entry per nonzero and grows on every fill-in event, and at this size
    it is a larger allocation than the matrix itself; sweeping instead keeps the memory at
    the size of the matrix.  The icosahedron at k = 4 (1.3 million cells) exhausts 7 GB
    through the heap and finishes in about three minutes this way.

    Returns the list of non-unit invariant factors (empty if H_d is torsion-free).
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
    """Build, Morse-reduce, Smith-normalise.  Returns (H, raw sizes, reduced sizes)."""
    t0 = time.time()
    raw, sz, B = reduce_top_down(G, k, verbose)
    t1 = time.time()
    H = homology(sz, B)
    if verbose:
        print("      raw %s -> %s [morse %.0fs snf %.0fs]"
              % (format(sum(raw.values()), ","), format(sum(sz.values()), ","),
                 t1 - t0, time.time() - t1), flush=True)
    return H, raw, sz
