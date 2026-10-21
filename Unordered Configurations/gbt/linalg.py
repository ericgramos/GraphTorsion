"""
Integer and modular linear algebra on sparse matrices stored by columns.

Throughout, a matrix is a list of columns, each column a dict {row index: coefficient}
with no zero values.  This is the format the three builders produce, and the format the
"adjoin and test" idiom needs.

WHAT IS HERE
------------
    indexed_columns   extract integer-indexed columns of d_d from a ChainComplex
    kernel_basis      a Z-basis of the kernel lattice of an integer matrix
    sparse_smith      exact (rank, invariant factors): unit-pivot peel, then exact SNF
    peel_sweeps       the same peel driven by sweeps, with no heap: much leaner on memory
    coker_torsion     the torsion of the cokernel (sparse_smith without the rank)
    rank_mod_p_sparse rank over F_p by Markowitz elimination -- runs to completion
    p_primary_count   how many invariant factors are divisible by p, from two ranks
    gf2_kernel        a basis of the kernel over F_2, as bitmasks of column indices
    rank_mod2         rank over F_2 by bitset elimination
    rank_mod_p        rank over F_p for an odd prime p
    betti_mod_p       Betti numbers of a ChainComplex over F_p

THE "ADJOIN AND TEST" IDIOM
---------------------------
To decide whether a cycle z represents a nonzero class in the torsion of H_d, compare

    coker_torsion(columns of d_{d+1})          and
    coker_torsion(columns of d_{d+1} + [z])

Adjoining z as an extra column adds it to the subgroup being quotiented by.  If the
torsion drops, [z] was a nonzero torsion class; if it does not move, z was already a
boundary.  (Precisely: [z] = 0 in H_d exactly when the cokernel is unchanged, rank and
torsion both, since a finitely generated abelian group is not isomorphic to a proper
quotient of itself.  If z is a cycle of infinite order the rank of the image goes up by
one instead -- `sparse_smith` reports both, so check the rank too when that matters.)

The same idiom with a whole family of columns (products, or classes pushed from
a subgraph) tests whether that family GENERATES the torsion.  Used in
`experiments/decomposability_K8.py`, `theta_class.py` and `octahedron_generator.py`.

A warning that cost real time during development: when hunting for a specific element of
ker(d) mod 2 whose halved boundary generates the torsion, do NOT scan the kernel basis in
order, and do not scan prefix sums either -- in the octahedron case the first fourteen
basis vectors and the first twelve prefix sums all landed in the wrong coset.  Random
elements each work with probability 1/2.  See `experiments/octahedron_generator.py`.

EXAMPLE
-------
    >>> from gbt import swiatkowski, complete
    >>> from gbt.linalg import indexed_columns, coker_torsion, kernel_basis
    >>> C = swiatkowski(complete(5), 2)
    >>> d2, rows, _ = indexed_columns(C, 2)
    >>> coker_torsion(d2, verbose=False)                # tors H_1(B_2 K_5)
    [2]
    >>> d1, _, _ = indexed_columns(C, 1)
    >>> len(kernel_basis(d1, len(C.gens[0])))           # rank of the 1-cycles
    166

COST
----
Every routine works on sparse columns in pure Python.  The integral peels and the modular
ranks are roughly linear in the number of nonzeros while the matrix stays sparse, and
degrade as fill-in grows; `kernel_basis` and the dense `smith` at the end of
`sparse_smith` are only meant for matrices of a few thousand columns.
"""

import heapq
from collections import defaultdict

from .smith import smith


# =========================================================================
# Getting columns out of a chain complex
# =========================================================================

def indexed_columns(C, d):
    """Columns of d_d of a ChainComplex, with integer row and column indices.

    C is a `gbt.chain.ChainComplex` and d a degree.  This is the bridge from the
    reference pipeline (generators are labels) to the routines below (generators are
    integers).

    Returns (cols, rowindex, colgens) where

        cols      list of dicts {row index: coefficient}, one per degree-d generator
        rowindex  {degree-(d-1) generator: row index}
        colgens   the degree-d generators, in the order used by `cols`

    Generators are sorted by repr, so the indexing is deterministic across runs.
    """
    rowgens = sorted(C.gens.get(d - 1, ()), key=repr)
    colgens = sorted(C.gens.get(d, ()), key=repr)
    rowindex = {g: i for i, g in enumerate(rowgens)}
    cols = [{rowindex[r]: v for r, v in C.B.get(d, {}).get(g, {}).items()}
            for g in colgens]
    return cols, rowindex, colgens


# =========================================================================
# Integer kernel lattice
# =========================================================================

def kernel_basis(cols, nrows):
    """A Z-basis of {x in Z^ncols : sum_j x_j * col_j = 0}.

    Integer Hermite-style elimination with a transform matrix carried alongside: for each
    row in turn, repeatedly reduce the live columns against the one of least absolute
    value in that row until a single column is left holding it.  Columns that end up zero
    give the kernel vectors, read off the transform.

    Returns a list of dicts {column index: coefficient}.  Exact, no modular lifting.

    `cols` is a list of dicts {row index: coefficient} and `nrows` an upper bound for the
    row indices (rows 0..nrows-1 are scanned).  The transform matrix is unimodular, so the
    columns that become zero form a Z-basis of the kernel lattice, not merely a Q-basis.
    Cost grows quickly with the number of columns (entries can grow and nothing is
    sparse-aware), so use it on matrices of at most a few thousand columns.
    """
    work = [dict(c) for c in cols]               # the columns being reduced
    trans = [{j: 1} for j in range(len(cols))]   # trans[j]: work[j] as a combination of cols
    pivot_of_row = {}                            # row r -> the column that owns its pivot
    for r in range(nrows):
        taken = set(pivot_of_row.values())
        live = [j for j in range(len(cols)) if j not in taken and work[j].get(r)]
        if not live:
            continue
        # Euclid on the entries of row r: reduce every live column by the one of least
        # absolute value until only one column has a nonzero entry in this row
        while len(live) > 1:
            live.sort(key=lambda j: abs(work[j][r]))
            a = live[0]
            for b in live[1:]:
                q = work[b][r] // work[a][r]
                if q:
                    for key, v in work[a].items():
                        nv = work[b].get(key, 0) - q * v
                        if nv:
                            work[b][key] = nv
                        else:
                            work[b].pop(key, None)
                    for key, v in trans[a].items():
                        nv = trans[b].get(key, 0) - q * v
                        if nv:
                            trans[b][key] = nv
                        else:
                            trans[b].pop(key, None)
            live = [j for j in live if work[j].get(r)]
        if live:
            pivot_of_row[r] = live[0]
    piv = set(pivot_of_row.values())
    return [trans[j] for j in range(len(cols)) if j not in piv and not work[j]]


# =========================================================================
# Torsion of a cokernel
# =========================================================================

def sparse_smith(cols, label="", verbose=True):
    """Exact (rank, non-unit invariant factors) of a sparse integer matrix.

    Two stages, for the same reason the homology pipelines have two stages: peel off every
    pair joined by a +-1 entry (a unimodular change of basis, so the invariant factors are
    unchanged apart from dropped 1's), then run the exact Smith normal form on what is
    left.  Without the peel the SNF is hopeless; with it the residual matrix is typically
    a few hundred rows across.

    Each peel step extracts a unit pivot by a unimodular change of basis, so it contributes
    exactly 1 to the rank; the rank returned is therefore (number peeled) + (rank of the
    residual), computed exactly, with no modular reconstruction.

    `cols` is a list of dicts {row index: coefficient}; rows need not be contiguous.
    Returns (rank, list of invariant factors > 1).  `label` names the matrix in the one
    progress line printed when verbose=True.  The input is not modified.

    The peel uses a Markowitz heap (one entry per unit nonzero); for matrices with millions
    of nonzeros use `peel_sweeps` first and call this on the residual.
    """
    m = {c: dict(v) for c, v in enumerate(cols) if v}   # working copy, by column
    rows = defaultdict(dict)                             # the same entries, by row
    for c, v in m.items():
        for r, x in v.items():
            rows[r][c] = x
    # candidate pivots: every +-1 entry, keyed by Markowitz cost (|col|-1)(|row|-1)
    heap = [((len(v) - 1) * (len(rows[r]) - 1), r, c)
            for c, v in m.items() for r, x in v.items() if x in (1, -1)]
    heapq.heapify(heap)
    killed = 0                                           # number of unit pivots used
    while heap:
        f, r, c = heapq.heappop(heap)
        v = m.get(c)
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
        # column operations c2 <- c2 - (v2/pv) c clear row r; then row r and column c
        # form a 1x1 unit block that splits off (a unimodular change of basis on both
        # sides), contributing a 1 to the Smith form and 1 to the rank
        for c2, v2 in prow.items():
            if c2 == c:
                continue
            col2 = m[c2]
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
        del m[c]
        del rows[r]
        killed += 1
        # fill-in can create new unit entries: queue them
        for c2 in prow:
            v3 = m.get(c2)
            if v3:
                for r2, x in v3.items():
                    if x in (1, -1):
                        heapq.heappush(
                            heap, ((len(v3) - 1) * (len(rows[r2]) - 1), r2, c2))

    # what is left has no unit entries; compress it to a dense matrix on the rows that
    # still occur and finish with the exact Smith form
    m = {c: v for c, v in m.items() if v}
    rs = sorted({r for v in m.values() for r in v})
    if not rs:
        if verbose:
            print("      %s: %d cols -> peeled %d, residual empty"
                  % (label or "matrix", len(cols), killed), flush=True)
        return killed, []
    ri = {r: i for i, r in enumerate(rs)}
    Mx = [[0] * len(m) for _ in rs]
    for j, c in enumerate(sorted(m)):
        for r, x in m[c].items():
            Mx[ri[r]][j] = x
    rank, dv = smith(Mx, len(rs), len(m))
    if verbose:
        print("      %s: %d cols -> peeled %d, residual %dx%d"
              % (label or "matrix", len(cols), killed, len(rs), len(m)), flush=True)
    return killed + rank, [d for d in dv if d not in (0, 1)]


def coker_torsion(cols, label="", verbose=True):
    """Non-unit invariant factors of coker(matrix), i.e. the torsion of Z^rows / image.

    A thin wrapper on `sparse_smith` that drops the rank.  This is the workhorse of the
    "adjoin and test" idiom described in the module docstring.

    `cols` are the columns (dicts {row index: coefficient}) of a map into Z^rows; the
    result is the sorted list of invariant factors > 1 of Z^rows / image.  With the columns
    of d_{d+1} this is exactly tors H_d.
    """
    return sparse_smith(cols, label=label, verbose=verbose)[1]


# =========================================================================
# Modular ranks
# =========================================================================

def gf2_kernel(cols):
    """A basis of ker(matrix mod 2), as integers whose bits are column indices.

    `cols` is a list of dicts {row index: coefficient}; only the parity of each entry is
    used.  Each column is turned into a bitmask of its odd rows and reduced against the
    pivots found so far (Gaussian elimination over F_2 with xor), while a second bitmask t
    records which original columns were combined.  A column that reduces to zero gives
    the kernel vector t.

    Each returned integer t encodes the column subset {j : bit j of t is set}; the
    corresponding mod-2 column sum is zero.  Columns are processed left to right, so the
    result is a basis in echelon form with respect to that order -- which is exactly why
    it is a bad place to look for a *generic* kernel element (see the module docstring).
    """
    piv = {}
    basis = []
    for j, col in enumerate(cols):
        v = 0
        for r, x in col.items():
            if x % 2:
                v ^= 1 << r
        t = 1 << j
        while v:
            b = v.bit_length() - 1
            if b not in piv:
                piv[b] = (v, t)
                v = 0
                break
            pv, pt = piv[b]
            v ^= pv
            t ^= pt
        else:
            basis.append(t)
    return basis


def rank_mod2(cols):
    """Rank over F_2 by bitset elimination on the rows.

    Packing each row of the matrix into one Python int and eliminating with xor is about
    as fast as pure Python gets: it is what makes the independent mod-2 cross-check of the
    RAW (unreduced) complexes affordable.
    """
    rows = {}                                  # rows[r] = bitmask of the columns odd in row r
    for j, col in enumerate(cols):
        for r, v in col.items():
            if v % 2:
                rows[r] = rows.get(r, 0) ^ (1 << j)
    piv = {}                                   # leading bit -> reduced row owning it
    rank = 0
    for r, v in rows.items():
        while v:
            b = v.bit_length() - 1
            if b not in piv:
                piv[b] = v
                rank += 1
                break
            v ^= piv[b]
    return rank


def peel_sweeps(cols, nrows=None, verbose=False):
    """Unit-pivot peel driven by cheap repeated sweeps instead of a priority queue.

    `sparse_smith` keeps every candidate pivot in a Markowitz heap.  That is the right
    strategy when the matrix fits comfortably in memory, but the heap holds one tuple per
    nonzero entry and grows further on every fill-in event, and for a matrix with a few
    million nonzeros it is the single largest allocation -- larger than the matrix.

    This routine keeps no queue at all.  It sweeps the remaining unit entries in rounds,
    accepting in each round only pivots whose Markowitz cost (|col|-1)(|row|-1) is within
    the current budget, and recomputing costs on the fly.  The budget starts at zero -- those
    pivots create no fill-in whatsoever -- and whenever a round finds nothing it doubles,
    so a handful of rounds covers any cost that occurs while still taking the cheap pivots
    first.  Each round is linear in the number of nonzeros, and the memory stays at the size
    of the matrix instead of the size of the matrix plus a queue entry per nonzero.

    Returns (residual columns, number of pivots used).  The residual has no +-1 entries
    left, and its invariant factors are those of the original with the units removed, so
    `sparse_smith` on the residual finishes the job.  (The budget sequence is 0, 1, 3, 7,
    15, ...: each round without progress sets budget = 2 * budget + 1, and any progress
    resets it to 0.)  `nrows` is accepted for symmetry with the other routines and is not
    used.  With verbose=True one line is printed per round.
    """
    m = {c: dict(v) for c, v in enumerate(cols) if v}   # working copy, by column
    rows = defaultdict(dict)                             # the same entries, by row
    for c, v in m.items():
        for r, x in v.items():
            rows[r][c] = x
    killed = 0                                           # unit pivots used so far
    budget = 0                                           # largest Markowitz cost accepted
    rounds = 0
    while True:
        progress = False
        rounds += 1
        for c in list(m):
            v = m.get(c)
            if not v:
                continue
            # the cheapest unit entry in this column, by Markowitz cost
            best = None
            for r, x in v.items():
                if x in (1, -1):
                    cost = (len(v) - 1) * (len(rows[r]) - 1)
                    if best is None or cost < best[0]:
                        best = (cost, r)
                        if cost == 0:
                            break
            if best is None or best[0] > budget:
                continue
            _, r = best
            pv = v[r]
            prow = dict(rows[r])
            # the same unit-pivot elimination as in `sparse_smith`: clear row r using
            # column c, then delete both
            for c2, v2 in prow.items():
                if c2 == c:
                    continue
                col2 = m[c2]
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
                if not col2:
                    del m[c2]
            for r2 in v:
                rows[r2].pop(c, None)
            del m[c]
            del rows[r]
            killed += 1
            progress = True
        if progress:
            budget = 0                     # cheap pivots may have reappeared
            if verbose:
                print("        round %d: %d pivots so far, %d columns left"
                      % (rounds, killed, len(m)), flush=True)
            continue
        # nothing at this budget.  If no unit entry is left anywhere we are done;
        # otherwise raise the budget and try again.
        if not any(x in (1, -1) for v in m.values() for x in v.values()):
            break
        budget = budget * 2 + 1
        if verbose:
            print("        round %d: raising the fill budget to %d"
                  % (rounds, budget), flush=True)
    return [v for v in m.values() if v], killed


def rank_mod_p_sparse(cols, p, verbose=False):
    """Rank over F_p by sparse elimination with a Markowitz pivot order.

    The same machinery as the integral unit-pivot peel in `sparse_smith`, with one crucial
    difference: over a field EVERY nonzero entry is a unit, so there is never a point at
    which the elimination has to stop.  It therefore runs all the way down and the number of
    pivots it performs IS the rank.

    This is the cheap way to detect torsion in a complex too large for the integral
    pipeline.  By the identity "torsion is a cokernel" the torsion of H_d is the torsion of
    coker d_{d+1}, and the number of invariant factors of a matrix divisible by a prime p is

        rank over Q  -  rank over F_p

    so two rank computations settle the p-primary part without any Smith normal form, any
    discrete Morse reduction, or any integer arithmetic beyond residues mod p.  Coefficients
    cannot grow, which is what makes it fit in memory where the integral peel does not.

    `cols` is a list of dicts {row index: coefficient} with integer entries (reduced mod p
    here) and p must be prime (the inverse is computed as x^(p-2)).  Returns the rank.
    With verbose=True a progress line is printed every 50,000 pivots.
    """
    m = {}
    for c, col in enumerate(cols):
        v = {r: x % p for r, x in col.items() if x % p}
        if v:
            m[c] = v
    rows = defaultdict(dict)
    for c, v in m.items():
        for r, x in v.items():
            rows[r][c] = x
    heap = [((len(v) - 1) * (len(rows[r]) - 1), r, c)
            for c, v in m.items() for r, x in v.items()]
    heapq.heapify(heap)
    rank = 0
    while heap:
        f, r, c = heapq.heappop(heap)
        v = m.get(c)
        if v is None or r not in v:
            continue
        rw = rows.get(r)
        if not rw:
            continue
        cur = (len(v) - 1) * (len(rw) - 1)
        if cur != f:
            heapq.heappush(heap, (cur, r, c))       # stale Markowitz cost
            continue
        inv = pow(v[r], p - 2, p)                    # p is prime, so this is 1/pivot
        prow = dict(rw)
        for c2, v2 in prow.items():
            if c2 == c:
                continue
            col2 = m[c2]
            fac = (-v2 * inv) % p
            for r2, v1 in v.items():
                if r2 == r:
                    continue
                nv = (col2.get(r2, 0) + fac * v1) % p
                if nv:
                    col2[r2] = nv
                    rows[r2][c2] = nv
                else:
                    col2.pop(r2, None)
                    rows[r2].pop(c2, None)
            col2.pop(r, None)
            if not col2:
                del m[c2]
        for r2 in v:
            rows[r2].pop(c, None)
        del m[c]
        del rows[r]
        rank += 1
        for c2 in prow:
            v3 = m.get(c2)
            if v3:
                for r2, x in v3.items():
                    heapq.heappush(heap, ((len(v3) - 1) * (len(rows[r2]) - 1), r2, c2))
        if verbose and rank % 50000 == 0:
            print("        %d pivots, %d columns left" % (rank, len(m)), flush=True)
    return rank


def p_primary_count(cols, p=2, big_prime=1000003, verbose=False):
    """How many invariant factors of a matrix are divisible by p.

    Equals rank over F_{big_prime} minus rank over F_p.  By "torsion is a cokernel"
    this is the number of p-primary summands of the torsion of H_d when `cols` are the
    columns of d_{d+1}: nonzero means there IS p-torsion.  It does not say which powers of
    p occur -- for that the integral pipeline is needed -- but it is the cheapest possible
    yes/no, and the count itself is the multiplicity in the primary decomposition.

    `big_prime` stands in for Q.  The rank over any F_q is at most the rational rank, with
    equality unless q divides some invariant factor; so the answer is exact as long as
    big_prime divides no invariant factor, which is overwhelmingly likely for a prime of
    this size but not guaranteed.  To be sure, run it again with a second large prime and
    check that the two answers agree.  `cols` is a list of dicts {row: coefficient}.
    """
    r_big = rank_mod_p_sparse(cols, big_prime, verbose=verbose)
    r_p = rank_mod_p_sparse(cols, p, verbose=verbose)
    if verbose:
        print("        rank over F_%d = %d, over F_%d = %d"
              % (big_prime, r_big, p, r_p), flush=True)
    return r_big - r_p


def rank_mod_p(cols, p):
    """Rank over F_p, p an odd prime, by sparse elimination with Markowitz-style ordering.

    Columns are processed from sparsest to densest, which keeps fill-in down on the
    matrices that arise here.  Each column is reduced against the stored pivot rows keyed
    by their smallest row index (a row-echelon form built column by column); a column that
    does not reduce to zero contributes 1 to the rank.  `cols` is a list of dicts
    {row: coefficient}; p must be prime (it also works for p = 2, but `rank_mod2` is
    faster there).
    """
    work = {c: {r: v % p for r, v in col.items() if v % p}
            for c, col in enumerate(cols)}
    work = {c: v for c, v in work.items() if v}
    piv = {}
    rank = 0
    for c in sorted(work, key=lambda c: len(work[c])):
        v = dict(work[c])
        while v:
            r = min(v)
            if r not in piv:
                piv[r] = v
                rank += 1
                break
            pr = piv[r]
            f = v[r] * pow(pr[r], p - 2, p) % p
            for key, val in pr.items():
                nv = (v.get(key, 0) - f * val) % p
                if nv:
                    v[key] = nv
                elif key in v:
                    del v[key]
    return rank


def betti_mod_p(C, p):
    """Betti numbers of a ChainComplex over F_p: {degree: dim H_d(C; F_p)}.

    Computed straight from the matrices, with no discrete Morse reduction and no Smith
    normal form, so this is an independent route to the same information: comparing
    dim H_d(F_2) with dim H_d(F_q) for a large prime q detects 2-torsion without any of
    the integral machinery being involved.

    C is a `gbt.chain.ChainComplex`; p a prime.  Uses `rank_mod2` for p = 2 and `rank_mod_p`
    otherwise.  By universal coefficients, dim H_d(C; F_p) = f_d + (number of invariant
    factors of tors H_d divisible by p) + (the same for tors H_{d-1}), where f_d is the
    free rank of H_d; `experiments/crosscheck_mod_p.py` tests exactly this.
    """
    degs = sorted(C.gens)
    rk = {}
    for d in degs + [max(degs) + 1]:
        cols, _, _ = indexed_columns(C, d)
        rk[d] = (rank_mod2(cols) if p == 2 else rank_mod_p(cols, p)) if cols else 0
    return {d: len(C.gens[d]) - rk.get(d, 0) - rk.get(d + 1, 0) for d in degs}
