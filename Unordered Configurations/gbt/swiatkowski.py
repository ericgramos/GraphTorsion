"""
The Swiatkowski chain complex of a graph, with optional partial reduction.

THE MODEL
---------
For a finite simple graph Gamma the Swiatkowski complex is

        S(Gamma)  =  Z[E]  (x)  (x)_{v in V(Gamma)}  S_v ,

bigraded by (homological degree, weight).  Here Z[E] is the polynomial ring on the edge
set -- the exponent of e counts the particles PARKED in the interior of e -- and each
vertex module S_v is free abelian on

    empty_v                 degree 0, weight 0    nothing happens at v
    v                       degree 0, weight 1    a particle SITS at v
    h_{v,e}  (e incident v) degree 1, weight 1    a particle at v MOVES out along e

with differential

        d h_{v,e}  =  e . empty_v  -  v ,

extended over the tensor product as a derivation with Koszul signs.  The weight-k
summand computes H_*(B_k Gamma) for EVERY graph, with no subdivision hypothesis
(Swiatkowski, Colloq. Math. 89 (2001) 69-79; An--Drummond-Cole--Knudsen, Doc. Math. 24
(2019) 1513-1583).  This is why the whole package can work with graphs that have very
few vertices and many edges -- exactly the regime where Abrams' discretised model
explodes.

ENCODING
--------
A basis element is a pair (mono, st):

    mono    a sorted tuple of edge indices, with repetition -- the parked particles
    st      a tuple of length N with
                st[v] = 0        v is empty
                st[v] = 1        a particle sits at v
                st[v] = 2 + e    a particle at v moves along the edge with index e

The homological degree is the number of v with st[v] >= 2, and the weight is
len(mono) + #{v : st[v] >= 1}.  Moving generators are always read in increasing order of
vertex, so the Koszul sign attached to the i-th moving slot (counting from 1) is
(-1)^(i-1).

PARTIAL REDUCTION
-----------------
`reduce_at` is a set R of vertices at which to reduce.  Reducing at v means:

  * choose a basepoint edge e0(v) at v (we take the least incident edge index);
  * change basis in the degree-1 part of S_v to the differences
        h~_{v,e} = h_{v,e} - h_{v,e0(v)},    e != e0(v);
  * cancel the pair (h_{v,e0(v)}, v), which is joined by a unit entry of d.

The sitting state at v then disappears, and

        d h~_{v,e}  =  e - e0(v)         (two parked terms, no sitting term).

Because this cancellation lives inside a single tensor factor, it can be performed at
any subset of vertices independently of the others:

    R = ()            the unreduced complex  (biggest, but every state is visible);
                      this is where the top differential of K_{n+1,b}
                      visibly contains the inclusion matrix W_{n-1,n}(b)
    R = all vertices  the fully reduced complex of An--Drummond-Cole--Knudsen
    R = one side      of a bipartite graph: harmless on the B side of K_{n+1,b}, but
                      reducing on the A side destroys the rows that form the inclusion
                      matrix W_{n-1,n}(b)

A vertex of degree 1 has no non-basepoint edge, so reducing at it deletes all of its
generators; `swiatkowski` therefore skips vertex subsets containing a reduced leaf.

THE DIFFERENTIAL ON A BASIS ELEMENT, WITH SIGNS
-----------------------------------------------
Write a basis element as  m . x_{v_1} x_{v_2} ... x_{v_N}  with x_v the state at v, taken in
increasing order of v.  Only moving states have degree 1, so by the Leibniz rule with
Koszul signs

    d(m . x_{v_1} ... x_{v_N})  =  sum over moving v_i of  (-1)^(number of moving v_j, j < i)
                                   m . x_{v_1} ... d(x_{v_i}) ... x_{v_N},

and d(x_{v_i}) is (e . empty) - (sit at v_i) when v_i is unreduced, or
(e . empty) - (e0(v_i) . empty) when it is reduced.  Multiplying by an edge e just adds e
to the monomial m (Z[E] is central and in degree 0, so it costs no sign).  This is the
loop at the bottom of `swiatkowski`, and `gbt.fast` / `gbt.lean` repeat it with integer
encodings.

EXAMPLE
-------
    >>> from gbt import swiatkowski, complete_bip
    >>> C = swiatkowski(complete_bip(3, 3), 2)
    >>> C.sizes()
    {0: 114, 1: 252, 2: 135}
    >>> C.B[1][((0,), (0, 0, 0, 5, 0, 0))]   # parked on e_0, vertex 3 moving along e_3
    {((0, 3), (0, 0, 0, 0, 0, 0)): 1, ((0,), (0, 0, 0, 1, 0, 0)): -1}

i.e. d(e_0 . h_{3,e_3}) = e_0 e_3 - e_0 . (sit at 3).  Building the unreduced complex is
practical only for small cases (it has `gbt.cell_count(G, k)` generators).
"""

import itertools
from collections import defaultdict

from .chain import ChainComplex


def edge_monomials(M, k):
    """{r: list of all edge monomials of weight r}, for r = 0..k.

    M is the number of edges.  A monomial is a sorted tuple of edge indices with repetition
    allowed (e.g. (0, 0, 3) = e_0^2 e_3: two particles parked on edge 0, one on edge 3), so
    there are C(M + r - 1, r) of weight r.  The lists are in lexicographic order, which
    fixes the generator order used by `gbt.fast` and `gbt.lean`.
    """
    return {r: list(itertools.combinations_with_replacement(range(M), r))
            for r in range(k + 1)}


def swiatkowski(G, k, reduce_at=()):
    """The weight-k Swiatkowski complex of G, reduced at the vertices in `reduce_at`.

    Parameters
    ----------
    G : Graph
        any simple graph; isolated vertices are allowed here (unlike `gbt.fast`).
    k : int
        the number of particles (the weight).
    reduce_at : iterable of vertices
        the set R above.  `()` gives the unreduced complex, `range(G.N)` the fully
        reduced one.

    Returns
    -------
    ChainComplex
        whose generators are the (mono, st) pairs of the module docstring, graded by
        homological degree.

    The construction is direct: enumerate the generators degree by degree, then write
    down each column of the differential from the local rule above.  It is the slowest
    of the three builders in this package but the only one that handles arbitrary
    `reduce_at`, so it is the reference implementation that `fast` and `lean` are
    validated against.
    """
    R = set(reduce_at)
    N, M = G.N, G.M
    # (an isolated vertex in R would make min() fail below; reduce only at vertices with
    # at least one edge)
    # Basepoint edge at each reduced vertex.  Any choice gives an isomorphic complex;
    # the least index makes runs reproducible.
    e0 = {v: min(G.inc[v]) for v in R}
    mon = edge_monomials(M, k)

    # ---- generators ------------------------------------------------------
    gens = defaultdict(set)
    for j in range(0, min(k, N) + 1):             # j = number of occupied vertices
        for S in itertools.combinations(range(N), j):
            choices = []
            for v in S:
                if v in R:
                    # reduced vertex: moving only, and only the non-basepoint differences
                    choices.append([2 + e for e in G.inc[v] if e != e0[v]])
                else:
                    # unreduced vertex: sitting, or moving along any incident edge
                    choices.append([1] + [2 + e for e in G.inc[v]])
            if any(len(c) == 0 for c in choices):
                continue                          # a reduced leaf contributes nothing
            for combo in itertools.product(*choices):
                st = [0] * N
                for v, x in zip(S, combo):
                    st[v] = x
                st = tuple(st)
                d = sum(1 for x in combo if x >= 2)       # degree = number moving
                # j occupied vertices (sitting or moving) carry weight j; the remaining
                # k - j particles are parked, in every possible way
                for m in mon[k - j]:                      # park the remaining weight
                    gens[d].add((m, st))

    # ---- differential ----------------------------------------------------
    B = {}
    for d in sorted(gens):
        if d == 0:
            continue
        mat = {}
        for g in gens[d]:
            m, st = g
            col = {}
            sgn = 1                               # Koszul sign of the current slot

            def add(key, c):
                """Add c to the coefficient of `key` in the current column, dropping zeros."""
                col[key] = col.get(key, 0) + c
                if col[key] == 0:
                    del col[key]

            for v in range(N):                    # slots in increasing vertex order
                if st[v] < 2:
                    continue
                e = st[v] - 2
                s1 = list(st)
                s1[v] = 0
                s1 = tuple(s1)                    # v becomes empty
                # "park on e": present for reduced and unreduced vertices alike
                add((tuple(sorted(m + (e,))), s1), sgn)
                if v in R:
                    # reduced: the second term parks on the basepoint edge
                    add((tuple(sorted(m + (e0[v],))), s1), -sgn)
                else:
                    # unreduced: the second term is "sit at v"
                    s2 = list(st)
                    s2[v] = 1
                    add((m, tuple(s2)), -sgn)
                # the next MOVING slot carries the opposite Koszul sign; sitting and empty
                # slots have degree 0 and are skipped above, so they do not flip it
                sgn = -sgn
            mat[g] = col
        B[d] = mat

    return ChainComplex({d: gens[d] for d in gens}, B)
