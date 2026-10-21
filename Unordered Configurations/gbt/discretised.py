"""
Abrams' discretised configuration space -- an independent model, used only to cross-check.

For a graph G, the discretised configuration space UD^k(G) is the cubical complex whose cells
are the k-element sets of pairwise disjoint closed cells of G (so each chosen cell is a
vertex or a whole edge, and no two share an endpoint).  The homological degree of a cell is
the number of edges it uses.

Abrams' theorem (PhD thesis, Berkeley 2000), as sharpened by Prue--Scrimshaw (Topology
Appl. 178 (2014) 136-145), says UD^k(G) is a deformation retract of UConf_k(G) provided G
is sufficiently subdivided:

    * every path between two essential vertices (degree != 2) has length at least k-1, and
    * every cycle has length at least k+1.

`gbt.subdivide` arranges that; the caller is responsible for checking the hypothesis, as
this module does not verify it.

This model is exponential in k and in the number of edges, so it is unusable for most of
the computations here -- the subdivided octahedron at k=4 already has 36 edges.  Its
value is as a completely independent implementation: it shares no code with the
Swiatkowski builders beyond the `Graph` class, so agreement between the two is strong
evidence that the sign conventions and the reduction are right.  See
`experiments/model_sanity.py`.

EXAMPLE
-------
    >>> from gbt import abrams, complete, subdivide
    >>> C = abrams(complete(5), 2)       # K_5 needs no subdivision for k = 2
    >>> C.sizes(), C.reduce_fully().homology()[1]
    ({0: 10, 1: 30, 2: 15}, (6, [2]))

UD^2(K_5) is the closed non-orientable surface of Euler characteristic 10 - 30 + 15 = -5,
so H_1 = Z^6 + Z/2 (Abrams).

PITFALL
-------
Without enough subdivision the answer is the homology of a DIFFERENT space and nothing
will warn you.  `experiments/verify_headlines.py` lists subdivisions that suffice for the
cases it checks.
"""

from collections import defaultdict

from .chain import ChainComplex


def abrams(G, k):
    """The discretised configuration space UD^k(G) as a ChainComplex.

    G is a `Graph` (already subdivided as needed: this is NOT checked) and k the number
    of particles.  Cells are frozensets of "token" indices, where tokens are the vertices
    (indices 0..N-1) followed by the edges (indices N..N+M-1).  A cell's degree is the number
    of edge tokens in it.  The boundary of a cell replaces one chosen edge e = (u, v), u < v,
    by each of its endpoints, with sign

        (-1)^j  *  (+1 for u, -1 for v),

    where j is the position of e among the cell's edges in increasing token order.  This is
    the cellular boundary of a product of intervals [u, v] (each oriented so that its
    boundary is u - v; any fixed orientation gives an isomorphic complex), with the
    Koszul sign of the factors in front of it.

    Returns a `gbt.chain.ChainComplex`.  The number of cells grows like C(N + M, k), so this
    is for small cases only.
    """
    N, M = G.N, G.M
    toks = [('v', u) for u in range(N)] + [('e', i) for i in range(M)]

    def support(t):
        """The vertices a token occupies: itself for a vertex, both ends for an edge."""
        return {t[1]} if t[0] == 'v' else set(G.edges[t[1]])

    sup = [support(t) for t in toks]
    tokidx = {t: i for i, t in enumerate(toks)}

    cells = defaultdict(set)

    def rec(start, chosen, used):
        """Depth-first enumeration of k-sets of tokens with pairwise disjoint supports:
        tokens are taken in increasing index order from `start`, `used` is the union of
        the supports chosen so far; each completed set is filed under its degree."""
        if len(chosen) == k:
            cells[sum(1 for i in chosen if toks[i][0] == 'e')].add(frozenset(chosen))
            return
        for i in range(start, len(toks)):
            if sup[i] & used:
                continue
            rec(i + 1, chosen + [i], used | sup[i])

    rec(0, [], set())

    B = {}
    for d in sorted(cells):
        if d == 0:
            continue
        mat = {}
        for c in cells[d]:
            col = {}
            es = sorted(i for i in c if toks[i][0] == 'e')
            for j, ei in enumerate(es):
                u, v = G.edges[toks[ei][1]]
                for w in (u, v):
                    face = frozenset((c - {ei}) | {tokidx[('v', w)]})
                    s = (-1) ** j * (1 if w == u else -1)
                    col[face] = col.get(face, 0) + s
                    if col[face] == 0:
                        del col[face]
            mat[c] = col
        B[d] = mat
    return ChainComplex(cells, B)
