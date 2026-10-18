"""
Abrams' discretised configuration space -- an independent model, used only to cross-check.

For a graph G, the discretised configuration space UD^k(G) is the cube complex whose cells
are the k-element sets of pairwise disjoint closed cells of G (so each chosen cell is a
vertex or a whole edge, and no two share an endpoint).  The homological degree of a cell is
the number of edges it uses.

Abrams' theorem, as sharpened by Prue--Scrimshaw, says UD^k(G) is a deformation retract of
UConf_k(G) provided G is sufficiently subdivided:

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
"""

from collections import defaultdict

from .chain import ChainComplex


def abrams(G, k):
    """The discretised configuration space UD^k(G) as a ChainComplex.

    Cells are frozensets of "token" indices, where tokens are the vertices followed by the
    edges.  The boundary of a cell replaces one chosen edge by each of its endpoints, with
    the sign (-1)^(position of the edge among the chosen edges) times +1 for the smaller
    endpoint and -1 for the larger.
    """
    N, M = G.N, G.M
    toks = [('v', u) for u in range(N)] + [('e', i) for i in range(M)]

    def support(t):
        return {t[1]} if t[0] == 'v' else set(G.edges[t[1]])

    sup = [support(t) for t in toks]
    tokidx = {t: i for i, t in enumerate(toks)}

    cells = defaultdict(set)

    def rec(start, chosen, used):
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
