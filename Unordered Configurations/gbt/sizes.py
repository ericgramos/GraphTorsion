"""
Exact cell counts, for ordering and budgeting jobs before running them.

A sweep has to know how big a complex will be BEFORE building it, both to sort the jobs
from cheapest to most expensive and to skip the ones that will not fit in memory.  Both
counts below are closed forms, so they cost microseconds.

FULLY REDUCED COMPLEX
---------------------
A degree-j cell of the fully reduced weight-k complex is: a choice of j occupied vertices;
for each, one of its deg(v) - 1 difference generators h~_{v,e}; and a monomial of weight
k-j in the M edges.  So

    reduced_cell_count(G, k)
        = sum_{j=0}^{min(k,N)}  e_j( deg(v) - 1 : v in V )  *  C(M + k - j - 1, k - j)

where e_j is the j-th elementary symmetric polynomial.  Specialising:

    K_n       : e_j = C(n, j) (n-2)^j
    K_{a,b}   : e_j = sum_{i+l=j} C(a,i)(b-1)^i C(b,l)(a-1)^l

which are the two formulas the sweeps are ordered by.  Using the symmetric
polynomial directly means one formula covers every graph, including the irregular ones in
the small-graph and Petersen-family sweeps.

UNREDUCED COMPLEX
-----------------
Same, with each occupied vertex contributing deg(v) + 1 states (sitting, or moving along
any of its deg(v) incident edges) instead of deg(v) - 1.
"""

from math import comb


def _elementary_symmetric(values):
    """[e_0, e_1, ..., e_n] of the multiset `values`, computed by the standard recursion."""
    e = [1]
    for x in values:
        e = [e[0]] + [e[i] + x * e[i - 1] for i in range(1, len(e))] + [x * e[-1]]
    return e


def reduced_cell_count(G, k):
    """Number of cells of the fully reduced weight-k Swiatkowski complex of G."""
    e = _elementary_symmetric([G.deg(v) - 1 for v in range(G.N)])
    return sum(e[j] * comb(G.M + (k - j) - 1, k - j)
               for j in range(0, min(k, G.N) + 1) if j < len(e))


def cell_count(G, k):
    """Number of cells of the UNREDUCED weight-k Swiatkowski complex of G."""
    e = _elementary_symmetric([G.deg(v) + 1 for v in range(G.N)])
    return sum(e[j] * comb(G.M + (k - j) - 1, k - j)
               for j in range(0, min(k, G.N) + 1) if j < len(e))


def reduced_cell_count_by_degree(G, k):
    """{degree: number of cells} for the fully reduced complex.

    The degree-j count needs the elementary symmetric polynomial restricted to subsets of
    size j with every factor nonzero, which is what e_j already is: a vertex of degree 1
    contributes a factor 0 and so kills exactly the subsets containing it.
    """
    e = _elementary_symmetric([G.deg(v) - 1 for v in range(G.N)])
    return {j: e[j] * comb(G.M + (k - j) - 1, k - j)
            for j in range(0, min(k, G.N) + 1) if j < len(e)}
