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
any of its deg(v) incident edges) instead of deg(v) - 1.  (Here j counts OCCUPIED vertices,
sitting or moving, so the sum is over occupation patterns, not over homological degree.)

WHERE THE PARKED FACTOR COMES FROM
----------------------------------
The k - j particles that are not at a vertex are parked in edge interiors, and only the
number on each edge matters: a monomial of degree k - j in M variables, of which there are
C(M + k - j - 1, k - j) (stars and bars).

EXAMPLE
-------
    >>> from gbt import complete, reduced_cell_count, cell_count
    >>> reduced_cell_count(complete(7), 3), cell_count(complete(7), 3)
    (25256, 46704)
    >>> from gbt.sizes import reduced_cell_count_by_degree
    >>> reduced_cell_count_by_degree(complete(7), 3)
    {0: 1771, 1: 8085, 2: 11025, 3: 4375}

These are the numbers to compare against the pipeline brackets in the README (about
6*10^5 cells for `gbt.fast`, 2*10^6 for `gbt.lean`, more for a single torsion degree).
"""

from math import comb


def _elementary_symmetric(values):
    """[e_0, e_1, ..., e_n] of the multiset `values`, computed by the standard recursion.

    e_j is the j-th elementary symmetric polynomial evaluated at `values`: the sum over all
    j-element sub-multisets of the product of their entries.  The recursion multiplies the
    generating polynomial prod (1 + x t) one factor at a time, so the cost is O(n^2).
    """
    e = [1]
    for x in values:
        e = [e[0]] + [e[i] + x * e[i - 1] for i in range(1, len(e))] + [x * e[-1]]
    return e


def reduced_cell_count(G, k):
    """Number of cells (all degrees together) of the fully reduced weight-k Swiatkowski
    complex of G, by the closed formula in the module docstring.

    G is a `Graph` and k the number of particles.  Costs microseconds, so sweeps call it
    before building anything, to sort jobs and to skip those above a memory cap.  This is
    the total of `gbt.fast.compute(G, k)[1]` (the raw sizes) over all degrees.
    """
    e = _elementary_symmetric([G.deg(v) - 1 for v in range(G.N)])
    return sum(e[j] * comb(G.M + (k - j) - 1, k - j)
               for j in range(0, min(k, G.N) + 1) if j < len(e))


def cell_count(G, k):
    """Number of cells (all degrees together) of the UNREDUCED weight-k Swiatkowski complex
    of G, i.e. sum(swiatkowski(G, k).sizes().values()).  Each occupied vertex contributes
    deg(v) + 1 states (sitting, or moving along one of its edges)."""
    e = _elementary_symmetric([G.deg(v) + 1 for v in range(G.N)])
    return sum(e[j] * comb(G.M + (k - j) - 1, k - j)
               for j in range(0, min(k, G.N) + 1) if j < len(e))


def reduced_cell_count_by_degree(G, k):
    """{degree: number of cells} for the fully reduced complex.

    In the fully reduced complex every occupied vertex is moving, so the homological degree
    j equals the number of occupied vertices and the degree-j count is the j-th term of the
    formula in the module docstring.  `gbt.lean.torsion_in_degree(G, k, d)` builds exactly
    the slices d and d+1, so their sum here is its memory bill.

    The degree-j count needs the elementary symmetric polynomial restricted to subsets of
    size j with every factor nonzero, which is what e_j already is: a vertex of degree 1
    contributes a factor 0 and so kills exactly the subsets containing it.
    """
    e = _elementary_symmetric([G.deg(v) - 1 for v in range(G.N)])
    return {j: e[j] * comb(G.M + (k - j) - 1, k - j)
            for j in range(0, min(k, G.N) + 1) if j < len(e)}
