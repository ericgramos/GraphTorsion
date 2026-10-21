"""
Printing and summarising homology, so that every experiment reports in the same format.

A homology dict here is always {degree: (free rank, list of invariant factors)}, which is
what all three pipelines return.

EXAMPLE
-------
    >>> from gbt import fast, complete_bip
    >>> from gbt.report import homology_line, torsion_line, has_odd_torsion
    >>> H, _, _ = fast.compute(complete_bip(4, 5), 3)
    >>> homology_line(H)
    'H_1=Z^12 + Z/2  H_2=Z^155 + (Z/2)^3 + Z/6  H_3=0'
    >>> torsion_line(H)
    'H_1: Z/2   H_2: (Z/2)^4 + Z/3'
    >>> has_odd_torsion(H)
    True

Note that `homology_line` prints invariant factors (Z/6) and `torsion_line` prints the
primary decomposition (Z/2 + Z/3); see `gbt.chain.fmt` and `gbt.chain.primary`.
"""

from .chain import fmt, odd_part, primary


def homology_line(H, skip_degree_zero=True):
    """One line of "H_1=... H_2=..." for a homology dict.

    Degree 0 is omitted by default: it is Z for every connected graph and carries no
    information.
    """
    return "  ".join("H_%d=%s" % (d, fmt(*H[d]))
                     for d in sorted(H) if d > 0 or not skip_degree_zero)


def torsion_summary(H):
    """{degree: primary decomposition string} for the degrees that have torsion, e.g.
    {1: 'Z/2', 2: '(Z/2)^4 + Z/3'}; degrees with torsion-free homology are omitted."""
    return {d: primary(H[d][1]) for d in sorted(H) if H[d][1]}


def torsion_line(H):
    """"none", or "H_2: (Z/2)^21" style summary of all the torsion."""
    t = torsion_summary(H)
    if not t:
        return "none"
    return "   ".join("H_%d: %s" % (d, s) for d, s in t.items())


def has_torsion(H, min_degree=0):
    """True if some H_d with d >= min_degree has torsion.  With min_degree=2 this ignores
    the Ko--Park torsion in H_1, which is the usual question in the sweeps."""
    return any(H[d][1] for d in H if d >= min_degree)


def has_odd_torsion(H):
    """True if some invariant factor has an odd prime divisor.

    This is the interesting question for most of the sweeps: 2-torsion is everywhere
    (the Ko--Park class alone supplies it), whereas odd torsion only showed up for
    bipartite graphs with a sparse large side.
    """
    return any(odd_part(x) > 1 for d in H for x in H[d][1])


def exponents(H):
    """The set of invariant factors appearing anywhere in H, as a sorted list (e.g. [2, 6]
    for K_{4,5} at k = 3).  These are invariant factors, not prime powers."""
    return sorted({x for d in H for x in H[d][1]})
