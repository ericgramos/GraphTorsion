"""unordered.py -- the reduced UNORDERED Swiatkowski complex (positive control for torsion detection).

Mathematics.  For the unordered configuration space UConf_n(G) = Conf_n(G)/S_n the reduced Swiatkowski
complex (Swiatkowski 2001; An-Drummond-Cole-Knudsen 2019) has as cells: a set of essential vertices, each
carrying a non-base half-edge h (generator h - h0), together with a number of particles on each edge
(a composition of n - k into E parts).  The differential is the ordered one with labels forgotten:
moving the particle at the i-th occupied vertex into edge(h) or edge(h0) with sign (-1)^i.

It is known that unordered configuration spaces of non-planar graphs have 2-torsion in H_1 (Ko-Park 2012),
e.g. H_1(UConf_2(K_5)) = Z^6 + Z/2.  Recovering these classes shows that the elimination machinery does
detect torsion when it is present.  (The ordered engine with -k "aa...a" computes the same groups; the
two independent implementations agree.)

Usage:  python3 unordered.py      prints UConf_n of K4, K5, K33 for n = 2, 3, 4.
Library: uhomology(edges, n) -> dict(cells, betti, torsion).
"""
import itertools
from gconf import normalize, sparse_snf


def compositions(total, parts):
    """All ordered ways to write `total` as a sum of `parts` non-negative integers."""
    if parts == 0:
        if total == 0:
            yield ()
        return
    for a in range(total + 1):
        for rest in compositions(total - a, parts - 1):
            yield (a,) + rest


def uhomology(edges, n):
    """Integral homology of the unordered configuration space UConf_n via the reduced unordered Swiatkowski complex."""
    V, edges = normalize(edges)
    E = len(edges)
    half = {v: [] for v in range(V)}
    for e, (u, w) in enumerate(edges):
        half[u].append((e, 0)); half[w].append((e, 1))
    base = {v: half[v][0] for v in range(V) if half[v]}
    choices = {v: [h for h in half[v] if h != base[v]] for v in range(V) if half[v]}
    ess = [v for v in range(V) if len(half[v]) >= 2]
    D = min(n, len(ess))
    cells = {}
    for k in range(D + 1):
        lst = []
        for vs in itertools.combinations(ess, k):
            for hs in itertools.product(*[choices[v] for v in vs]):
                for comp in compositions(n - k, E):
                    lst.append((tuple(zip(vs, hs)), comp))
        cells[k] = lst
    index = {k: {c: i for i, c in enumerate(cells[k])} for k in cells}
    rank, divs = {}, {}
    for k in range(1, D + 1):
        trip = []
        for i, (vh, comp) in enumerate(cells[k]):
            for pos, (v, h) in enumerate(vh):
                sgn = (-1) ** pos
                rest = vh[:pos] + vh[pos + 1:]
                for hh, cf in ((h, sgn), (base[v], -sgn)):
                    c2 = list(comp); c2[hh[0]] += 1
                    trip.append((i, index[k - 1][(rest, tuple(c2))], cf))
        rank[k], divs[k] = sparse_snf(len(cells[k]), len(cells[k - 1]), trip, tag="u")
    sizes = [len(cells[k]) for k in range(D + 1)]
    betti = [sizes[k] - rank.get(k, 0) - rank.get(k + 1, 0) for k in range(D + 1)]
    torsion = [divs.get(k + 1, []) for k in range(D + 1)]
    return dict(cells=sizes, betti=betti, torsion=torsion)


if __name__ == "__main__":
    K = lambda m: list(itertools.combinations(range(m), 2))
    Kb = lambda a, b: [(i, a + j) for i in range(a) for j in range(b)]
    for name, g in [("K4 (planar)", K(4)), ("K5", K(5)), ("K33", Kb(3, 3))]:
        for n in (2, 3, 4):
            r = uhomology(g, n)
            print(f"UConf_{n}({name}): betti={r['betti']} torsion={r['torsion']}", flush=True)
