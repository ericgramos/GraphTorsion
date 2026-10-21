"""explode.py -- vertex-explosion rank bookkeeping (rational).

Given Betti numbers of Conf_n(G), Conf_n(G_v) and Conf_{n-1}(G_v) (G_v = G exploded at v), the long exact
sequence of the ordered vertex explosion (the complex of G is the mapping cone of a map Delta between
complexes of G_v; validation/check_explosion.py checks this cell by cell)
    ... -> A_i --delta_i--> B_i -> H_i(Conf_n G) -> A_{i-1} --delta_{i-1}--> B_{i-1} -> ...
with A_i = (d - 1) n b_i(Conf_{n-1} G_v) and B_i = b_i(Conf_n G_v) determines the ranks of the
connecting maps: r_0 = B_0 - b_0(G),  r_i = B_i + A_{i-1} - r_{i-1} - b_i(G).
analyse(edges, v, n) prints these ranks with kernel and cokernel dimensions.  It was used to see how
classes of G arise from G_v (cokernel part) and from the boundary (kernel part).
"""
import itertools, json, sys
from collections import Counter
from gconf import homology, normalize

_cache = {}


def betti(edges, n):
    """Betti numbers of Conf_n of the graph (rational, through gconf.homology), cached by (edges, n)."""
    key = (tuple(sorted(map(tuple, edges))), n)
    if key not in _cache:
        if n == 0:
            _cache[key] = [1]
        else:
            r = homology(edges, n, tag="ex")
            assert not any(r["torsion"])
            _cache[key] = r["betti"]
    return _cache[key]


def explode(edges, v):
    """The exploded graph Gamma_v: every half-edge at v is replaced by a new leaf (a loop at v gives two leaves)."""
    V, edges = normalize(edges)
    out, nxt = [], V
    for (a, b) in edges:
        if a == v and b == v:        # loop: both ends become separate leaves
            out.append((nxt, nxt + 1)); nxt += 2
        elif a == v:
            out.append((nxt, b)); nxt += 1
        elif b == v:
            out.append((a, nxt)); nxt += 1
        else:
            out.append((a, b))
    return out


def analyse(edges, v, n):
    """Rank bookkeeping of the vertex explosion at v for Conf_n: compares the Betti numbers of Gamma and Gamma_v."""
    V, edges = normalize(edges)
    d = sum((a == v) + (b == v) for a, b in edges)
    Gv = explode(edges, v)
    bG = betti(edges, n)
    B = betti(Gv, n)
    Ap = betti(Gv, n - 1)
    D = max(len(bG), len(B), len(Ap) + 1)
    get = lambda L, i: L[i] if 0 <= i < len(L) else 0
    A = [(d - 1) * n * get(Ap, i) for i in range(D)]
    r = []
    for i in range(D):
        ri = get(B, i) + (A[i - 1] if i else 0) - (r[i - 1] if i else 0) - get(bG, i)
        r.append(ri)
    rows = []
    for i in range(D):
        rows.append(dict(i=i, A=A[i], B=get(B, i), rank=r[i], coker=get(B, i) - r[i], ker=A[i] - r[i], b=get(bG, i)))
    return rows


if __name__ == "__main__":
    K = lambda m: list(itertools.combinations(range(m), 2))
    Kb = lambda a, b: [(i, a + j) for i in range(a) for j in range(b)]
    pet = [(i, (i + 1) % 5) for i in range(5)] + [(i, i + 5) for i in range(5)] + [(5 + i, 5 + (i + 2) % 5) for i in range(5)]
    graphs = {"theta3": [(0, 1)] * 3, "K4": K(4), "K5": K(5), "K33": Kb(3, 3), "Petersen": pet,
              "prism": [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)]}
    maxn = {"theta3": 6, "K4": 5, "K5": 4, "K33": 4, "Petersen": 4, "prism": 4}
    for name, g in graphs.items():
        for n in range(1, maxn[name] + 1):
            rows = analyse(g, 0, n)
            s = "  ".join(f"i={x['i']}: rk {x['rank']}/{x['A']}->{x['B']} (ker {x['ker']}, coker {x['coker']})" for x in rows)
            print(f"{name} n={n}: {s}", flush=True)
