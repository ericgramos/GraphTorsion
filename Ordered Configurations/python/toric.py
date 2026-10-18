"""toric.py -- the "toric bound": the largest homological degree reachable by products of disjoint basic
classes.  A basic class is a cycle class (one particle running around an embedded cycle) or a star class
(two particles exchanging at an essential vertex).  For n particles

    tor(G, n) = max over vertex-disjoint packings P of cycles  of  |P| + min(#essential vertices off P,
                                                                         (n - |P|) // 2),   |P| <= n.

The survey found that the top nonvanishing degree of H_*(Conf_n G) is always >= tor(G, n), and exceeds it
(by exactly one) in 16 cases ("exotic" classes, e.g. H_3(Conf_4(K_5)) = Z^3).  Requires networkx.
"""
import itertools
from collections import Counter
import networkx as nx
from gconf import normalize


def cycle_vertex_sets(edges):
    V, edges = normalize(edges)
    sets = set()
    cnt = Counter(tuple(sorted(e)) for e in edges)
    for (u, v), m in cnt.items():
        if u == v:
            sets.add(frozenset([u]))
        elif m >= 2:
            sets.add(frozenset([u, v]))
    G = nx.Graph()
    G.add_nodes_from(range(V))
    G.add_edges_from((u, v) for u, v in edges if u != v)
    for c in nx.simple_cycles(G, length_bound=V):
        sets.add(frozenset(c))
    # keep inclusion-minimal sets: any packing can use a minimal one instead
    mins = [s for s in sets if not any(t < s for t in sets)]
    return V, edges, mins


def toric_bound(edges, n):
    V, edges, cycles = cycle_vertex_sets(edges)
    deg = Counter()
    for u, v in edges:
        deg[u] += 1; deg[v] += 1
    ess = {v for v in range(V) if deg[v] >= 3}
    best = 0
    def rec(i, used, a):
        nonlocal best
        if a <= n:
            s = len(ess - used)
            best = max(best, a + min(s, (n - a) // 2))
        if a >= n:
            return
        for j in range(i, len(cycles)):
            c = cycles[j]
            if not (c & used):
                rec(j + 1, used | c, a + 1)
    rec(0, frozenset(), 0)
    return best


if __name__ == "__main__":
    # Compare the toric bound with the top nonvanishing degree for every case of the exact survey.
    import json, os
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "survey", "survey_exact.jsonl")
    below = exotic = 0
    seen = set()
    for line in open(fn):
        d = json.loads(line)
        top = max(i for i, b in enumerate(d["betti"]) if b > 0)
        tb = toric_bound(d["edges"], d["n"])
        if top < tb:
            below += 1
            print(f"BELOW BOUND  {d['name']} n={d['n']}: top degree {top}, toric bound {tb}")
        elif top > tb:
            exotic += 1
            print(f"exotic       {d['name']} n={d['n']}: top degree {top} = toric bound {tb} + {top - tb}, "
                  f"H_{top} = Z^{d['betti'][top]}")
    print(f"cases below the bound: {below};  cases above the bound: {exotic}")
