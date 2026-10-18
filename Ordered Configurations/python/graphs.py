"""graphs.py -- enumeration of connected multigraphs with all degrees >= 3, by genus.

Connected graphs without leaves and without vertices of degree 2, considered up to homeomorphism, are the
same as connected multigraphs (loops and multi-edges allowed) with all degrees >= 3, up to isomorphism.
For first Betti number (genus) g >= 2 such a graph has V <= 2g - 2 vertices and E = V + g - 1 edges.
multigraphs(g) enumerates them by backtracking over edge multiplicities and removes isomorphic copies with
a canonical form (canon) computed by brute force over vertex permutations that respect an invariant
colouring (degree, number of loops, neighbour degrees).

Usage:  python3 graphs.py [G]  prints, for g = 2..G (default 4), the number of classes and the number of simple
        ones (3, 15, 111 classes for g = 2, 3, 4; checked against an independent networkx enumeration).  G = 5
        works but takes a long time.
Used by surveys/search.py for the survey of all genus 3 and 4 graphs (10 + 103 classes with at least
three essential vertices).
"""
import itertools, sys
from functools import lru_cache


def _perms_by_class(V, edges):
    """permutations respecting an invariant vertex colouring (degree, loops, multiset of
    neighbour degrees); vertices are mapped to positions sorted by colour."""
    deg = [0] * V; loops = [0] * V
    for u, v in edges:
        deg[u] += 1; deg[v] += 1
        if u == v:
            loops[u] += 1
    nb = [[] for _ in range(V)]
    for u, v in edges:
        nb[u].append(deg[v]); nb[v].append(deg[u])
    nb = [sorted(x) for x in nb]
    col = [(deg[x], loops[x], tuple(nb[x])) for x in range(V)]
    classes = {}
    for x in range(V):
        classes.setdefault(col[x], []).append(x)
    keys = sorted(classes)
    slots = []
    pos = 0
    for k in keys:
        slots.append((classes[k], list(range(pos, pos + len(classes[k])))))
        pos += len(classes[k])
    def rec(i, perm):
        if i == len(slots):
            yield perm
            return
        verts, targets = slots[i]
        for p in itertools.permutations(targets):
            for v, t in zip(verts, p):
                perm[v] = t
            yield from rec(i + 1, perm)
    yield from rec(0, [0] * V)


def canon(V, edges):
    best = None
    for perm in _perms_by_class(V, edges):
        key = tuple(sorted(tuple(sorted((perm[u], perm[v]))) for u, v in edges))
        if best is None or key < best:
            best = key
    return best


def connected(V, edges):
    adj = {i: set() for i in range(V)}
    for u, v in edges:
        adj[u].add(v); adj[v].add(u)
    seen = {0}; st = [0]
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y); st.append(y)
    return len(seen) == V


def multigraphs(g, min_deg=3, max_V=None):
    out = {}
    Vmax = 2 * g - 2 if max_V is None else max_V
    for V in range(1, Vmax + 1):
        E = V + g - 1
        pairs = [(i, j) for i in range(V) for j in range(i, V)]
        # backtracking over multiplicities of each pair
        deg = [0] * V
        mult = [0] * len(pairs)

        def rec(k, remaining):
            if remaining == 0:
                if min(deg) >= min_deg:
                    edges = [p for p, m in zip(pairs, mult) for _ in range(m)]
                    if connected(V, edges):
                        c = canon(V, edges)
                        out.setdefault(c, (V, list(c)))
                return
            if k == len(pairs):
                return
            i, j = pairs[k]
            # vertex i is finished after its last pair (i, V-1); prune if it can't reach min_deg
            for m in range(remaining, -1, -1):
                add = 2 * m if i == j else m
                deg[i] += add if i == j else m
                if i != j:
                    deg[j] += m
                mult[k] = m
                ok = True
                if j == V - 1 and deg[i] < min_deg:
                    ok = False
                if ok:
                    rec(k + 1, remaining - m)
                deg[i] -= add if i == j else m
                if i != j:
                    deg[j] -= m
                mult[k] = 0

        rec(0, E)
    return list(out.values())


if __name__ == "__main__":
    for g in range(2, int(sys.argv[1]) + 1 if len(sys.argv) > 1 else 5):   # genus 5 takes a long time
        gs = multigraphs(g)
        print(g, len(gs), sum(1 for V, e in gs if all(u != v for u, v in e)))
