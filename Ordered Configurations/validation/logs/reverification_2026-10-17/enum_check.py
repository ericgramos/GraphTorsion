# Independent enumeration of connected multigraphs (loops allowed) with all degrees >= 3, by genus,
# using networkx isomorphism, to check python/graphs.py.
import itertools, sys, networkx as nx
from networkx.algorithms.isomorphism import MultiGraphMatcher
def gen(g):
    out = []
    for Vn in range(1, 2 * g - 1):
        E = Vn + g - 1
        slots = [(i, j) for i in range(Vn) for j in range(i, Vn)]
        def rec(k, left, deg, chosen):
            if left == 0:
                if all(d >= 3 for d in deg):
                    out.append((Vn, list(chosen)))
                return
            if k == len(slots): return
            # remaining capacity check: degree deficit must be coverable
            deficit = sum(max(0, 3 - d) for d in deg)
            if deficit > 2 * left: return
            i, j = slots[k]
            for m in range(left, -1, -1):
                d2 = deg[:]
                if i == j: d2[i] += 2 * m
                else: d2[i] += m; d2[j] += m
                rec(k + 1, left - m, d2, chosen + [(i, j)] * m)
        rec(0, E, [0] * Vn, [])
    # keep connected, dedupe by isomorphism
    reps = {}
    for Vn, e in out:
        G = nx.MultiGraph(); G.add_nodes_from(range(Vn)); G.add_edges_from(e)
        if not nx.is_connected(G): continue
        inv = (Vn, tuple(sorted(dict(G.degree()).values())), sum(1 for u, v in e if u == v))
        lst = reps.setdefault(inv, [])
        if not any(MultiGraphMatcher(G, H).is_isomorphic() for H in lst):
            lst.append(G)
    return [H for l in reps.values() for H in l]
for g in map(int, sys.argv[1:]):
    L = gen(g)
    ess3 = sum(1 for H in L if H.number_of_nodes() >= 3)
    print(f"genus {g}: {len(L)} classes, {ess3} with >= 3 vertices", flush=True)
