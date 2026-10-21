"""walks.py -- split a signed 1-chain with coefficients +-1 (a 1-cycle of a configuration space) into closed walks.

A 1-cell "particle p at vertex v with half-edge h" with coefficient c is an arrow from the 0-cell with p at the end of
the base edge of v to the 0-cell with p at the end of edge(h) if c = +1, and the reverse arrow if c = -1.  A 1-cycle is
then a balanced directed multigraph and splits into closed walks (Hierholzer).  Each walk is printed as a sequence of
moves "particle p passes through vertex v from edge e to edge e'".
"""
from collections import defaultdict


def closed_walks(chain, base, edges):
    """Split a balanced +-1 chain into closed walks; returns the walks and their text."""
    def push(words, p, h):
        """Put particle p into edge(h) at the end of h."""
        e = h // 2
        w = list(words)
        w[e] = ((p,) + w[e]) if h % 2 == 0 else (w[e] + (p,))
        return tuple(w)
    out = defaultdict(list)
    for (v, p, h, words), c in chain.items():
        if c == 0:
            continue
        assert abs(c) == 1, "coefficients must be +-1"
        a, b = push(words, p, base[v]), push(words, p, h)
        src, dst, s = (a, b, +1) if c > 0 else (b, a, -1)
        out[src].append((dst, v, p, h, s))
    walks = []
    while any(out.values()):
        start = next(k for k, v in out.items() if v)
        path, cur = [], start
        while out[cur]:
            dst, v, p, h, s = out[cur].pop()
            path.append((cur, dst, v, p, h, s))
            cur = dst
            if cur == start and not out[cur]:
                break
        walks.append(path)
    text = []
    for path in walks:
        conf = lambda w: " ".join(f"{edges[e]}#{e}:{list(x)}" for e, x in enumerate(w) if x)
        text.append(f"closed walk of {len(path)} moves, starting from {conf(path[0][0])}")
        for cur, dst, v, p, h, s in path:
            eb, eh = base[v] // 2, h // 2
            fr, to = (eb, eh) if s > 0 else (eh, eb)
            text.append(f"    {p} through vertex {v}: edge {edges[fr]}#{fr} -> edge {edges[to]}#{to}")
    return walks, text
