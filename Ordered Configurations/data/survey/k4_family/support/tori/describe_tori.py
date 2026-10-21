"""describe_tori.py -- read a torus chain back as two closed walks of configurations.

A torus of analysis/decomp is a product z_x * z_P of a 1-cycle z_x of Conf_A(star of x) and a 1-cycle z_P of
Conf_B(Gamma minus x).  Each of its 2-cells has the vertex x occupied (by a particle of A) and one other vertex v
occupied (by a particle of B).  This script splits every cell into its two factors (the particles of A are the star
factor, those of B the piece factor; on an edge at x the particles of A sit next to x), collects the two 1-chains,
and walks each one as a closed sequence of 0-cells.  A 1-cell "particle p at vertex v with half-edge h" joins the
0-cell where p sits at the end of the base edge of v to the 0-cell where p sits at the end of edge(h), so a walk
reads as "p passes through v from the base edge into edge(h)" (or back).

Usage: python3 describe_tori.py tori_file graph_file id [id ...]
       tori_file: the DUMP output of analysis/decomp ("Z id 0 x=.. A=.. i1=.. i2=.." + "coef cell");
       graph_file: "V E n" + edges (the subgraph H the tori live on).
"""
import re
import sys
from collections import defaultdict

CELL = re.compile(r"V\[(.*)\] E\[(.*)\]$")


def parse(cell, V, E):
    """Parse the text form of a cell (edge words padded to E slots)."""
    m = CELL.match(cell)
    vs = m.group(1).split(",")
    es = m.group(2).split("|")
    verts = [None if x == "-" else tuple(map(int, x.split("/"))) for x in vs]
    eds = [tuple(map(int, x.split())) if x.strip() else () for x in es] + [()] * (E - len(es))
    return verts, eds


def load_graph(fn):
    """Read 'V E n' and the edges from a graph file."""
    nums = list(map(int, open(fn).read().split()))
    V, E = nums[0], nums[1]
    edges = [(nums[3 + 2 * i], nums[4 + 2 * i]) for i in range(E)]
    return V, E, edges


def walk(chain, base, edges, label):
    """chain: {1-cell (vertex, particle, h, words): coef}.  Returns the cyclic walk as text lines."""
    def push(words, p, h):
        """Put particle p into edge(h) at the end of h (prepend at a tail end, append at a head end)."""
        e = h // 2
        w = list(words)
        w[e] = ((p,) + w[e]) if h % 2 == 0 else (w[e] + (p,))
        return tuple(w)
    adj = defaultdict(list)
    for (v, p, h, words), c in chain.items():
        a = push(words, p, base[v])     # p on the base edge
        b = push(words, p, h)           # p on edge(h)
        adj[a].append((b, v, p, h, c, +1))     # traversed from a to b: p goes from the base edge into edge(h)
        adj[b].append((a, v, p, h, c, -1))
    if not adj:
        return [f"{label}: empty"]
    start = min(adj)
    lines, cur, used, signs = [], start, set(), set()
    while True:
        nxt = [t for t in adj[cur] if (t[1], t[2], t[3], t[0] if t[5] < 0 else cur) not in used]
        if not nxt:
            break
        b, v, p, h, c, d = nxt[0]
        used.add((v, p, h, b if d < 0 else cur))
        signs.add(c * d)                 # along a cycle, coefficient times direction is constant
        e_base, e_h = base[v] // 2, h // 2
        src, dst = (e_base, e_h) if d > 0 else (e_h, e_base)
        lines.append(f"particle {p} through vertex {v}, from edge {edges[src]}#{src} to edge {edges[dst]}#{dst}")
        cur = b
        if cur == start:
            break
    def conf(words):
        """Short text form of a 0-cell: the nonempty edge words."""
        return " ".join(f"#{e}:{list(w)}" for e, w in enumerate(words) if w)
    ok = cur == start and len(signs) == 1 and len(lines) == len(chain)
    return [f"{label}: start at {conf(start)}"] + ["  " + l for l in lines] + [f"  ({len(lines)} moves; a single closed walk with consistent signs: {ok})"]


def main(tfn, gfn, ids):
    """Print both factors of the requested tori as closed walks."""
    V, E, edges = load_graph(gfn)
    halves = defaultdict(list)
    for e, (a, b) in enumerate(edges):
        halves[a].append(2 * e); halves[b].append(2 * e + 1)
    base = {v: min(hs) for v, hs in halves.items()}
    want = set(ids)
    cur, hdr, cells = None, {}, defaultdict(list)
    for line in open(tfn):
        if line.startswith("Z "):
            cur = line.split()[1]
            if cur in want:
                hdr[cur] = line.strip()
            continue
        if cur in want and line.strip():
            c, cell = line.rstrip("\n").split(" ", 1)
            cells[cur].append((int(c), cell))
    for i in ids:
        h = hdr[i]
        x = int(re.search(r"x=(\d+)", h).group(1))
        A = set(map(int, re.search(r"A=([\d,]+)", h).group(1).split(",")))
        print(h)
        star, piece = factor_signs(cells[i], V, E, x, A)
        for l in walk(star, base, edges, f"  star factor at vertex {x}, particles {sorted(A)}"):
            print("   " + l)
        B = sorted(set(range(5)) - A)
        for l in walk(piece, base, edges, f"  factor on Gamma minus {x}, particles {B}"):
            print("   " + l)


def factor_signs(cells, V, E, x, A):
    """Recover the two signed 1-chains from the product chain (coefficient of (a, b) = s * c_a * d_b)."""
    prod = {}
    for c, cell in cells:
        verts, eds = parse(cell, V, E)
        wA = tuple(tuple(q for q in w if q in A) for w in eds)
        wB = tuple(tuple(q for q in w if q not in A) for w in eds)
        v = next(u for u in range(V) if u != x and verts[u] is not None)
        a = (x,) + verts[x] + (wA,)
        b = (v,) + verts[v] + (wB,)
        prod[(a, b)] = c
    As = sorted({a for a, _ in prod}); Bs = sorted({b for _, b in prod})
    a0, b0 = As[0], Bs[0]
    star = {a: prod[(a, b0)] * prod[(a0, b0)] for a in As if (a, b0) in prod}
    # the product carries the Koszul sign +1 if the piece vertex is below x, else -1 (analysis/decomp)
    piece = {b: prod[(a0, b)] * (1 if b[0] < x else -1) for b in Bs if (a0, b) in prod}
    star = {(a[0], a[1], a[2], a[3]): s for a, s in star.items()}
    piece = {(b[0], b[1], b[2], b[3]): s for b, s in piece.items()}
    return star, piece


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
