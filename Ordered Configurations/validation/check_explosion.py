"""Chain-level check of the ordered vertex explosion.

Independent implementation of the reduced ordered Swiatkowski complex (graphs may have
leaves, loops, multi-edges). For a vertex v we verify, cell by cell and with signs, that

    S_n(Gamma)  ==  Cone( Delta : (+)_{q, h != h0} S_{[n]-q}(Gamma_v)  ->  S_n(Gamma_v) )

under Phi(q, h, c') = eps(c) * c, where c is c' with particle q placed at v with generator
h - h0, and eps(c) = (-1)^{# occupied vertices of c preceding v}.  Cone convention:
d(b, a) = (d b + Delta a, -d a), i.e.  d(q,h,c') = Delta(c') - (q,h,d c').
Also checks that each sigma_h^q is a chain map S_{[n]-q}(Gamma_v) -> S_n(Gamma_v).
"""
import itertools, sys
from collections import defaultdict


class Complex:
    def __init__(self, V, edges, labels):
        """The reduced ordered Swiatkowski complex of a graph for the given particle labels (an independent implementation)."""
        self.V, self.edges, self.labels = V, list(edges), tuple(sorted(labels))
        self.E = len(edges)
        self.halves = defaultdict(list)                 # vertex -> half-edges (2e = tail end, 2e+1 = head end)
        for e, (a, b) in enumerate(edges):
            self.halves[a].append(2 * e); self.halves[b].append(2 * e + 1)
        self.base = {v: min(hs) for v, hs in self.halves.items()}

    def cells(self):
        """all cells: (verts, words) with verts[v] = None or (particle, half-edge != base)."""
        L, V = self.labels, self.V
        occ_choices = [[None] + [(h,) for h in self.halves.get(v, []) if h != self.base.get(v)] for v in range(V)]
        out = []
        for occ in itertools.product(*occ_choices):
            occupied = [v for v in range(V) if occ[v] is not None]
            if len(occupied) > len(L): continue
            for parts in itertools.permutations(L, len(occupied)):
                verts = [None] * V
                for v, p in zip(occupied, parts): verts[v] = (p, occ[v][0])
                rest = [p for p in L if p not in parts]
                for perm in itertools.permutations(rest):
                    for cuts in itertools.combinations(range(len(rest) + self.E - 1), self.E - 1):
                        words, prev, i = [], -1, 0
                        for cpos in list(cuts) + [len(rest) + self.E - 1]:
                            m = cpos - prev - 1; words.append(tuple(perm[i:i + m])); i += m; prev = cpos
                        out.append((tuple(verts), tuple(words)))
        return out

    def push(self, cell, v, h):
        """Remove the particle at vertex v and push it into edge(h) at the end of h (prepend at a tail end, append at a head end)."""
        verts, words = cell
        p = verts[v][0]
        nv = list(verts); nv[v] = None
        nw = list(words); e = h // 2
        nw[e] = (p,) + words[e] if h % 2 == 0 else words[e] + (p,)
        return (tuple(nv), tuple(nw))

    def d(self, cell):
        """Boundary of a cell as a dictionary {face: coefficient}."""
        verts, _ = cell
        out = defaultdict(int); j = 0
        for v in range(self.V):
            if verts[v] is None: continue
            s = 1 if j % 2 == 0 else -1; j += 1
            out[self.push(cell, v, verts[v][1])] += s
            out[self.push(cell, v, self.base[v])] -= s
        return {k: c for k, c in out.items() if c}


def explode(V, edges, v):
    """Gamma_v: each half-edge at v gets its own new leaf; edge indices and orientations kept."""
    nxt = V; out = []
    for (a, b) in edges:
        if a == v: a = nxt; nxt += 1
        if b == v: b = nxt; nxt += 1
        out.append((a, b))
    return nxt, out


def check(V, edges, v, n):
    """Check the ordered vertex explosion for (Gamma, v, n) cell by cell and sign by sign; returns the number of failures."""
    G = Complex(V, edges, range(n))
    Vx, ex = explode(V, edges, v)
    hv = [h for h in G.halves[v] if h != G.base[v]]
    h0 = G.base[v]
    # sigma_h^q on a Gamma_v cell: insert q at the v-end of edge(h); the Gamma_v vertex set is 0..Vx-1
    def sigma(cell, q, h):
        """sigma_h^q: insert particle q at the v-end of edge(h) (a cell of Gamma_v becomes a cell of Gamma)."""
        verts, words = cell
        nw = list(words); e = h // 2
        nw[e] = (q,) + words[e] if h % 2 == 0 else words[e] + (q,)
        return (verts, tuple(nw))
    def lift(cellx):   # Gamma_v cell -> Gamma cell with v empty (drop the extra leaf slots)
        """A cell of Gamma_v as a cell of Gamma with v empty (drop the extra leaf slots)."""
        verts, words = cellx
        return (tuple(verts[:V]), words)
    nbad = 0; ncells = 0
    # (i) subcomplex: Gamma_v complex equals the v-empty part of Gamma's complex
    Gx = Complex(Vx, ex, range(n))
    for cx in Gx.cells():
        ncells += 1
        lhs = {lift(k): c for k, c in Gx.d(cx).items()}
        rhs = G.d(lift(cx))
        if lhs != rhs: nbad += 1
    # (ii)+(iii): v-occupied cells, cone formula
    for q in range(n):
        Gq = Complex(Vx, ex, [p for p in range(n) if p != q])
        for cx in Gq.cells():
            # sigma is a chain map: d sigma = sigma d
            for h in hv + [h0]:
                a = {k: c for k, c in Gx.d(sigma(cx, q, h)).items()}
                b = defaultdict(int)
                for k, c in Gq.d(cx).items(): b[sigma(k, q, h)] += c
                b = {k: c for k, c in b.items() if c}
                if a != b: nbad += 1
            for h in hv:
                verts, words = lift(cx)
                nv = list(verts); nv[v] = (q, h); c = (tuple(nv), words)
                eps = (-1) ** sum(1 for w in range(v) if nv[w] is not None)
                # d_Gamma(Phi(q,h,c')) where Phi(q,h,c') = eps * c
                lhs = defaultdict(int)
                for k, coef in G.d(c).items(): lhs[k] += eps * coef
                # Phi(cone differential): Delta(c') - (q,h, d c')
                rhs = defaultdict(int)
                rhs[lift(sigma(cx, q, h))] += 1
                rhs[lift(sigma(cx, q, h0))] -= 1
                for k, coef in Gq.d(cx).items():
                    kv, kw = lift(k); nk = list(kv); nk[v] = (q, h); kc = (tuple(nk), kw)
                    epsk = (-1) ** sum(1 for w in range(v) if nk[w] is not None)
                    rhs[kc] += -coef * epsk
                lhs = {k: c for k, c in lhs.items() if c}; rhs = {k: c for k, c in rhs.items() if c}
                ncells += 1
                if lhs != rhs: nbad += 1
    # (iv) cell counts: |S(Gamma)| = |S(Gamma_v)| + (d-1) * n * |S(Gamma_v; n-1)|
    tot = len(G.cells())
    pred = len(Gx.cells()) + len(hv) * n * len(Complex(Vx, ex, range(n - 1)).cells())
    return nbad, ncells, tot == pred


if __name__ == "__main__":
    tests = [
        ("Y graph, center", 4, [(0, 1), (0, 2), (0, 3)], 0, 3),
        ("theta_3, vertex 1", 2, [(0, 1), (0, 1), (0, 1)], 1, 3),
        ("T(1,2,2), middle vertex 1", 3, [(0, 1), (0, 2), (0, 2), (1, 2), (1, 2)], 1, 3),
        ("T(2,2,2), last vertex 2", 3, [(0, 1)] * 2 + [(0, 2)] * 2 + [(1, 2)] * 2, 2, 3),
        ("K4, vertex 2", 4, list(itertools.combinations(range(4), 2)), 2, 3),
        ("loop + bridge (lollipop), loop vertex", 2, [(0, 1), (1, 1)], 1, 3),
        ("figure eight, the vertex", 1, [(0, 0), (0, 0)], 0, 3),
        ("loop and multi-edge, vertex 0", 2, [(0, 0), (0, 1), (1, 0), (0, 1)], 0, 3),
        ("tripled C4, vertex 3", 4, [(0, 1), (1, 2), (2, 3), (3, 0)] * 3, 3, 2),
        ("K4 + leaf at v, vertex 0", 5, list(itertools.combinations(range(4), 2)) + [(0, 4)], 0, 3),
    ]
    allok = True
    for name, V, ed, v, n in tests:
        nbad, nc, cnt = check(V, ed, v, n)
        ok = nbad == 0 and cnt
        allok &= ok
        print(f"{name:42s} n={n}: {nc} identities checked, mismatches {nbad}, cell count {'ok' if cnt else 'WRONG'}  {'OK' if ok else 'FAIL'}", flush=True)
    print("ALL OK" if allok else "FAILURES")
