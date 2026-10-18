"""G2_x of the triangle T(a,b,c): remove vertex x, its edges become hairs (leaves)."""
import sys
def T(a,b,c): return [(0,1)]*a+[(0,2)]*b+[(1,2)]*c
def g2(edges, x):
    verts = sorted({v for e in edges for v in e} - {x})
    idx = {v: i for i, v in enumerate(verts)}; nxt = len(verts); out = []
    for (s, t) in edges:
        if s == x: out.append((nxt, idx[t])); nxt += 1
        elif t == x: out.append((idx[s], nxt)); nxt += 1
        else: out.append((idx[s], idx[t]))
    return nxt, out
if __name__ == "__main__":
    a, b, c, x = map(int, sys.argv[1:5])
    V, ed = g2(T(a, b, c), x)
    print(V, len(ed)); [print(s, t) for s, t in ed]
