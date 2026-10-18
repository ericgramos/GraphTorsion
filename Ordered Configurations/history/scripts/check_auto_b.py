"""Independent check of oconf19's coinvariant complexes for label x graph-automorphism groups.
Builds the full ordered reduced Swiatkowski complex (with the same invariant base half-edges), the signed
permutation action of each group element, and computes dim H_k(C^G; F_p) via the averaging idempotent
(p prime to |G|). Compares with dim H_k(C_G; F_p) from oconf19 (stage-0 ranks)."""
import itertools, json, subprocess, sys
import flint
from check_explosion import Complex
S = "/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad"

def closure(gens, n):
    idp = tuple(range(n)); G = {idp}; fr = [idp]
    while fr:
        nf = []
        for x in fr:
            for g in gens:
                y = tuple(g[x[t]] for t in range(n))
                if y not in G: G.add(y); nf.append(y)
        fr = nf
    return sorted(G)

def run(V, edges, n, p, labgens, autos, bases):
    C = Complex(V, edges, range(n)); C.base = dict(bases)
    E = len(edges)
    Lall = closure([tuple(g) + tuple(range(len(g), n)) for g in labgens], n) if labgens else [tuple(range(n))]
    def comb(vm, em):
        out = list(vm)
        for x in em:
            f = isinstance(x, str); e = int(x[1:]) if f else x
            out += [V + 2 * e + f, V + 2 * e + 1 - f]
        return tuple(out)
    Aall = closure([comb(vm, em) for vm, em in autos], V + 2 * E) if autos else [tuple(range(V + 2 * E))]
    def act(L, A, cell):
        verts, words = cell
        nv = [None] * V; nw = [None] * E
        occ = [v for v in range(V) if verts[v] is not None]
        for v in occ:
            q, h = verts[v]; nv[A[v]] = (L[q], A[V + h] - V)
        for e in range(E):
            x = A[V + 2 * e] - V; w = tuple(L[q] for q in words[e])
            nw[x // 2] = w[::-1] if x % 2 else w
        w = [A[v] for v in occ]; inv = sum(1 for i in range(len(w)) for j in range(i + 1, len(w)) if w[i] > w[j])
        return (tuple(nv), tuple(nw)), (-1) ** inv
    cells = C.cells(); bydeg = {}
    for c in cells: bydeg.setdefault(sum(1 for x in c[0] if x is not None), []).append(c)
    D = max(bydeg); idx = {k: {c: i for i, c in enumerate(bydeg[k])} for k in bydeg}
    G = [(L, A) for L in Lall for A in Aall]; inv_order = pow(len(G), -1, p)
    # chain-map check d g = g d, and build averaging matrices
    Emat, Dmat = {}, {}
    bad = 0
    for k in range(D + 1):
        m = len(bydeg[k]); M = flint.nmod_mat(m, m, p)
        for j, c in enumerate(bydeg[k]):
            for L, A in G:
                gc, s = act(L, A, c); M[idx[k][gc], j] += s * inv_order
                if k >= 1:   # d(g c) vs g(d c)
                    lhs = {kk: s * vv for kk, vv in C.d(gc).items()}
                    rhs = {}
                    for f, cf in C.d(c).items():
                        gf, sf = act(L, A, f); rhs[gf] = rhs.get(gf, 0) + cf * sf
                    if {a: b for a, b in lhs.items() if b} != {a: b for a, b in rhs.items() if b}: bad += 1
        Emat[k] = M
        if k >= 1:
            Dk = flint.nmod_mat(len(bydeg[k - 1]), m, p)
            for j, c in enumerate(bydeg[k]):
                for f, cf in C.d(c).items(): Dk[idx[k - 1][f], j] += cf
            Dmat[k] = Dk
    dimInv = [Emat[k].rank() for k in range(D + 1)]
    rk = [0] + [(Dmat[k] * Emat[k]).rank() for k in range(1, D + 1)] + [0]
    return [dimInv[k] - rk[k] - rk[k + 1] for k in range(D + 1)], bad, len(G)

def oconf(V, edges, n, p, labgens, autos):
    inp = f"{V} {len(edges)} {n}\n" + "".join(f"{a} {b}\n" for a, b in edges)
    args = ["./oconf22b", "-F", "-q", str(p), "2", "-o", S + "/ca"]
    if labgens: args += ["-g", ";".join(",".join(map(str, g)) for g in labgens)]
    if autos: args += ["-a", ";".join(",".join(map(str, vm)) + "|" + ",".join(map(str, em)) for vm, em in autos)]
    out = subprocess.run(args, input=inp, capture_output=True, text=True)
    if out.returncode: raise RuntimeError(out.stderr)
    d = json.loads(out.stdout); c = d["cells"]
    bases = [int(x) for x in out.stderr.split("bases:")[1].split("\n")[0].split()]
    r = [0] + [m["stages"][0] for m in d["maps"]] + [0]
    tors = {m["k"] - 1: m["stages"][1:] for m in d["maps"] if len(m["stages"]) > 1}
    return [c[k] - r[k] - r[k + 1] for k in range(len(c))], bases, tors

def em_from(V, edges, vm):
    used = set(); em = []
    for (a, b) in edges:
        for j, (c, d) in enumerate(edges):
            if j not in used and (c, d) == (vm[a], vm[b]): em.append(j); used.add(j); break
        else:
            for j, (c, d) in enumerate(edges):
                if j not in used and (d, c) == (vm[a], vm[b]): em.append("~%d" % j); used.add(j); break
            else: raise ValueError("not an automorphism")
    return em

tests = []
K4 = [(0, 1), (0, 2), (0, 3), (1, 2), (2, 3), (3, 1)]
k1 = [1, 0, 3, 2]; k2 = [2, 3, 0, 1]
tests.append(("K4 Klein four (reversals), n=3, p=3", 4, K4, 3, 3, [], [(k1, em_from(4, K4, k1)), (k2, em_from(4, K4, k2))]))
tests.append(("K4 Klein four x labels S3, n=3, p=5", 4, K4, 3, 5, [[1, 0, 2], [1, 2, 0]], [(k1, em_from(4, K4, k1)), (k2, em_from(4, K4, k2))]))
th4 = [(0, 1)] * 4
tests.append(("theta4 edge swap, n=3, p=3", 2, th4, 3, 3, [], [([0, 1], [1, 0, 2, 3])]))
tests.append(("theta4 vertex swap (reversal), n=3, p=3", 2, th4, 3, 3, [], [([1, 0], ["~0", "~1", "~2", "~3"])]))
tests.append(("theta4 edge swap x S3 labels, n=3, p=5", 2, th4, 3, 5, [[1, 0, 2], [1, 2, 0]], [([0, 1], [1, 0, 2, 3])]))
K33 = [(i, 3 + j) for i in range(3) for j in range(3)]
d3 = [1, 2, 0, 4, 5, 3]; sw = [1, 0, 2, 3, 4, 5]; side = [3, 4, 5, 0, 1, 2]
tests.append(("K33 diagonal Z3, n=3, p=2", 6, K33, 3, 2, [], [(d3, em_from(6, K33, d3))]))
tests.append(("K33 diagonal Z3 x labels Z3, n=3, p=2", 6, K33, 3, 2, [[1, 2, 0]], [(d3, em_from(6, K33, d3))]))
tests.append(("K33 swap left 0,1, n=3, p=3", 6, K33, 3, 3, [], [(sw, em_from(6, K33, sw))]))
tests.append(("K33 side swap (reversal), n=3, p=3", 6, K33, 3, 3, [], [(side, em_from(6, K33, side))]))
K5 = [(i, (i + 1) % 5) for i in range(5)] + [(i, (i + 2) % 5) for i in range(5)]
r5 = [1, 2, 3, 4, 0]
tests.append(("K5 rotation Z5, n=3, p=2", 5, K5, 3, 2, [], [(r5, em_from(5, K5, r5))]))
tests.append(("K5 rotation x labels Z3, n=3, p=2", 5, K5, 3, 2, [[1, 2, 0]], [(r5, em_from(5, K5, r5))]))
T = [(0, 1), (0, 2), (0, 2), (1, 2), (1, 2), (1, 2)]
tests.append(("T(1,2,3) parallel swap, n=3, p=3", 3, T, 3, 3, [], [([0, 1, 2], [0, 2, 1, 3, 4, 5])]))
dK4 = K4 * 2
tests.append(("doubled K4 Klein four, n=3, p=3", 4, dK4, 3, 3, [], [(k1, em_from(4, dK4, k1)), (k2, em_from(4, dK4, k2))]))
allok = True
for name, V, ed, n, p, lg, au in tests:
    mine, bases, tors = oconf(V, ed, n, p, lg, au)
    ref, bad, order = run(V, ed, n, p, lg, au, {v: bases[v] for v in range(V) if bases[v] >= 0})
    ok = mine == ref and bad == 0
    allok &= ok
    print(f"{name:42s} |G|={order:3d} oconf19 {mine} independent {ref} chainmap-fail {bad}  {'OK' if ok else 'FAIL'}", flush=True)
print("ALL OK" if allok else "FAILURES")
