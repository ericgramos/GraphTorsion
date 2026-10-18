"""Independent verifier for a p-torsion certificate in H_i(Conf_n(G); Z) (any prime p, default 2).

Shares no code with the C++ engine. It rebuilds the reduced ordered Swiatkowski
complex from its definition (Hainaut-Knudsen-Wawrykow, Def. 3.12):

  * half-edges: edge e = (tail, head) has half-edges 2e (tail end) and 2e+1 (head end);
    the base half-edge h0(v) of a vertex is its lowest-numbered half-edge;
  * a cell: each vertex is empty or carries one particle p together with a non-base
    half-edge h (the generator h - h0(v)); every edge carries an ordered list of particles
    (listed from tail to head); every particle 1..n appears exactly once;
  * degree = number of occupied vertices;
  * boundary: if v_0 < v_1 < ... are the occupied vertices and v_j carries (p, h), then
      d = sum_j (-1)^j ( [p pushed into edge(h) at the end of h]
                         - [p pushed into edge(h0) at the end of h0] ),
    where pushing in at a tail end prepends and at a head end appends.

Certificate (file lines):  PRIME p (optional, default 2)
                           W [coef] <cell>  (a chain w in degree i+1, coefficients 0..p-1)
                           F [coef] <cell>  (a cochain f in degree i, values mod p)
Checks:
  (a) every cell named is a genuine cell of the right degree;
  (b) the integral boundary of w is divisible by p; put z = d(w)/p, an integral i-chain;
  (c) d(z) = 0, so z is an integral i-cycle (and pz = d(w) is a boundary);
  (d) f is a mod-p cocycle: f(d sigma) = 0 mod p for EVERY (i+1)-cell sigma
      (all (i+1)-cells are enumerated independently and counted against a closed formula);
  (e) f(z) is nonzero mod p.
Then z is not a boundary over Z (else f(z) = f(d u) = 0 mod p), but pz is: [z] has order p.
"""
import sys, re, math, itertools, time, gzip
from collections import defaultdict


def load(fn):
    W, F = [], []
    edges = []
    p = 2
    opener = gzip.open if fn.endswith(".gz") else open
    for line in opener(fn, "rt"):
        tag, _, rest = line.rstrip("\n").partition(" ")
        if tag == "GRAPH":
            V, E, n = map(int, rest.split())
        elif tag == "EDGE":
            edges.append(tuple(map(int, rest.split())))
        elif tag == "DEGREE":
            deg = int(rest)
        elif tag == "PRIME":
            p = int(rest)
        elif tag in ("W", "F"):
            first, _, cell = rest.partition(" ")
            coef, cell = (int(first), cell) if first.isdigit() else (1, rest)
            (W if tag == "W" else F).append((coef, cell))
    return V, edges, n, deg, W, F, p


CELL_RE = re.compile(r"V\[(.*)\] E\[(.*)\]$")


def parse_cell(s, V, E):
    m = CELL_RE.match(s)
    vs = m.group(1).split(",")
    es = m.group(2).split("|")
    assert len(vs) == V and len(es) == E, s
    verts = tuple(None if x == "-" else tuple(map(int, x.split("/"))) for x in vs)
    eds = tuple(tuple(map(int, x.split())) if x.strip() else () for x in es)
    return (verts, eds)


def main(fn):
    t0 = time.time()
    V, edges, n, deg, Wraw, Fraw, P = load(fn)
    E = len(edges)
    end_vertex = {}
    for e, (a, b) in enumerate(edges):
        end_vertex[2 * e] = a
        end_vertex[2 * e + 1] = b
    halves = defaultdict(list)
    for h in range(2 * E):
        halves[end_vertex[h]].append(h)
    base = {v: min(hs) for v, hs in halves.items()}

    def valid(cell, k):
        verts, eds = cell
        seen = []
        occ = 0
        for v, x in enumerate(verts):
            if x is not None:
                p, h = x
                if end_vertex.get(h) != v or h == base[v]:
                    return False
                seen.append(p); occ += 1
        for seq in eds:
            seen.extend(seq)
        return occ == k and sorted(seen) == list(range(n))

    def boundary(cell):
        verts, eds = cell
        out = []
        j = 0
        for v, x in enumerate(verts):
            if x is None:
                continue
            p, h = x
            sign = 1 if j % 2 == 0 else -1
            j += 1
            nv = list(verts); nv[v] = None; nv = tuple(nv)
            for hh, c in ((h, sign), (base[v], -sign)):
                e = hh // 2
                ne = list(eds)
                ne[e] = (p,) + eds[e] if hh % 2 == 0 else eds[e] + (p,)
                out.append(((nv, tuple(ne)), c))
        return out

    W = {}
    for a, s_ in Wraw:
        c = parse_cell(s_, V, E); assert c not in W, "duplicate W cell"; W[c] = a
    F = {}
    for a, s_ in Fraw:
        c = parse_cell(s_, V, E); assert c not in F, "duplicate F cell"; F[c] = a % P
    assert all(valid(c, deg + 1) for c in W), "(a) invalid W cell"
    assert all(valid(c, deg) for c in F), "(a) invalid F cell"
    print(f"(a) ok: p = {P}; {len(W)} cells in w (degree {deg+1}), {len(F)} cells in f (degree {deg})", flush=True)

    dw = defaultdict(int)
    for c, a in W.items():
        for f, s in boundary(c):
            dw[f] += a * s
    dw = {k: v for k, v in dw.items() if v}
    assert all(v % P == 0 for v in dw.values()), "(b) boundary of w not divisible by p"
    z = {k: v // P for k, v in dw.items()}
    print(f"(b) ok: d(w) = {P}z with z supported on {len(z)} cells, max |coefficient| {max(abs(v) for v in z.values())}", flush=True)

    dz = defaultdict(int)
    for c, a in z.items():
        for f, s in boundary(c):
            dz[f] += a * s
    assert all(v == 0 for v in dz.values()), "(c) z is not a cycle"
    print("(c) ok: d(z) = 0", flush=True)

    # (d) enumerate every (deg+1)-cell independently and test the cocycle condition
    k = deg + 1
    hv = {v: [h for h in halves[v] if h != base[v]] for v in range(V)}
    count = 0
    bad = 0
    def compositions(m, parts):
        # all ways to write m as an ordered sum of `parts` non-negative integers
        for cuts in itertools.combinations(range(m + parts - 1), parts - 1):
            prev = -1; comp = []
            for c in cuts:
                comp.append(c - prev - 1); prev = c
            comp.append(m + parts - 1 - prev - 1)
            yield comp
    comps = list(compositions(n - k, E))
    for occ in itertools.combinations(range(V), k):
        if any(not hv[v] for v in occ):
            continue
        for parts in itertools.permutations(range(n), k):
            rest = [p for p in range(n) if p not in parts]
            for hs in itertools.product(*[hv[v] for v in occ]):
                verts = [None] * V
                for v, p, h in zip(occ, parts, hs):
                    verts[v] = (p, h)
                verts = tuple(verts)
                for perm in itertools.permutations(rest):
                    for comp in comps:
                        eds = []; i = 0
                        for m in comp:
                            eds.append(tuple(perm[i:i + m])); i += m
                        cell = (verts, tuple(eds))
                        count += 1
                        acc = 0
                        for f, s in boundary(cell):
                            a = F.get(f)
                            if a:
                                acc += s * a
                        if acc % P:
                            bad += 1
    # closed-form count: n! * e_k(d_v - 1) * C(n-k+E-1, E-1)
    dm1 = [len(hv[v]) for v in range(V)]
    ek = sum(math.prod(c) for c in itertools.combinations(dm1, k))
    expected = math.factorial(n) * ek * math.comb(n - k + E - 1, E - 1)
    assert count == expected, f"(d) enumerated {count} cells, formula says {expected}"
    assert bad == 0, f"(d) f fails the cocycle condition on {bad} cells"
    print(f"(d) ok: f is a mod-{P} cocycle (checked on all {count} cells of degree {k})", flush=True)

    pairing = sum(v * F[c] for c, v in z.items() if c in F)
    assert pairing % P != 0, "(e) f(z) = 0 mod p"
    print(f"(e) ok: f(z) = {pairing} = {pairing % P} mod {P}", flush=True)
    print(f"VERIFIED: H_{deg}(Conf_{n}(G); Z) contains an element of order {P}, G = {edges}  [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main(sys.argv[1])
