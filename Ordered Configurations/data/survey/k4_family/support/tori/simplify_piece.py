"""simplify_piece.py -- write the second factors of the two tori in terms of basic exchanges.

Each torus of the order-16 class is s_A x c_A: s_A is an exchange of the two particles A at vertex 1 and c_A is a
1-cycle of the other three particles B on H minus vertex 1.  H minus vertex 1 is the edge 23 with two hairs at 2 (the
copies of 12, cut near 1), two hairs at 3 (the copies of 13) and, separately, the leaf edge 01.  Changing c_A by a
boundary of Conf_B(H minus 1) changes the torus by a boundary (s_A x d e = +-d(s_A x e)), so only the class of c_A in
H_1(Conf_B(H minus 1)) matters, and that group is torsion-free.

The complex of Conf_B(H minus 1) is written in the cell language of H with vertex 1 never occupied (vertex 0 is a leaf
and is never occupied either).  A BASIC EXCHANGE is the cycle in which two particles a, b of B pass each other at a
vertex v in {2, 3} (the generator of H_1(Conf_2 of the star of v), six moves), while the third particle q stays at a
fixed spot: a fixed edge and a fixed number of particles between q and the end of that edge away from v.

Method: (1) choose basic exchanges that are independent modulo boundaries (rank mod a large prime); (2) solve
c_A = sum lam_k g_k + boundary for these lam modulo two large primes and lift them by rational reconstruction;
(3) check exactly that the lam are integers and that c_A - sum lam_k g_k is a rational combination of boundaries.
Since H_1 is torsion-free (B_1 is saturated in Z_1), (3) shows the identity holds over Z.

Outcome for the two tori of the order-16 class: H_1(Conf_3(H minus 1)) has rank 40, the basic exchanges span rank
24 of it, and neither second factor lies in their span (the script then prints these ranks instead of coefficients).
So the second factors genuinely involve all three particles; shortest_loops.py describes them by loops instead.
Usage: python3 simplify_piece.py tori_file exact_check.json graph_file > description.txt
"""
import itertools
import json
import math
import os
import sys
from collections import defaultdict
from fractions import Fraction

import flint

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from describe_tori import factor_signs, load_graph, walk

PRIMES = (1000003, 998244353)


def piece_complex(edges, B, frozen=1):
    """All cells (degrees 0, 1, 2) of Conf_B(H) with vertex `frozen` never occupied, and the base half-edges."""
    V = 1 + max(max(e) for e in edges)
    E = len(edges)
    halves = defaultdict(list)
    for e, (a, b) in enumerate(edges):
        halves[a].append(2 * e)
        halves[b].append(2 * e + 1)
    base = {v: min(h) for v, h in halves.items()}
    movable = [v for v in range(V) if v != frozen and len(halves[v]) > 1]
    cells = {0: set(), 1: set(), 2: set()}
    for k in (0, 1, 2):
        for occ in itertools.combinations(movable, k):
            for parts in itertools.permutations(B, k):
                rest = [p for p in B if p not in parts]
                for hs in itertools.product(*[[h for h in halves[v] if h != base[v]] for v in occ]):
                    verts = [None] * V
                    for v, p, h in zip(occ, parts, hs):
                        verts[v] = (p, h)
                    for perm in itertools.permutations(rest):
                        for es in itertools.combinations_with_replacement(range(E), len(rest)):
                            words = [[] for _ in range(E)]
                            for p, e in zip(perm, es):
                                words[e].append(p)
                            cells[k].add((tuple(verts), tuple(tuple(w) for w in words)))
    key = lambda c: (tuple((-1, -1) if x is None else x for x in c[0]), c[1])   # None is not comparable
    return base, {k: sorted(v, key=key) for k, v in cells.items()}


def bd(cell, base):
    """Boundary of a cell of the piece complex, as a list of (face, sign)."""
    verts, words = cell
    out, j = [], 0
    for v, x in enumerate(verts):
        if x is None:
            continue
        p, h = x
        s = 1 if j % 2 == 0 else -1
        j += 1
        nv = list(verts)
        nv[v] = None
        for hh, c in ((h, s), (base[v], -s)):
            w = list(words)
            e = hh // 2
            w[e] = ((p,) + w[e]) if hh % 2 == 0 else (w[e] + (p,))
            out.append(((tuple(nv), tuple(w)), c))
    return out


def basic_exchanges(pc, base, edges, B):
    """(description, {1-cell: coef}) for every basic exchange."""
    out = []
    for v in (2, 3):
        at_v = {e for e, (x, y) in enumerate(edges) if v in (x, y)}
        for a, b in itertools.combinations(B, 2):
            q = [p for p in B if p not in (a, b)][0]
            by_spot = defaultdict(list)
            for c in pc[1]:
                verts, words = c
                occ = [u for u, x in enumerate(verts) if x is not None]
                if occ != [v] or verts[v][0] not in (a, b):
                    continue
                qe = next(e for e, w in enumerate(words) if q in w)
                w = words[qe]
                # the spot of q: its edge and how many particles separate it from the end of that edge AWAY from v
                # (for an edge not at v, from its head end); a and b stay between v and q
                if qe in at_v and edges[qe][1] == v:
                    k = w.index(q)                    # v is the head: the far end is the tail end, index from the left
                else:
                    k = len(w) - 1 - w.index(q)       # count from the head end
                by_spot[(qe, k)].append(c)
            for (qe, qk), cs in sorted(by_spot.items()):
                nodes = list({f for c in cs for f, _ in bd(c, base)})
                ni = {n: i for i, n in enumerate(nodes)}
                M = flint.fmpz_mat(len(nodes), len(cs))
                for j, c in enumerate(cs):
                    for f, s in bd(c, base):
                        M[ni[f], j] += s
                K, nk = M.nullspace()
                if nk == 0:      # one of a, b shares the edge of q and cannot reach v: no exchange here
                    continue
                assert nk == 1, f"expected one exchange cycle, found {nk}"
                vec = [int(K[j, 0]) for j in range(len(cs))]
                g = 0
                for x in vec:
                    g = math.gcd(g, abs(x))
                cyc = {c: x // g for c, x in zip(cs, vec) if x}
                end = "tail" if (qe in at_v and edges[qe][1] == v) else "head"
                out.append((f"exchange of {a} and {b} at vertex {v}, with {q} parked on edge {edges[qe]}#{qe}, {qk} particles from its {end} end", cyc))
    return out


def rank_mod(rows, n, p):
    """Rank mod p of a list of sparse rows."""
    M = flint.nmod_mat(max(1, len(rows)), n, p)
    for i, r in enumerate(rows):
        for j, x in r.items():
            M[i, j] = x % p
    return M.rank()


def ratrec(a, p):
    """Rational reconstruction of a residue mod p."""
    B = math.isqrt(p // 2)
    r0, r1, s0, s1 = p, a % p, 0, 1
    while r1 > B:
        q = r0 // r1
        r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)


def main(tfn, ecfn, gfn):
    """For both tori: the class of the second factor in terms of basic exchanges, if it lies in their span."""
    V, E, edges = load_graph(gfn)
    co = {k: int(v) for k, v in json.load(open(ecfn))["best"][1].items()}
    cells, hdr, cur = defaultdict(list), {}, None
    for line in open(tfn):
        if line.startswith("Z "):
            cur = line.split()[1]
            hdr[cur] = line.strip()
            continue
        if cur in co:
            a, c = line.rstrip("\n").split(" ", 1)
            cells[cur].append((int(a), c))
    groups = {}
    for t, c in co.items():
        A = tuple(sorted(map(int, hdr[t].split("A=")[1].split()[0].split(","))))
        star, piece = factor_signs(cells[t], V, E, 1, set(A))
        if A not in groups:
            groups[A] = [star, defaultdict(int)]
        assert groups[A][0] == star
        for k, v in piece.items():
            groups[A][1][k] += c * v
    halves = defaultdict(list)
    for e, (a, b) in enumerate(edges):
        halves[a].append(2 * e)
        halves[b].append(2 * e + 1)
    baseH = {v: min(h) for v, h in halves.items()}
    print("The cycle y* of Conf_5(H) is the sum of two tori s_A x c_A (H = edges", edges, ").")
    for A, (star, piece) in sorted(groups.items()):
        B = [p for p in range(5) if p not in A]
        base, pc = piece_complex(edges, B)
        i1 = {c: i for i, c in enumerate(pc[1])}
        n = len(pc[1])
        target = {}
        for (v, p, h, words), x in piece.items():
            if x:
                verts = [None] * 4
                verts[v] = (p, h)
                target[i1[(tuple(verts), tuple(words))]] = x
        bdrows = []
        for c2 in pc[2]:
            r = defaultdict(int)
            for f, s in bd(c2, base):
                r[i1[f]] += s
            bdrows.append(dict(r))
        gens = basic_exchanges(pc, base, edges, B)
        gvec = [{i1[c]: x for c, x in g.items()} for _, g in gens]
        p0 = PRIMES[0]
        rB = rank_mod(bdrows, n, p0)
        # independent exchanges modulo boundaries, greedily
        chosen, rows, r = [], list(bdrows), rB
        for k, g in enumerate(gvec):
            if rank_mod(rows + [g], n, p0) > r:
                rows.append(g); chosen.append(k); r += 1
        i0 = {c: i for i, c in enumerate(pc[0])}
        d1 = []
        for c in pc[1]:
            r1 = defaultdict(int)
            for f, s_ in bd(c, base):
                r1[i0[f]] += s_
            d1.append(dict(r1))
        h1 = n - rank_mod(d1, len(pc[0]), p0) - rB      # rank of H_1(Conf_B(H minus 1))
        # solve target = sum lam_k g_k + boundaries, mod two primes (unique in lam since chosen are independent)
        lams = []
        for p in PRIMES:
            m = len(bdrows) + len(chosen)
            M = flint.nmod_mat(n, m + 1, p)
            for i, row in enumerate(bdrows):
                for j, x in row.items():
                    M[j, i] = x % p
            for k, ci in enumerate(chosen):
                for j, x in gvec[ci].items():
                    M[j, len(bdrows) + k] = x % p
            for j, x in target.items():
                M[j, m] = x % p
            Er, rk = M.rref()
            # read off lam: pivot columns among the generator columns
            lam = [0] * len(chosen)
            piv, jj = [], 0
            for i in range(rk):
                while int(Er[i, jj]) == 0:
                    jj += 1
                piv.append(jj)
                jj += 1
            if m in piv:
                print(f"  [B = {B}: H_1 rank {h1}, exchanges span rank {r - rB}; the factor is not in their span]")
                break
            for i, pj in enumerate(piv):
                if pj >= len(bdrows):
                    lam[pj - len(bdrows)] = int(Er[i, m])
            lams.append([ratrec(x, p) for x in lam])
        if len(lams) < 2:
            continue
        assert lams[0] == lams[1], "the two primes disagree"
        lam = lams[0]
        assert all(x.denominator == 1 for x in lam), "non-integral coefficients"
        lam = [int(x) for x in lam]
        # exact check: target - sum lam g is in the rational span of the boundaries (rank does not go up)
        diff = defaultdict(int, target)
        for l, ci in zip(lam, chosen):
            for j, x in gvec[ci].items():
                diff[j] -= l * x
        diff = {j: x for j, x in diff.items() if x}
        assert rank_mod(bdrows + [diff], n, PRIMES[1]) == rank_mod(bdrows, n, PRIMES[1])
        print(f"\n=== torus with the star at vertex 1 for particles {list(A)} ===")
        for l in walk(star, baseH, edges, f"s_{''.join(map(str, A))}: exchange of {A[0]} and {A[1]} at vertex 1"):
            print("  " + l)
        print(f"  c_{''.join(map(str, A))} (particles {B} on H minus vertex 1), up to boundaries:")
        print(f"    H_1 has rank {h1}; the {len(gens)} basic exchanges span rank {r - rB} of it; the factor is")
        for l, ci in zip(lam, chosen):
            if l:
                print(f"    {l:+d} x {gens[ci][0]}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
