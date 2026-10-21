"""shortest_loops.py -- find short loops homologous to the second factors c_A of the two tori.

The 1-skeleton of the complex of Conf_B(H minus vertex 1) is a graph (0-cells and 1-cells).  Its simple cycles of
length at most L are enumerated (networkx), each written as a signed 1-chain, and compared with c_A modulo the
boundaries of 2-cells: c_A ~ sum e_i gamma_i holds over Q iff appending c_A - sum e_i gamma_i to the boundary rows does
not raise the rank (mod two large primes), and then over Z because H_1 is torsion-free.  Single loops are tried first,
then pairs (by hashing the normal forms, so the pair search is linear in the number of loops).  If no single loop or pair is found and the environment variable PER_WALK is set, each closed walk of the raw factor
(from walks.py) is shortened separately instead.
Usage: [PER_WALK=1] python3 shortest_loops.py tori_file exact_check.json graph_file L
       exact_check.json: the output of exact_check.py (EXACT_OUT=...), whose "best" vector gives the tori and their
       coefficients; L: the largest loop length to enumerate (12 gives about 18,000 loops, 16 about 425,000).
"""
import itertools
import json
import os
import sys
from collections import defaultdict

import flint
import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from describe_tori import factor_signs, load_graph
from simplify_piece import piece_complex, bd
from walks import closed_walks

P = (1000003, 998244353)


def main(tfn, ecfn, gfn, L):
    """Search for each torus family; prints the loops found as closed walks."""
    V, E, edges = load_graph(gfn)
    co = {k: int(v) for k, v in json.load(open(ecfn))["best"][1].items()}
    cells, hdr, cur = defaultdict(list), {}, None
    for line in open(tfn):
        if line.startswith("Z "):
            cur = line.split()[1]; hdr[cur] = line.strip(); continue
        if cur in co:
            a, c = line.rstrip("\n").split(" ", 1); cells[cur].append((int(a), c))
    groups = {}
    for t, c in co.items():
        A = tuple(sorted(map(int, hdr[t].split("A=")[1].split()[0].split(","))))
        star, piece = factor_signs(cells[t], V, E, 1, set(A))
        if A not in groups:
            groups[A] = [star, defaultdict(int)]
        for k, v in piece.items():
            groups[A][1][k] += c * v
    result = {}
    for A, (star, piece) in sorted(groups.items()):
        B = [p for p in range(5) if p not in A]
        base, pc = piece_complex(edges, B)
        i1 = {c: i for i, c in enumerate(pc[1])}
        n = len(pc[1])
        def key_of(v, p, h, words):
            """The 1-cell 'particle p at vertex v with half-edge h' with the given edge words, in the piece's cell format."""
            verts = [None] * 4; verts[v] = (p, h); return (tuple(verts), tuple(words))
        target = defaultdict(int)
        for (v, p, h, words), x in piece.items():
            if x:
                target[i1[key_of(v, p, h, words)]] += x
        rowsB = []
        for c2 in pc[2]:
            r = defaultdict(int)
            for f, s in bd(c2, base):
                r[i1[f]] += s
            rowsB.append(r)
        # canonical normal form modulo B_1 (mod a large prime): reduce by the reduced row echelon form of the
        # boundary rows.  Two 1-chains are homologous over Q iff their normal forms agree (checked for two primes),
        # and then over Z, because H_1 is torsion-free.
        nforms = {}
        for p in P:
            M = flint.nmod_mat(len(rowsB), n, p)
            for i, r in enumerate(rowsB):
                for j, x in r.items():
                    M[i, j] = x % p
            Er, rk = M.rref()
            piv, jj, R = [], 0, []
            for i in range(rk):
                while int(Er[i, jj]) == 0:
                    jj += 1
                piv.append(jj)
                R.append({j: int(Er[i, j]) for j in range(jj, n) if int(Er[i, j])})
                jj += 1
            nforms[p] = (piv, R)

        def nf(vec, p):
            """Normal form of a 1-chain modulo the boundaries, mod p (reduction by the reduced row echelon form)."""
            piv, R = nforms[p]
            v = defaultdict(int)
            for j, x in vec.items():
                v[j] = x % p
            for pc_, row in zip(piv, R):
                a_ = v.get(pc_, 0)
                if a_:
                    for j, x in row.items():
                        v[j] = (v[j] - a_ * x) % p
            return tuple(sorted((j, x) for j, x in v.items() if x))

        def neg(t, p):
            """Negative of a normal form."""
            return tuple((j, (-x) % p) for j, x in t)

        def add(t1, t2, p):
            """Sum of two normal forms."""
            d = defaultdict(int)
            for j, x in t1:
                d[j] = (d[j] + x) % p
            for j, x in t2:
                d[j] = (d[j] + x) % p
            return tuple(sorted((j, x) for j, x in d.items() if x))

        Gs = nx.Graph()                         # the 1-skeleton
        for c in pc[1]:
            (f1, s1), (f2, s2) = bd(c, base)     # f2: p on the base edge, f1: p on edge(h)
            Gs.add_edge(f2, f1, cell=c)
        cycles = []
        for cyc in nx.simple_cycles(Gs, length_bound=L):
            vec = defaultdict(int)
            for a2, b2 in zip(cyc, cyc[1:] + cyc[:1]):
                c = Gs[a2][b2]["cell"]
                (f1, s1), (f2, s2) = bd(c, base)
                vec[i1[c]] += 1 if (a2 == f2 and b2 == f1) else -1
            cycles.append((len(cyc), dict(vec)))
        cycles.sort(key=lambda t: t[0])
        p0, p1 = P
        tN, tN1 = nf(target, p0), nf(target, p1)
        print(f"A = {A}: {len(cycles)} simple loops of length <= {L}", flush=True)
        forms = [nf(v, p0) for _, v in cycles]
        found = None
        for i, f in enumerate(forms):
            for e in (1, -1):
                if (f if e > 0 else neg(f, p0)) == tN and nf({j: e * x for j, x in cycles[i][1].items()}, p1) == tN1:
                    found = [(e, cycles[i][1])]
                    break
            if found:
                break
        if found:
            print(f"   a single loop of {len(next(iter([found[0][1]])))} one-cells", flush=True)
        else:
            index = defaultdict(list)
            for i, f in enumerate(forms):
                index[f].append((i, 1))
                index[neg(f, p0)].append((i, -1))
            best = None
            for i, f in enumerate(forms):
                for e1 in (1, -1):
                    rest = add(tN, neg(f, p0) if e1 > 0 else f, p0)
                    for k, e2 in index.get(rest, []):
                        tot = cycles[i][0] + cycles[k][0]
                        if best is None or tot < best[0]:
                            s2 = defaultdict(int)
                            for j, x in cycles[i][1].items():
                                s2[j] += e1 * x
                            for j, x in cycles[k][1].items():
                                s2[j] += e2 * x
                            if nf(dict(s2), p1) == tN1:
                                best = (tot, [(e1, cycles[i][1]), (e2, cycles[k][1])])
            found = best[1] if best else None
            print("   no single loop of this length" + (f"; best pair has {best[0]} moves in total" if best else "; no pair either"), flush=True)
        if not found and os.environ.get("PER_WALK"):
            # shorten each closed walk of the raw factor separately (each is a cycle by itself)
            chain0 = {}
            for (v, p, h, words), x in piece.items():
                if x:
                    chain0[(v, p, h, words)] = x
            walks0, _ = closed_walks(chain0, base, edges)
            found = []
            for path in walks0:
                wv = defaultdict(int)
                for cur_, dst, v, p, h, s_ in path:
                    # the 1-cell of this move is "p at v with half-edge h", the other particles as in either end
                    # configuration with p removed; it is traversed with sign s_ (+1: from the base edge into edge(h))
                    ww = tuple(tuple(q for q in w if q != p) for w in cur_)
                    wv[i1[key_of(v, p, h, ww)]] += s_
                wN = nf(dict(wv), p0)
                best1 = None
                for i, f in enumerate(forms):
                    for e in (1, -1):
                        if (f if e > 0 else neg(f, p0)) == wN and nf({j: e * x for j, x in cycles[i][1].items()}, p1) == nf(dict(wv), p1):
                            if best1 is None or cycles[i][0] < best1[0]:
                                best1 = (cycles[i][0], e, cycles[i][1])
                if best1 and best1[0] < len(path):
                    print(f"   walk of {len(path)} moves ~ loop of {best1[0]} moves", flush=True)
                    found.append((best1[1], best1[2]))
                else:
                    print(f"   walk of {len(path)} moves kept", flush=True)
                    found.append((1, dict(wv)))
        result[A] = found
        if found:
            for e, vec in found:
                chain = {}
                for j, x in vec.items():
                    verts, words = pc[1][j]
                    v = next(u for u, y in enumerate(verts) if y is not None)
                    chain[(v, verts[v][0], verts[v][1], words)] = e * x
                w, text = closed_walks(chain, base, edges)
                print(f"   {'+' if e > 0 else '-'} loop:")
                print("\n".join("      " + t for t in text))
        else:
            print("   not found with loops of this length")
    return result


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
