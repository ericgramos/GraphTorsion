"""gconf.py -- Python driver for the engine ``engine/oconf``.

This module is the convenient entry point for computing the homology of ordered configuration spaces of
graphs from Python.  It writes the graph to the engine's standard input, parses the JSON the engine prints,
finishes the computation where the engine leaves work (exact Smith normal form of the small non-unit
residual blocks, via FLINT), and returns plain Python dictionaries.

Requirements:  Python 3.8+, ``python-flint`` (``pip install python-flint``), and a compiled engine.
The engine binary is looked up in this order:
    1. the environment variable OCONF_BIN (a path), if set;
    2. ../engine/oconf relative to this file.
Temporary residual files go to the directory GCONF_SCRATCH (default: the system temporary directory).

Main functions
--------------
homology(edges, n, prime=0, colors=None)
    Betti numbers and torsion of H_*(Conf_n(G)) over Z (prime = 0) or dimensions over F_prime.
local_torsion(edges, n, p, k=None, colors=None)
    Exact p-primary torsion using the local Smith form modulo p^k (engine option -q p k).
gal_euler(V, edges, n)
    Euler characteristic of Conf_n(G) from Gal's formula (an independent check of the cell counts).
abrams_homology(edges, n)
    Homology of Abrams' discretized configuration space (an independent model, small cases only).

Graph format: a list of pairs (u, v); vertices may be any hashable labels; loops (u, u) and repeated
pairs (multi-edges) are allowed.  Edge (u, v) is oriented u -> v; the orientation does not change the
homology but does fix the cell conventions of the engine.

Example:
    >>> from gconf import homology
    >>> T133 = [(0, 1)] + [(0, 2)] * 3 + [(1, 2)] * 3
    >>> r = homology(T133, 6)          # about 8 minutes and 4 GB
    >>> r["betti"], r["torsion"]
    ([1, 195, 224725, 43091], [[], [], [3, 3, 3, 3, 3, 3], []])
"""
import json, os, subprocess, tempfile, math
from collections import Counter

import flint  # python-flint: exact Smith normal form (fmpz_mat.snf) and ranks mod p (nmod_mat)

HERE = os.path.dirname(os.path.abspath(__file__))
ENGINE = os.environ.get("OCONF_BIN", os.path.join(HERE, "..", "engine", "oconf"))
SCRATCH = os.environ.get("GCONF_SCRATCH", tempfile.gettempdir())


def normalize(edges):
    """Relabel the vertices as 0..V-1 (in sorted order of the original labels).

    Returns (V, list of edges as integer pairs).  The engine requires vertices 0..V-1."""
    verts = sorted({x for e in edges for x in e})
    idx = {v: i for i, v in enumerate(verts)}
    return len(verts), [(idx[u], idx[v]) for u, v in edges]


def smooth(edges):
    """Remove vertices of degree 2 (this does not change the topological space, hence not Conf_n).

    A degree-2 vertex w with two distinct incident edges a-w, w-b is replaced by one edge a-b.  A degree-2
    vertex carrying a single loop is a circle component and is kept as it is."""
    edges = [tuple(e) for e in edges]
    changed = True
    while changed:
        changed = False
        deg = Counter()
        for u, v in edges:
            deg[u] += 1
            deg[v] += 1
        for w, d in deg.items():
            if d != 2:
                continue
            inc = [i for i, (u, v) in enumerate(edges) if u == w or v == w]
            if len(inc) == 1:          # a loop at w: w is a circle component
                continue
            i, j = inc
            a = edges[i][0] if edges[i][1] == w else edges[i][1]
            b = edges[j][0] if edges[j][1] == w else edges[j][1]
            edges = [e for k, e in enumerate(edges) if k not in (i, j)] + [(a, b)]
            changed = True
            break
    return edges


def snf_divisors(fn):
    """Read a residual matrix file written by the engine ("rows cols" then "r c v" lines) and return
    (rank, list of nonzero elementary divisors) computed exactly over Z with FLINT."""
    with open(fn) as f:
        r, c = map(int, f.readline().split())
        if r == 0 or c == 0:
            return 0, []
        M = [[0] * c for _ in range(r)]
        for line in f:
            i, j, v = map(int, line.split())
            M[i][j] = v
    S = flint.fmpz_mat(M).snf()
    d = [abs(int(S[i, i])) for i in range(min(r, c)) if S[i, i] != 0]
    return len(d), d


def snf_divisors_modp(fn, p):
    """Rank modulo p of a residual matrix file (used when the engine ran with prime coefficients)."""
    with open(fn) as f:
        r, c = map(int, f.readline().split())
        if r == 0 or c == 0:
            return 0, []
        M = [[0] * c for _ in range(r)]
        for line in f:
            i, j, v = map(int, line.split())
            M[i][j] = v
    return flint.nmod_mat(M, p).rank(), []


def gal_euler(V, edges, n):
    """Euler characteristic of Conf_n(G) by Gal's formula (Gal 2001):

        sum_n chi(Conf_n G) t^n / n!  =  (1 - t)^(-E)  *  prod_v (1 + (1 - d_v) t),

    d_v the valence of v.  Returns chi(Conf_n G)."""
    deg = Counter()
    for u, v in edges:
        deg[u] += 1
        deg[v] += 1
    E = len(edges)
    poly = [1]                                   # coefficients of prod_v (1 + (1 - d_v) t)
    for v in range(V):
        a = 1 - deg[v]
        poly = [(poly[i] if i < len(poly) else 0) + (a * poly[i - 1] if i >= 1 else 0) for i in range(len(poly) + 1)]
    coef = 0                                     # coefficient of t^n in the product with (1-t)^(-E)
    for i, p in enumerate(poly):
        j = n - i
        if j >= 0:
            coef += p * math.comb(E + j - 1, j) if E > 0 else (p if j == 0 else 0)
    return coef * math.factorial(n)


def homology(edges, n, prime=0, tag=None, keep=False, timeout=None, colors=None):
    """Homology of Conf_n(G) (or of the quotient by a Young subgroup, with colors).

    Parameters
        edges   list of (u, v) pairs
        n       number of particles (<= 10)
        prime   0 for integer coefficients (exact; torsion is reported), otherwise a prime p: the result
                is then dim H_i(Conf_n; F_p) in "betti" (no torsion list)
        colors  optional word like "aabc": particles with equal letters become indistinguishable
                (quotient by the Young subgroup; "aaa..." gives the unordered configuration space)
        keep    keep the residual files instead of deleting them
    Returns a dict with
        cells    [|C_0|, |C_1|, ...]           betti   [b_0, b_1, ...]
        torsion  list per degree i of the elementary divisors > 1 of d_{i+1}, i.e. the torsion
                 subgroup of H_i is the product of Z/d over that list
        rank     ranks of the boundary maps    chi, gal  Euler characteristic of the complex and from
                                               Gal's formula (they must agree)
        residual (rows, cols, max |entry|) of each residual block
    """
    V, edges = normalize(edges)
    E = len(edges)
    tag = tag or f"g{os.getpid()}"
    prefix = os.path.join(SCRATCH, tag)
    inp = f"{V} {E} {n}\n" + "".join(f"{u} {v}\n" for u, v in edges)
    args = [ENGINE, "-o", prefix] + (["-p", str(prime)] if prime else []) + (["-k", colors] if colors else [])
    out = subprocess.run(args, input=inp, capture_output=True, text=True, timeout=timeout)
    if out.returncode != 0:
        raise RuntimeError(out.stderr)
    info = json.loads(out.stdout)
    cells = info["cells"]
    D = len(cells) - 1
    rank = [0] * (D + 2)
    divisors = [[] for _ in range(D + 2)]       # nonunit elementary divisors of d_k
    resid = {}
    for m in info["maps"]:
        k = m["k"]
        if prime:
            rr, _ = snf_divisors_modp(m["res_file"], prime)
            rank[k] = m["pivots"] + rr
        else:
            rr, d = snf_divisors(m["res_file"])
            rank[k] = m["pivots"] + rr          # unit pivots + rank of the residual block
            divisors[k] = [x for x in d if x != 1]
        resid[k] = (m["res_rows"], m["res_cols"], m["maxabs"])
        if not keep:
            os.remove(m["res_file"])
    betti = [cells[k] - rank[k] - rank[k + 1] for k in range(D + 1)]
    torsion = [divisors[k + 1] for k in range(D + 1)]
    chi = sum((-1) ** k * c for k, c in enumerate(cells))
    young = math.prod(math.factorial(c) for c in Counter(colors).values()) if colors else 1
    return dict(V=V, E=E, n=n, cells=cells, betti=betti, torsion=torsion, rank=rank[:D + 2],
                chi=chi, gal=gal_euler(V, edges, n) // young, residual=resid, prime=prime)


def sparse_snf(nr, nc, triplets, tag="m"):
    """Rank and nonunit elementary divisors of an arbitrary sparse integer matrix, given as (row, col,
    value) triplets.  Uses the engine's generic mode (-m) for unit-pivot elimination and FLINT for the
    residual block."""
    prefix = os.path.join(SCRATCH, f"{tag}{os.getpid()}")
    inp = f"{nr} {nc}\n" + "".join(f"{r} {c} {v}\n" for r, c, v in triplets)
    out = subprocess.run([ENGINE, "-m", "-o", prefix], input=inp, capture_output=True, text=True)
    info = json.loads(out.stdout)
    rr, d = snf_divisors(info["res_file"])
    os.remove(info["res_file"])
    return info["pivots"] + rr, [x for x in d if x != 1]


# ---------------- independent check: Abrams' discretized model ----------------

def subdivide(edges, m):
    """Subdivide every edge into m edges (Abrams' model needs a sufficiently subdivided graph)."""
    V, edges = normalize(edges)
    nxt = V
    out = []
    for u, v in edges:
        prev = u
        for _ in range(m - 1):
            out.append((prev, nxt))
            prev = nxt
            nxt += 1
        out.append((prev, v))
    return out


def abrams_homology(edges, n):
    """Integral homology of Abrams' discretized configuration space D_n(G) (Abrams 2000).

    D_n(G) is the cube complex of ordered n-tuples of cells (vertices and closed edges) of G with pairwise
    disjoint closures; for a sufficiently subdivided simple graph it is homotopy equivalent to Conf_n(G).
    The boundary of an edge factor (a, b) is b - a.  Dense Smith normal form via FLINT (through
    sparse_snf), so only small cases are feasible.  Used to validate the Swiatkowski model."""
    from itertools import permutations, combinations
    V, edges = normalize(edges)
    cells_all = [("v", v) for v in range(V)] + [("e", i) for i in range(len(edges))]

    def closure(c):
        return {c[1]} if c[0] == "v" else set(edges[c[1]])

    by_dim = {}

    def rec(chosen, used):                     # choose n cells with pairwise disjoint closures, in order
        if len(chosen) == n:
            d = sum(1 for c in chosen if c[0] == "e")
            by_dim.setdefault(d, []).append(tuple(chosen))
            return
        for c in cells_all:
            cl = closure(c)
            if cl & used:
                continue
            rec(chosen + [c], used | cl)

    rec([], set())
    D = max(by_dim)
    index = {d: {c: i for i, c in enumerate(by_dim[d])} for d in by_dim}
    ranks, divs = {}, {}
    for d in range(1, D + 1):
        rows = len(by_dim[d - 1])
        cols = len(by_dim[d])
        trip = []
        for j, cell in enumerate(by_dim[d]):
            s = 0                              # number of edge factors before the current one (sign)
            for pos, c in enumerate(cell):
                if c[0] != "e":
                    continue
                a, b = edges[c[1]]
                sign = (-1) ** s
                s += 1
                for vtx, cf in ((b, sign), (a, -sign)):
                    f = list(cell)
                    f[pos] = ("v", vtx)
                    trip.append((j, index[d - 1][tuple(f)], cf))
        ranks[d], divs[d] = sparse_snf(cols, rows, trip)
    cells = [len(by_dim.get(d, [])) for d in range(D + 1)]
    betti = [cells[k] - ranks.get(k, 0) - ranks.get(k + 1, 0) for k in range(D + 1)]
    torsion = [divs.get(k + 1, []) for k in range(D + 1)]
    return dict(cells=cells, betti=betti, torsion=torsion)


def local_torsion(edges, n, p, k=None, tag=None, colors=None, extra=()):
    """Exact p-primary torsion via the local Smith form over Z_(p), computed modulo p^k (engine -q p k).

    k defaults to the largest exponent with p^k < 32 (16 for p = 2, 27 for p = 3, ...).  Returns a dict
    with "torsion": per degree i the list of p-power orders of the cyclic p-primary summands of H_i that
    have order < p^k, "rank_local" (per map, the number of elementary divisors of order < p^k, i.e. the sum
    of the stages) and "maps" (raw engine output: "stages" per map).
    IMPORTANT: summands of order >= p^k are invisible mod p^k.  They are absent exactly when rank_local
    equals the rank over a large prime (homology(edges, n, prime=1000003)["rank"]); check this.  extra: further engine options,
    e.g. ("-g", "1,2,3,4,0", "-F") for a quotient by a label group of order prime to p."""
    if k is None:
        k = 1
        while p ** (k + 1) < 32:
            k += 1
    V, edges = normalize(edges)
    tag = tag or f"loc{os.getpid()}"
    prefix = os.path.join(SCRATCH, tag)
    inp = f"{V} {len(edges)} {n}\n" + "".join(f"{u} {v}\n" for u, v in edges)
    args = [ENGINE, "-o", prefix, "-q", str(p), str(k)] + (["-k", colors] if colors else []) + list(extra)
    out = subprocess.run(args, input=inp, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(out.stderr)
    info = json.loads(out.stdout)
    D = len(info["cells"]) - 1
    tors = [[] for _ in range(D + 1)]
    rank = [0] * (D + 2)
    for m in info["maps"]:
        if os.path.exists(m["res_file"]):
            os.remove(m["res_file"])
        st = m["stages"]                       # st[s] = number of elementary divisors equal to p^s
        rank[m["k"]] = sum(st)
        tors[m["k"] - 1] = [p ** s for s, c in enumerate(st) for _ in range(c) if s > 0]
    return dict(cells=info["cells"], torsion=tors, rank_local=rank, maps=info["maps"])


def summarize(json_text_or_dict):
    """Summarise one engine JSON output (string or dict): dimensions of H_i with the engine's coefficients
    (computed from the stage-0 ranks for local runs, from the pivots otherwise) and any p-torsion found.

    Returns (dims, torsion) where torsion maps i -> list of stage counts [stage1, stage2, ...] of H_i."""
    d = json.loads(json_text_or_dict) if isinstance(json_text_or_dict, str) else json_text_or_dict
    c = d["cells"]
    r = [0] * (len(c) + 1)
    tors = {}
    for m in d["maps"]:
        st = m.get("stages")
        r[m["k"]] = st[0] if st else m["pivots"]
        if st and (len(st) > 1 or m.get("leftover_nnz")):
            tors[m["k"] - 1] = st[1:]
    dims = [c[i] - r[i] - r[i + 1] for i in range(len(c))]
    return dims, tors
