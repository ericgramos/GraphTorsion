#!/usr/bin/env python3
"""
A large search for exotic torsion in H_2: odd torsion in planar graphs, and torsion of order >= 4.

THE QUESTIONS
-------------
  (1) Does any PLANAR graph have torsion other than 2-torsion in H_2(B_k)?  (Every planar
      example so far has exactly Z/2 summands.)
  (2) Which graphs have "strange" torsion in H_2: a primary summand Z/p^s with p^s >= 4 (Z/4,
      Z/5, Z/8, Z/9, ...)?  The first such example, (Z/4)^2 in H_2(B_3) of K_{4,5} plus a
      triangle on the large side, was found by supergraph_walks.py.

THE METHOD
----------
For each graph Gamma and particle number k, the torsion of H_2 is the torsion of coker d_3
("torsion is a cokernel"), so only the degree-2 and degree-3 slices of the fully reduced
complex are built (gbt.lean).  Then:

  1. Unit pivots are peeled with gbt.linalg.peel_sweeps.  This is a unimodular change of basis,
     so the invariant factors other than 1 are untouched; what is left (the residual) is
     usually tens to hundreds of columns.
  2. r = rank of the residual over F_1000003 (standing in for Q), by dense elimination with
     numpy (`dense_rank`): the residual has few columns but can be dense, and there a dense
     elimination is many times faster than the sparse Markowitz routine.
  3. For every prime p <= 97: the number of invariant factors divisible by p is r minus the
     rank over F_p.  If it is positive, the local Smith form at p (gbt.linalg.local_smith,
     computed modulo p^K) gives the exact p-primary part: stage s counts the invariant factors
     of valuation exactly s.  K starts at 6 for p = 2 and 3 otherwise, and is doubled until the
     stage counts add up to r -- that is, until nothing of order >= p^K is hidden.

The exact Smith form over Z is NOT used: on these residuals it can stall on entry growth (the
K_{4,5}-plus-triangle residual has 13 columns and still defeated it), whereas the local Smith
form never takes a gcd.  The price is that primes above 97 are not examined; a summand Z/p with
p > 97 would be missed.  Each record says which primes were checked.

THE FAMILIES (each graph is simple; multigraphs enter through subdivision, which does not change
the configuration space)

    P8       all 385 connected planar graphs on 8 vertices with minimum degree 3 (nauty), k = 4, 5
    P9       all 3,898 connected planar graphs on 9 vertices with minimum degree 3, k = 4
    P10      a random sample of the 48,754 planar graphs on 10 vertices, minimum degree 3, k = 4
    K45AB    K_{4,5} plus any graph on each side, up to isomorphism, k = 3
    K46B     K_{4,6} plus any graph on the 6-vertex side, up to isomorphism, k = 3
    K55B     K_{5,5} plus any graph on one side, k = 3
    R9       a random sample of the 84,242 connected graphs on 9 vertices, minimum degree 3, k = 3
    MULTI    multigraphs made simple by subdividing repeated edges: K_4, K_5 and K_{3,3} with
             some edges doubled or tripled, the "theta" multigraphs, k = 3, 4, 5

  second phase (run with --family K45AB4 K55B4 K56B K57B K47B MP R10D):

    K45AB4   the K45AB graphs whose H_2 at k = 3 has torsion other than Z/2, again at k = 4
             (needs the K45AB results in the data file; it is empty otherwise)
    K55B4    the K55B graphs again, at k = 4 (only 14 of the 34 fit under --cap 900000)
    K56B     K_{5,6} plus any graph on the 6-vertex side, k = 3
    K57B     K_{5,7} plus a random sample of graphs on the 7-vertex side, k = 3
    K47B     K_{4,7} plus the same sample of 7-vertex graphs, k = 3
    MP       complete multipartite graphs K_{3,3,3}, K_{2,2,2,2}, ..., k = 3, 4
    R10D     random dense graphs on 10 vertices (24 to 34 edges, minimum degree 3), k = 3

USAGE
-----
    python3 experiments/torsion_hunt.py [--family F ...] [--sample 600] [--cap 400000]
                                        [--timeout 600] [--workers 2]

    --family    which families to run, in the order given (default: the first-phase families,
                cheap first: K45AB K46B MULTI K55B R9 P8 P9 P10; the second phase must be named)
    --sample    size of the random samples R9 and P10 (default 600; fixed seed, reproducible)
    --cap       skip a job whose degree-2 and degree-3 slices exceed this many cells (skipped
                jobs are not recorded; rerun with a larger cap to include them)
    --timeout   give up on a single job after this many seconds (it is recorded as such)
    --workers   parallel processes (default 2)

Results are appended to data/torsion_hunt.jsonl, one JSON object per line:
{"family", "name", "k", "N", "M", "planar", "edges", "rank", "primes_checked",
 "torsion": {"p": [stage_0, stage_1, ...]}, "summary", "seconds"} where torsion[p][s] is the
number of invariant factors of d_3 of p-adic valuation exactly s, and "summary" is the torsion of
H_2 as a primary decomposition (e.g. "(Z/2)^11 + Z/3 + (Z/4)^2").  Jobs already in the file are
skipped, so the search resumes.  Lines marked "!!" in the output are hits: odd torsion in a
planar graph, or a summand of order >= 4.  The run ends with a list of all hits.  A job that
exceeds --timeout is recorded with "summary": null, "timed_out" and "cells" (none did).

WHAT THE OUTPUT LOOKS LIKE
--------------------------
    6332 jobs to run
    !! K45AB    K45+A0+B5        k=3 N= 9 M=23        H_2: (Z/2)^11 + Z/3 + (Z/4)^2  [0.2s]
       P9       P9 HCQuTnl       k=4 N= 9 M=19 planar H_2: Z/2                      [6.5s]
    ...
    7291 jobs finished.  HITS (odd torsion in a planar graph, or a summand of order >= 4):
      K45AB    K45+A0+B5             k=3 N=9 M=23 planar=False  H_2: (Z/2)^11 + Z/3 + ...

(Columns narrowed here.)  Columns: hit marker, family, job name (for nauty families the
graph6 string of the graph, so the graph can be rebuilt with networkx.from_graph6_bytes),
k, vertices, edges, planarity, the torsion of H_2 as a primary decomposition, and the time.  K45+Ai+Bj means K_{4,5} plus the i-th
graph on the 4-vertex side and the j-th graph on the 5-vertex side, numbered as `side_graphs`
lists them (A0 = no edges, A1 = one edge, A3 = two disjoint edges; B5 = a triangle, B9 = a
triangle with a pendant edge, B15 = K_4 minus an edge, B22 = K_4; see `side_graphs`).

COST
----
Building the families takes about a minute and a half (nauty enumerates 84,242 graphs on 9
vertices and the planar graphs on 10).  A job at k = 3 takes 0.1 to 30 seconds; a planar graph
on 9 vertices at k = 4 a few seconds; K_{4,5} supergraphs at k = 4 one to three minutes; K_{5,5}
supergraphs at k = 4 three minutes and about 2 GB each.  The stored data (7,291 jobs) took about
twelve hours on two cores.  Exit code: always 0; the findings are in the output.

PITFALLS
--------
* Only H_2 is examined, and only primes p <= 97.
* The planar families are exhaustive only for minimum degree 3: graphs with leaves are left out,
  and leaves can matter (see minimal_carriers.py).
* nauty is called as `nauty-geng` and `nauty-planarg` (the Debian/Ubuntu names).  On macOS
  (Homebrew) the programs are `geng` and `planarg`: change the two names in `geng` and
  `planar_g6`.
* A job that crashes (for example by running out of memory) leaves no record and is retried on
  the next run.  Never stop the search with `pkill -f torsion_hunt`: that pattern also matches
  the shell that started it.  Use the process number.
"""

import argparse
import itertools
import json
import multiprocessing as mp
import os
import random
import re
import subprocess
import sys
import time
import warnings

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import networkx as nx                                                              # noqa: E402

from gbt import Graph                                                              # noqa: E402
from gbt.lean import boundary, slice_cells                                         # noqa: E402
from gbt.linalg import local_smith, peel_sweeps                                    # noqa: E402
from gbt.sizes import reduced_cell_count_by_degree                                 # noqa: E402
from gbt.swiatkowski import edge_monomials                                         # noqa: E402

OUT = os.path.join(ROOT, "data", "torsion_hunt.jsonl")
#: a large prime standing in for Q: the rank over F_BIG equals the rank over Q unless BIG divides
#: some invariant factor, which would itself be a (very) strange torsion summand
BIG = 1000003
#: the primes whose torsion is examined: 2, 3, 5, ..., 97
PRIMES = [p for p in range(2, 98) if all(p % q for q in range(2, int(p ** 0.5) + 1))]


# --------------------------------------------------------------------------------- the engine

def h2_torsion(G, k, d=2):
    """The p-primary parts of tors H_d(B_k G) for every prime p <= 97.

    G is a gbt.Graph with no isolated vertex, k the number of particles, d the degree (2 by
    default; any d >= 0 works).  Returns (tors, r): tors is {p: stages} for the primes that
    occur, where stages[s] is the number of invariant factors of d_{d+1} of p-adic valuation
    exactly s (stages[0] counts the units); r is the rank of the residual over F_BIG.  The
    p-primary part of tors H_d is then the sum over s >= 1 of (Z/p^s)^{stages[s]}.  An empty
    dict means no torsion at any prime <= 97.  The method is described in the module docstring.
    """
    mon = edge_monomials(G.M, k)
    # Index the degree-d cells (the rows of d_{d+1}); only this slice and the degree-(d+1)
    # slice are ever built, which is what makes large k and dense graphs affordable.
    low = {key: i for i, key in enumerate(slice_cells(G, k, d, mon))}
    cols = boundary(G, k, d + 1, mon, low)
    del low
    # Peel unit pivots: a unimodular change of basis, so every invariant factor > 1 survives in
    # the residual and the peeled pivots only account for invariant factors equal to 1.
    res, _ = peel_sweeps(cols)
    del cols
    if not res:
        return {}, 0
    r = dense_rank(res, BIG)
    out = {}
    for p in PRIMES:
        if dense_rank(res, p) == r:
            continue                       # no invariant factor divisible by p
        # Work modulo p^K and make K large enough that nothing is hidden: summands of order
        # >= p^K would be zero modulo p^K, and then the stage counts would fall short of r.
        K = 6 if p == 2 else 3
        while True:
            st = local_smith(res, p, K)
            if sum(st) == r:               # nothing of order >= p^K is hidden
                break
            K *= 2
        out[p] = st
    return out, r


def dense_rank(cols, p):
    """Rank over F_p of a (small) residual matrix, by dense elimination with numpy.

    The residual left by the unit-pivot peel has few columns but can be fairly dense, and the
    sparse Markowitz routine spends most of its time on heap operations there; a dense
    elimination on the rows that actually occur is much faster.  Entries are reduced mod p in
    Python first (they can exceed 64 bits), and p < 2^31 keeps every product inside int64.
    """
    import numpy as np
    rows = sorted({r for c in cols for r in c})
    index = {r: i for i, r in enumerate(rows)}
    A = np.zeros((len(cols), len(rows)), dtype=np.int64)     # one row per column of the matrix
    for j, c in enumerate(cols):
        for r, v in c.items():
            A[j, index[r]] = v % p
    # Row reduction over F_p: for each column of A in turn, find a row (below the rows already
    # used as pivots) with a nonzero entry, scale it to 1, and clear that column in every row
    # below it.  The number of pivots found is the rank.
    rank = 0
    nr, nc = A.shape
    for col in range(nc):
        if rank == nr:
            break
        piv = np.nonzero(A[rank:, col])[0]
        if piv.size == 0:
            continue
        i = rank + piv[0]
        if i != rank:
            A[[rank, i]] = A[[i, rank]]
        inv = pow(int(A[rank, col]), -1, p)
        A[rank] = (A[rank] * inv) % p
        below = np.nonzero(A[rank + 1:, col])[0] + rank + 1
        if below.size:
            A[below] = (A[below] - np.outer(A[below, col], A[rank])) % p
        rank += 1
    return rank


def summary(tors):
    """'(Z/2)^11 + Z/3 + (Z/4)^2' from {p: stages}, ordered by the order of the summand."""
    parts = []
    for p, st in tors.items():
        for s, c in enumerate(st):
            if s and c:
                parts.append((p ** s, c))
    parts.sort()
    return " + ".join(("Z/%d" % q) if c == 1 else "(Z/%d)^%d" % (q, c) for q, c in parts) or "0"


def strange(tors, planar):
    """True for a hit: odd torsion in a planar graph, or any summand of order >= 4.

    `tors` is {p: stages} as returned by h2_torsion (keys may be int); `planar` a bool.  A
    summand Z/p^s with s >= 1 counts when p^s >= 4 (Z/4, Z/5, Z/7, Z/8, Z/9, ...), or when the
    graph is planar and p is odd (so Z/3 in a planar graph is a hit too).
    """
    for p, st in tors.items():
        for s, c in enumerate(st):
            if s and c and (p ** s >= 4 or (planar and p != 2)):
                return True
    return False


# ------------------------------------------------------------------------------- the families

def geng(n, extra=()):
    """graph6 strings of all connected graphs on n vertices with minimum degree 3, one per
    isomorphism class, as printed by nauty-geng -c -d3 n (extra: further geng options)."""
    out = subprocess.run(["nauty-geng", "-c", "-d3", "-q", *extra, str(n)], capture_output=True,
                         text=True, check=True).stdout.split()
    return out


def planar_g6(n):
    """graph6 strings of the connected planar graphs on n vertices with minimum degree 3."""
    p1 = subprocess.Popen(["nauty-geng", "-c", "-d3", "-q", str(n)], stdout=subprocess.PIPE)
    out = subprocess.run(["nauty-planarg", "-q"], stdin=p1.stdout, capture_output=True,
                         text=True, check=True).stdout.split()
    p1.wait()
    return out


def from_g6(s):
    """networkx graph from a graph6 string."""
    return nx.from_graph6_bytes(s.encode())


def side_graphs(n):
    """One representative of every graph on n vertices (labelled 0..n-1), up to isomorphism.

    The order is fixed and reproducible: by number of edges, then by the order in which
    itertools.combinations meets edge sets, keeping the first representative of each class.
    Job names such as K45+A3+B5 refer to positions in these lists (11 graphs on 4 vertices, 34
    on 5, 156 on 6).  Brute force, so only for n <= 6.
    """
    reps = []
    pairs = list(itertools.combinations(range(n), 2))
    for r in range(len(pairs) + 1):
        for S in itertools.combinations(pairs, r):
            g = nx.Graph()
            g.add_nodes_from(range(n))
            g.add_edges_from(S)
            if not any(nx.is_isomorphic(g, h) for h in reps):
                reps.append(g)
    return reps


def bipartite_plus(a, b, sideA, sideB):
    """Edges of K_{a,b} (A = 0..a-1, B = a..a+b-1) plus graphs on the two sides."""
    E = [(i, a + j) for i in range(a) for j in range(b)]
    E += [tuple(e) for e in sideA.edges()]
    E += [(a + u, a + v) for u, v in sideB.edges()]
    return E


def subdivide_multi(V, medges):
    """A simple graph homeomorphic to the multigraph on 0..V-1 with the given edge list
    (repeats allowed): every repeat of an edge after the first is subdivided once."""
    seen = set()
    out = []
    nxt = V
    for u, v in medges:
        key = (min(u, v), max(u, v))
        if key not in seen:
            seen.add(key)
            out.append(key)
        else:
            out += [(u, nxt), (nxt, v)]
            nxt += 1
    return out


def families(sample):
    """{name: list of (job name, edge list, list of k)} in the cheap-first default order."""
    F = {}

    # K_{4,5} plus a graph on each side.  Different (A, B) pairs can give isomorphic graphs (for
    # instance through the symmetry of the 5-vertex side), so keep one job per isomorphism class.
    jobs = []
    seen = []
    A4, B5 = side_graphs(4), side_graphs(5)
    for i, ga in enumerate(A4):
        for j, gb in enumerate(B5):
            E = bipartite_plus(4, 5, ga, gb)
            g = nx.Graph(E)
            if any(nx.is_isomorphic(g, h) for h in seen):
                continue
            seen.append(g)
            jobs.append(("K45+A%d+B%d" % (i, j), E, [3]))
    F["K45AB"] = jobs

    F["K46B"] = [("K46+B%d" % j, bipartite_plus(4, 6, nx.empty_graph(4), gb), [3])
                 for j, gb in enumerate(side_graphs(6))]

    multi = []
    K4 = list(itertools.combinations(range(4), 2))
    for r in range(1, 7):
        for S in itertools.combinations(range(6), r):
            m = K4 + [K4[i] for i in S]
            multi.append(("K4 doubled %s" % (S,), subdivide_multi(4, m), [3, 4, 5]))
    for S in [(0,), (0, 1), (0, 1, 2), (0, 1, 2, 3, 4, 5)]:
        m = K4 + [K4[i] for i in S] * 2
        multi.append(("K4 tripled %s" % (S,), subdivide_multi(4, m), [3, 4]))
    K5 = list(itertools.combinations(range(5), 2))
    for r in (1, 2, 3, 5, 10):
        m = K5 + K5[:r]
        multi.append(("K5 doubled first %d" % r, subdivide_multi(5, m), [3, 4]))
    K33 = [(a, b) for a in range(3) for b in range(3, 6)]
    for r in (1, 3, 9):
        m = K33 + K33[:r]
        multi.append(("K33 doubled first %d" % r, subdivide_multi(6, m), [3, 4]))
    for t in range(3, 8):
        multi.append(("theta multigraph %d" % t, subdivide_multi(2, [(0, 1)] * t), [3, 4, 5]))
    # remove isomorphic repeats
    uniq, seen = [], []
    for name, E, ks in multi:
        g = nx.Graph(E)
        if any(nx.is_isomorphic(g, h) for h in seen):
            continue
        seen.append(g)
        uniq.append((name, E, ks))
    F["MULTI"] = uniq

    F["K55B"] = [("K55+B%d" % j, bipartite_plus(5, 5, nx.empty_graph(5), gb), [3, 4])
                 for j, gb in enumerate(B5)]

    rng = random.Random(2026)
    r9 = geng(9)
    F["R9"] = [("R9 %s" % s, list(from_g6(s).edges()), [3]) for s in rng.sample(r9, sample)]

    F["P8"] = [("P8 %s" % s, list(from_g6(s).edges()), [4, 5]) for s in planar_g6(8)]
    F["P9"] = [("P9 %s" % s, list(from_g6(s).edges()), [4]) for s in planar_g6(9)]
    p10 = planar_g6(10)
    F["P10"] = [("P10 %s" % s, list(from_g6(s).edges()), [4]) for s in rng.sample(p10, sample)]

    # ---- second phase: larger bipartite cores with graphs on their sides, and dense graphs
    # At k = 4 only the K45AB graphs whose H_2 at k = 3 has torsion other than Z/2 are rerun
    # (edge stabilization is injective, so whatever torsion they have at k = 3 persists; the
    # question is whether larger orders appear).  Requires the K45AB results in OUT.
    exotic = set()
    if os.path.exists(OUT):
        for line in open(OUT):
            rec = json.loads(line)
            if rec["family"] == "K45AB" and rec.get("summary") and \
                    any(q != "2" for q in re.findall(r"Z/(\d+)", rec["summary"])):
                exotic.add(rec["name"])
    F["K45AB4"] = [(name, E, [4]) for name, E, _ in F["K45AB"] if name in exotic]
    F["K55B4"] = [(name, E, [4]) for name, E, _ in F["K55B"]]
    B6 = side_graphs(6)
    F["K56B"] = [("K56+B%d" % j, bipartite_plus(5, 6, nx.empty_graph(5), gb), [3])
                 for j, gb in enumerate(B6)]
    B7 = [nx.convert_node_labels_to_integers(g) for g in nx.graph_atlas_g()
          if g.number_of_nodes() == 7]
    rng2 = random.Random(7)
    sample7 = rng2.sample(range(len(B7)), min(sample // 4, len(B7)))
    F["K57B"] = [("K57+B7_%d" % j, bipartite_plus(5, 7, nx.empty_graph(5), B7[j]), [3])
                 for j in sample7]
    F["K47B"] = [("K47+B7_%d" % j, bipartite_plus(4, 7, nx.empty_graph(4), B7[j]), [3])
                 for j in sample7]
    mp_jobs = []
    for parts in [(3, 3, 3), (2, 2, 2, 2), (2, 2, 2, 2, 2), (3, 3, 2), (4, 4, 1), (4, 4, 2),
                  (3, 3, 3, 1), (2, 3, 4), (1, 4, 5), (2, 4, 5)]:
        g = nx.complete_multipartite_graph(*parts)
        mp_jobs.append(("K%s" % ",".join(map(str, parts)), list(g.edges()), [3, 4]))
    F["MP"] = mp_jobs
    dense = []
    for i in range(sample // 2):
        while True:
            g = nx.gnm_random_graph(10, rng2.randint(24, 34), seed=rng2.randrange(10 ** 9))
            if nx.is_connected(g) and min(d for _, d in g.degree()) >= 3:
                break
        dense.append(("R10dense %d" % i, list(g.edges()), [3]))
    F["R10D"] = dense
    return F


# ---------------------------------------------------------------------------------- the runner

def work(job, q):
    """Run one job in a child process and put the record on the queue.

    job = (family, name, edge list, k).  Each job runs in its own process so that the parent
    can kill it after --timeout and so that its memory is returned to the system when it ends.
    """
    family, name, edges, k = job
    t0 = time.time()
    # Relabel the vertices 0..N-1 (subdivision vertices and graph6 labels need not be).
    g = nx.convert_node_labels_to_integers(nx.Graph(edges))
    G = Graph(g.number_of_nodes(), sorted((min(u, v), max(u, v)) for u, v in g.edges()))
    tors, r = h2_torsion(G, k)
    q.put({"family": family, "name": name, "k": k, "N": G.N, "M": G.M,
           "planar": nx.check_planarity(g)[0], "edges": [list(e) for e in G.edges],
           "rank": r, "primes_checked": "<= 97", "torsion": {str(p): st for p, st in tors.items()},
           "summary": summary(tors), "seconds": round(time.time() - t0, 1)})


def main():
    """Build the job list, skip stored jobs, run the rest in parallel with a timeout each."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", nargs="*", default=None)
    ap.add_argument("--sample", type=int, default=600)
    ap.add_argument("--cap", type=int, default=400_000)
    ap.add_argument("--timeout", type=int, default=600)
    ap.add_argument("--workers", type=int, default=2)
    args = ap.parse_args()

    # Resume: every (job name, k) already in the data file is skipped, and earlier hits are kept
    # for the final list.
    done = set()
    hits = []
    if os.path.exists(OUT):
        for line in open(OUT):
            rec = json.loads(line)
            done.add((rec["name"], rec["k"]))
            if rec.get("summary") and strange({int(p): st for p, st in rec["torsion"].items()},
                                              rec["planar"]):
                hits.append(rec)
    F = families(args.sample)
    order = args.family or ["K45AB", "K46B", "MULTI", "K55B", "R9", "P8", "P9", "P10"]
    queue = []
    for fam in order:
        for name, E, ks in F[fam]:
            for k in ks:
                if (name, k) in done:
                    continue
                g = nx.Graph(E)
                G = Graph(g.number_of_nodes(), [(min(u, v), max(u, v)) for u, v in
                          nx.convert_node_labels_to_integers(g).edges()])
                # The closed formula gives the size of the two slices before anything is built.
                s = reduced_cell_count_by_degree(G, k)
                cells = s.get(2, 0) + s.get(3, 0)
                if cells > args.cap:
                    continue
                queue.append((fam, name, E, k, cells))
    print("%d jobs to run" % len(queue), flush=True)

    # A small scheduler: keep --workers child processes busy, collect their records from the
    # queue as they arrive (writing each to the data file at once, so nothing is lost if the
    # search is interrupted), and kill any child that exceeds --timeout.
    q = mp.Queue()
    running = {}
    out = open(OUT, "a")
    qi = 0
    finished = 0
    while qi < len(queue) or running:
        while qi < len(queue) and len(running) < args.workers:
            fam, name, E, k, cells = queue[qi]
            qi += 1
            p = mp.Process(target=work, args=((fam, name, E, k), q))
            p.start()
            running[p.pid] = (p, time.time(), fam, name, E, k, cells)
        time.sleep(0.05)
        while not q.empty():
            rec = q.get()
            out.write(json.dumps(rec) + "\n")
            out.flush()
            finished += 1
            hit = strange({int(p): st for p, st in rec["torsion"].items()}, rec["planar"])
            if hit:
                hits.append(rec)
            print("%s %-8s %-34s k=%d N=%2d M=%2d %s H_2: %-28s [%ss]"
                  % ("!!" if hit else "  ", rec["family"], rec["name"][:34], rec["k"], rec["N"],
                     rec["M"], "planar" if rec["planar"] else "      ", rec["summary"],
                     rec["seconds"]), flush=True)
        for pid, (p, t0, fam, name, E, k, cells) in list(running.items()):
            if not p.is_alive():
                p.join()
                del running[pid]
            elif time.time() - t0 > args.timeout:
                p.kill()
                p.join()
                del running[pid]
                g = nx.Graph(E)
                rec = {"family": fam, "name": name, "k": k, "N": g.number_of_nodes(),
                       "M": g.number_of_edges(), "planar": nx.check_planarity(g)[0],
                       "edges": [list(e) for e in E], "cells": cells, "torsion": {},
                       "summary": None, "timed_out": args.timeout}
                out.write(json.dumps(rec) + "\n")
                out.flush()
                print("   %-8s %-34s k=%d  TIMED OUT after %ds" % (fam, name[:34], k, args.timeout),
                      flush=True)
    print("\n%d jobs finished.  HITS (odd torsion in a planar graph, or a summand of order >= 4):"
          % finished)
    for rec in hits:
        print("  %-8s %-40s k=%d N=%d M=%d planar=%s  H_2: %s"
              % (rec["family"], rec["name"], rec["k"], rec["N"], rec["M"], rec["planar"],
                 rec["summary"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
