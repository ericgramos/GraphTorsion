#!/usr/bin/env python3
"""
Can adding edges kill 2-torsion?  Random supergraph walks from known carriers.

THE QUESTION
------------
Is "H_i(B_k Gamma) has no 2-torsion" (i >= 2, k fixed) a minor-closed property?  Equivalently:
if a graph H carries 2-torsion in H_i(B_k H), does every graph containing H as a minor carry
2-torsion in H_i(B_k) too?  On graphs with at most 7 vertices this holds exhaustively (every
single edge deletion and contraction was checked against the stored small-graph sweeps).  For
odd torsion the analogous statement is FALSE: the Z/3 of K_{4,5} dies in K_9
(density_interpolation.py).  This script looks for the first failure for 2-torsion beyond seven
vertices.

THE METHOD
----------
Start from a known carrier H and k.  Repeatedly enlarge the graph at random, keeping H as a
subgraph (so H stays a minor of every graph on the walk):

    * with probability 3/4, add a random edge between two existing non-adjacent vertices;
    * with probability 1/4 (while the graph has fewer than --nmax vertices), add a NEW vertex
      joined to 1, 2 or 3 random existing vertices (a pendant edge when it is 1).

After each enlargement compute the torsion of H_2(B_k) with `lean.torsion_in_degree` (only the
degree-2 and degree-3 slices are built).  If it has no 2-primary summand, compute the torsion of
H_3, ..., H_k as well; if none of them has 2-torsion either, the graph is a COUNTEREXAMPLE to
minor-closedness, and it is printed with its full edge list.  A walk stops after --steps
enlargements, when the graph becomes complete with --nmax vertices, or when the two slices would
exceed --cap cells.

The starting carriers (name: graph, k):

    octahedron       K_{2,2,2}, k = 4         Z/2 in H_2; 6 vertices, 12 edges
    c7a, c7b, c7c    the three 7-vertex, 12-edge carriers at k = 4 (read from
                     data/small_graphs_n7_k4_d1.json), each with Z/2 in H_2
    leaf             K_5 + a vertex joined to two of its vertices + a pendant edge, k = 4
    antiprism_4      k = 4, 8 vertices, 16 edges
    prism_5          k = 4, 10 vertices, 15 edges
    triangulation_9  the stacked maximal planar graph on 9 vertices of planar_odd_search.py, k = 4
    K_{4,5}          k = 3; carries (Z/2)^4 + Z/3, and only the 2-primary part is tracked
    K_7              k = 3, Z/2 in H_2

OUTPUT
------
One line per step:

    octahedron  walk 0 step 3  +edge (2,7)      N= 8 M=16  slices   62,208  H_2: Z/2        [1.2s]

where "+edge (u,v)" or "+vertex 8 -> [1, 5]" says what was added, N and M are the vertex and
edge counts, "slices" is the number of cells in the degree-2 and degree-3 slices, and the
torsion of H_2 is printed as a primary decomposition.  A line marked "<== NO 2-TORSION IN H_2"
is followed by the check of the higher degrees and, if that also fails, by
"COUNTEREXAMPLE" and the edge list.  The run ends with a summary per start: the number of
graphs computed, the largest one, and the number of counterexamples.

USAGE
-----
    python3 experiments/supergraph_walks.py [--start NAME ...] [--walks 4] [--steps 8]
                                            [--cap 600000] [--nmax 10] [--seed 0]

    --start   which carriers to start from (default: all of them, in the order above)
    --walks   random walks per start (default 4)
    --steps   enlargements per walk (default 8)
    --cap     skip a graph whose degree-2 and degree-3 slices together exceed this many cells
              (default 600,000: at most a minute or two each, under 2 GB)
    --nmax    largest number of vertices (default 10)
    --seed    seed of the random generator; walk w of start s uses seed (seed, s, w), so runs
              are reproducible and a rerun with more --walks extends the earlier walks

The results are stored in data/supergraph_walks.json, key "START|walk|step", record
{"k": k, "N": N, "M": M, "edges": [[u, v], ...], "slices": cells, "torsion": {"2": [...],
...}} (invariant factors of the torsion of each degree computed).  Stored graphs are not
recomputed.  Exit code: 1 if a counterexample was found, 0 otherwise.  The default run takes
about half an hour; the stored run (--walks 6 --steps 10) about two hours.

WHAT THE STORED RUN FOUND
-------------------------
470 graphs from nine of the ten starts, up to 10 vertices and 32 edges.  The 2-torsion of H_2
never died: no counterexample.  The number of Z/2 summands goes up and down along a walk (for
instance (Z/2)^5 -> (Z/2)^3 -> (Z/2)^2 as edges are added) but never reaches zero.  The walk
K_{4,5} walk 3 found the first summand of order 4 in H_2: after the edges (5,6), (4,8), (5,8),
(4,5) were added on the 5-vertex side,

    K_{4,5}          walk 3 step 4  +edge (4,5)  N= 9 M=24  ...  H_2: (Z/2)^14 + Z/3 + (Z/4)^2

which led to the systematic search torsion_hunt.py.  The run was stopped at K_{4,5} walk 3 step
7, where the exact Smith form inside lean.torsion_in_degree stalled (see PITFALLS), so the
K_{4,5} walks 4-5 and the K_7 walks are not in the stored data.

PITFALLS
--------
* Only the stated k is tested.  Torsion persists as k grows (edge stabilization is injective,
  An--Drummond-Cole--Knudsen, Subdivisional spaces and graph braid groups, Proposition 5.21), so
  a counterexample at this k would remain one only if the bigger graph never acquires 2-torsion
  at larger k; the script reports the fixed-k statement.
* Random walks are a search, not a proof: "no counterexample" means none was met.
* The torsion is computed with lean.torsion_in_degree, whose final exact Smith form can stall
  for hours on a tiny residual with growing entries (it did on K_{4,5} plus a triangle).  If a
  step hangs, stop the run (by process number) and check that graph with the local Smith form,
  gbt.linalg.local_smith, as torsion_hunt.py does; it handles such cases in a fraction of a
  second.
"""

import argparse
import json
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import (Graph, antiprism, complete, complete_bip, from_edges, lean,       # noqa: E402
                 octahedron, paths, primary, prism_graph)
from gbt.sizes import reduced_cell_count_by_degree                                 # noqa: E402

DATA_FILE = "supergraph_walks.json"

#: the leaf example: K_5 on 0..4, vertex 5 joined to 3 and 4, pendant edge 5-6
LEAF_EDGES = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4),
              (3, 5), (4, 5), (5, 6)]


def stacked_triangulation(n, seed):
    """The stacked maximal planar graph of planar_odd_search.py: K_4, then each new vertex is
    put inside a face chosen with random.Random(seed) and joined to its three corners."""
    edges = {(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)}
    rng = random.Random(seed)
    faces = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    for v in range(4, n):
        f = faces.pop(rng.randrange(len(faces)))
        for x in f:
            edges.add((min(v, x), max(v, x)))
        faces += [(v, f[0], f[1]), (v, f[0], f[2]), (v, f[1], f[2])]
    return Graph(n, sorted(edges))


def starts():
    """[(name, edge list, k)] of the starting carriers, in the order documented above."""
    out = [("octahedron", list(octahedron().edges), 4)]
    with open(os.path.join(paths.DATA, "small_graphs_n7_k4_d1.json")) as fh:
        rows = json.load(fh)
    c12 = [r for r in rows if r["torsion_above_H1"] and r["edges"] == 12]
    for name, r in zip(("c7a", "c7b", "c7c"), c12):
        out.append((name, [tuple(e) for e in r["edgelist"]], 4))
    out.append(("leaf", LEAF_EDGES, 4))
    out.append(("antiprism_4", list(antiprism(4).edges), 4))
    out.append(("prism_5", list(prism_graph(5).edges), 4))
    # planar_odd_search.py builds triangulation_9 with seed 2 (see its families())
    out.append(("triangulation_9", list(stacked_triangulation(9, 2).edges), 4))
    out.append(("K_{4,5}", list(complete_bip(4, 5).edges), 3))
    out.append(("K_7", list(complete(7).edges), 3))
    return out


def two_torsion(factors):
    """True if the invariant factors contain an even number, i.e. there is 2-torsion."""
    return any(f % 2 == 0 for f in factors)


def slice_cells(G, k, d):
    """Cells in the degree-d and degree-(d+1) slices of the fully reduced complex."""
    s = reduced_cell_count_by_degree(G, k)
    return s.get(d, 0) + s.get(d + 1, 0)


def enlarge(edges, N, nmax, rng):
    """One random enlargement of the graph with vertices 0..N-1 and the given edge set.

    Returns (new edge set, new N, description) or None if the graph is already complete on
    nmax vertices.  New vertices get the next free number, so the starting graph keeps its
    labels and stays a subgraph.
    """
    present = set(edges)
    missing = [(u, v) for u in range(N) for v in range(u + 1, N) if (u, v) not in present]
    if N < nmax and (not missing or rng.random() < 0.25):
        r = rng.choice((1, 2, 3))
        nbrs = sorted(rng.sample(range(N), min(r, N)))
        return present | {(x, N) for x in nbrs}, N + 1, "+vertex %d -> %s" % (N, nbrs)
    if not missing:
        return None
    e = rng.choice(missing)
    return present | {e}, N, "+edge (%d,%d)" % e


def main():
    """Run the walks, store every graph computed, and report counterexamples."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", nargs="*", default=None)
    ap.add_argument("--walks", type=int, default=4)
    ap.add_argument("--steps", type=int, default=8)
    ap.add_argument("--cap", type=int, default=600_000)
    ap.add_argument("--nmax", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    res = paths.load_json(DATA_FILE)       # stored steps are reused, so a rerun is cheap
    found = []
    summary = []
    for si, (name, edges0, k) in enumerate(starts()):
        if args.start and name not in args.start:
            continue
        computed, largest = 0, (0, 0)
        for w in range(args.walks):
            # One generator per (seed, start, walk): walk w is the same in every run, so a rerun
            # with more --walks or --steps extends the stored walks instead of replacing them.
            rng = random.Random("%d|%d|%d" % (args.seed, si, w))
            edges = {tuple(sorted(e)) for e in edges0}
            N = 1 + max(max(e) for e in edges)
            for step in range(1, args.steps + 1):
                nxt = enlarge(edges, N, args.nmax, rng)
                if nxt is None:
                    break
                edges, N, what = nxt
                G = Graph(N, sorted(edges))
                cells = slice_cells(G, k, 2)
                if cells > args.cap:
                    print("%-16s walk %d step %d  %-18s N=%2d M=%2d  slices %9s  (over the cap;"
                          " walk ends)" % (name, w, step, what, N, G.M, format(cells, ",")))
                    break
                key = "%s|%d|%d" % (name, w, step)
                rec = res.get(key)
                # Recompute unless the stored record is for exactly this graph (it always is
                # when the seed is unchanged; the check guards against a changed --seed).
                if rec is None or rec["edges"] != [list(e) for e in sorted(edges)]:
                    t0 = time.time()
                    tors = {"2": lean.torsion_in_degree(G, k, 2)}
                    rec = {"k": k, "N": N, "M": G.M, "edges": [list(e) for e in sorted(edges)],
                           "slices": cells, "torsion": tors, "seconds": round(time.time() - t0, 1)}
                    res[key] = rec
                    paths.save_json(DATA_FILE, res)
                computed += 1
                largest = max(largest, (N, G.M))
                t2 = rec["torsion"]["2"]
                flag = "" if two_torsion(t2) else "   <== NO 2-TORSION IN H_2"
                print("%-16s walk %d step %d  %-18s N=%2d M=%2d  slices %9s  H_2: %-16s [%ss]%s"
                      % (name, w, step, what, N, G.M, format(cells, ","),
                         primary(t2) if t2 else "0", rec.get("seconds", "cached"), flag),
                      flush=True)
                if two_torsion(t2):
                    continue
                # H_2 lost its 2-torsion.  The property under test is "no 2-torsion in ANY
                # degree >= 2", so look at H_3, ..., H_k before calling it a counterexample.
                elsewhere = False
                for d in range(3, k + 1):
                    if str(d) not in rec["torsion"]:
                        if slice_cells(G, k, d) > args.cap:
                            print("      H_%d: too large to check" % d)
                            elsewhere = None
                            continue
                        rec["torsion"][str(d)] = lean.torsion_in_degree(G, k, d)
                        paths.save_json(DATA_FILE, res)
                    td = rec["torsion"][str(d)]
                    print("      H_%d: %s" % (d, primary(td) if td else "0"))
                    elsewhere = elsewhere or two_torsion(td)
                if elsewhere is False:
                    print("      COUNTEREXAMPLE: no 2-torsion in any degree >= 2 at k = %d,"
                          " yet it contains %s" % (k, name))
                    print("      edges:", sorted(edges))
                    found.append((name, k, sorted(edges)))
        summary.append((name, k, computed, largest))

    print("\nsummary")
    for name, k, computed, (N, M) in summary:
        print("  %-16s k=%d  %3d graphs computed, largest N=%d M=%d" % (name, k, computed, N, M))
    print("counterexamples to minor-closedness of 'no 2-torsion above H_1': %d" % len(found))
    print("wrote data/%s" % DATA_FILE)
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
