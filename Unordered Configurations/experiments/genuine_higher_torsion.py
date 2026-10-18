#!/usr/bin/env python3
"""
Indecomposable torsion above H_2 in a planar graph: a search where products are impossible.

WHY THE VERTEX COUNT IS THE WHOLE POINT
---------------------------------------
higher_degree_torsion.py produces 2-torsion in H_i of a planar graph for every i >= 2, but it does it by
taking products with disjoint cycles, so the classes are K\"unneth classes -- a torsion class in
H_2 tensored with the fundamental classes of i-2 disjoint circles.  Those are not interesting:
they say nothing that H_2 did not already say.  The question is whether a planar graph can carry
torsion in H_3 that is NOT of that form.

The cleanest way to be sure is to look where no product can reach.  A product class in
H_3(B_k Gamma) comes from two VERTEX-DISJOINT subgraphs G_1, G_2 with k_1 + k_2 = k, and by
K\"unneth the only sources of torsion in total degree 3 are

    H_p(B_{k_1}G_1) (x) H_q(B_{k_2}G_2),  p + q = 3, with torsion in one factor, or
    Tor(H_p, H_q),                        p + q = 2.

Now count vertices.  Torsion in H_1 needs a non-planar factor (Ko--Park), which a subgraph of a
planar graph cannot be.  Torsion in H_2 needs at least 12 edges (true of every graph on at most
7 vertices, by the exhaustive sweeps in sweep_small_graphs.py), so at least 6 vertices.  And the other factor must have nonzero H_q with
q >= 1, so it needs a cycle: at least 3 more vertices.  Hence

    any planar graph on at most 8 vertices whose H_3 has torsion has INDECOMPOSABLE H_3 torsion.

So the search below is exhaustive over the planar graphs on 6 and 7 vertices with minimum degree
3 -- all 55 of them, from the graph atlas -- together with a SAMPLE of 10 eight-vertex planar
graphs (cubic ones found by random sampling, the 4-antiprism, some maximal planar ones); the
eight-vertex part is not exhaustive.  Minimum degree 3
loses only subdivisions and graphs with leaves; leaves are NOT covered by this search.

Only the torsion of H_3 is computed, via tors H_d = tors coker d_{d+1} and `lean.torsion_in_degree`, so only the
degree-3 and degree-4 slices of each complex are built.

WHAT THE SEARCH FOUND
---------------------
Nothing.  Of the 65 graphs in the pool, those whose two slices fit under the default cap were
computed: 59 at k = 5 (all 9 + 46 on 6 and 7 vertices, 4 of the 10 on 8 vertices) and 52 at
k = 6 (9 + 40 + 3), as stored in data/genuine_higher_torsion.json.  Every one has torsion-free
H_3.  The octahedron has torsion-free H_3 at k = 7 and H_4 at k = 6 and 7, and the
4-antiprism torsion-free H_3 at k = 6 and H_4 at k = 5 (`--deep`).  Every planar example of torsion
above H_2 that we have is still a Kunneth class.

There is a conceptual reason to expect this.  An--Knudsen's three generators of planar H_2 all
have homological dimension at most two: a cycle and a star have free fundamental group and no
homology above degree one, and the theta graph has two essential vertices, so Swiatkowski's
dimension bound gives H_i(B_k Theta) = 0 for i >= 3 (which the computation confirms through
k = 6).  If their generation statement extends to higher degrees with the same atoms, then every
class in H_i for i >= 3 is built from atom classes living in degrees <= 2, hence is a product --
and so is any torsion in it.  Degree two would then be the only degree in which a planar graph
*creates* torsion -- which is the conjecture this search supports.

    python3 experiments/genuine_higher_torsion.py [--k 5] [--cap CELLS] [--degree 3]
    python3 experiments/genuine_higher_torsion.py --deep

k = 5 is a few minutes for the whole sweep; k = 6 is an hour or so.  `--deep` pushes the
octahedron to k = 7 in degrees 3 and 4 (and the 4-antiprism to k = 6), which is the furthest this implementation
reaches (a million cells per slice, a few minutes each).  Exits 0 whether or not anything is
found -- the finding is the output, not the status.
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import (antiprism, from_networkx, lean, octahedron, paths,            # noqa: E402
                 primary)
from gbt.sizes import reduced_cell_count_by_degree as slices_by_degree         # noqa: E402

DATA_FILE = "genuine_higher_torsion.json"


def atlas_planar(nvert, mindeg=3):
    """All connected planar graphs on `nvert` vertices with minimum degree >= mindeg."""
    import networkx as nx
    from networkx.generators.atlas import graph_atlas_g
    if nvert > 7:
        raise ValueError("the atlas stops at 7 vertices")
    out = []
    for Gx in graph_atlas_g():
        if Gx.number_of_nodes() != nvert or not nx.is_connected(Gx):
            continue
        if min(d for _, d in Gx.degree()) < mindeg:
            continue
        if not nx.check_planarity(Gx)[0]:
            continue
        out.append(Gx)
    return out


def eight_vertex_planar():
    """The 8-vertex planar graphs we can get at without a graph generator.

    Cubic and 4-regular ones by random sampling up to isomorphism, plus the 4-antiprism and
    maximal planar graphs grown by stacking.
    """
    import random
    import networkx as nx
    out = []

    def add(Gx):
        if not nx.is_connected(Gx) or min(d for _, d in Gx.degree()) < 3:
            return
        if not nx.check_planarity(Gx)[0]:
            return
        if any(nx.is_isomorphic(Gx, H) for H in out):
            return
        out.append(Gx)

    rng = random.Random(0)
    for _ in range(3000):
        add(nx.random_regular_graph(3, 8, seed=rng.randrange(10 ** 9)))
    for _ in range(3000):
        add(nx.random_regular_graph(4, 8, seed=rng.randrange(10 ** 9)))
    add(nx.Graph(antiprism(4).edges))
    # maximal planar: stack a vertex into a random face of a triangulation
    for seed in range(40):
        G = nx.complete_graph(4)
        r = random.Random(seed)
        faces = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
        for v in range(4, 8):
            f = faces.pop(r.randrange(len(faces)))
            for x in f:
                G.add_edge(v, x)
            faces += [(v, f[0], f[1]), (v, f[0], f[2]), (v, f[1], f[2])]
        add(G)
    return out


def deep_probe():
    """The octahedron and the 4-antiprism as far as the machine reaches: degrees 3 and 4.

    Both have at most 8 vertices, so nothing here can be a product.
    """
    res = paths.load_json(DATA_FILE)
    jobs = []
    for name, G in (("octahedron", octahedron()),):
        for k, d in ((7, 3), (6, 4), (7, 4)):
            jobs.append((name, G, k, d))
    jobs.append(("antiprism_4", antiprism(4), 6, 3))
    jobs.append(("antiprism_4", antiprism(4), 5, 4))

    hdr = ("%-14s %3s %3s %4s %4s %13s %-14s %s"
           % ("graph", "N", "M", "k", "deg", "slice cells", "torsion", "time"))
    print(hdr)
    print("-" * len(hdr))
    found = []
    for name, G, k, d in jobs:
        sl = slices_by_degree(G, k)
        cells = sl.get(d, 0) + sl.get(d + 1, 0)
        if cells > 2_200_000:
            print("%-14s %3d %3d %4d %4d %13s  (skipped: too large)"
                  % (name, G.N, G.M, k, d, format(cells, ",")))
            continue
        key = "deep|%s|%d|%d" % (name, k, d)
        if key in res:
            t = res[key]["torsion"]
            took = "cached"
        else:
            t0 = time.time()
            t = lean.torsion_in_degree(G, k, d)
            took = "%.0fs" % (time.time() - t0)
            res[key] = {"torsion": t, "N": G.N, "M": G.M, "k": k, "degree": d}
            paths.save_json(DATA_FILE, res)
        if t:
            found.append((name, k, d, t))
        print("%-14s %3d %3d %4d %4d %13s %-14s %s%s"
              % (name, G.N, G.M, k, d, format(cells, ","),
                 primary(t) if t else "none", took,
                 "   <== INDECOMPOSABLE" if t else ""), flush=True)
    print("-" * len(hdr))
    if found:
        print("\nindecomposable higher torsion: %s" % (found,))
    else:
        print("\nThe octahedron has torsion-free H_3 and H_4 at every k <= 7.")
        print("Its H_2 torsion does not propagate upward on its own: the higher-degree")
        print("classes that do carry torsion all need an extra disjoint cycle.")
    print("\nwrote data/%s" % DATA_FILE)
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--degree", type=int, default=3)
    ap.add_argument("--cap", type=int, default=400_000,
                    help="skip a graph if its two relevant slices exceed this")
    ap.add_argument("--deep", action="store_true",
                    help="instead: push the octahedron to k=7 and the 4-antiprism to k=6, degrees 3 and 4")
    args = ap.parse_args()
    k, d = args.k, args.degree

    if args.deep:
        return deep_probe()

    import networkx as nx
    pool = []
    for n in (6, 7):
        for Gx in atlas_planar(n):
            pool.append((n, Gx))
    for Gx in eight_vertex_planar():
        pool.append((8, Gx))
    print("%d planar graphs with minimum degree >= 3: %s on 6, %s on 7, %s on 8 vertices"
          % (len(pool), sum(1 for n, _ in pool if n == 6),
             sum(1 for n, _ in pool if n == 7),
             sum(1 for n, _ in pool if n == 8)))
    print("looking for torsion in H_%d at k = %d, where no Kunneth product can contribute"
          % (d, k))
    print("(a product would need at least 9 vertices -- see the docstring)\n")

    res = paths.load_json(DATA_FILE)
    hdr = "%4s %4s %4s %7s %12s %-14s %s" % ("#", "N", "M", "girth", "slice cells",
                                             "H_%d torsion" % d, "time")
    print(hdr)
    print("-" * len(hdr))
    hits = []
    skipped = 0
    for i, (n, Gx) in enumerate(pool):
        G = from_networkx(Gx)
        sl = slices_by_degree(G, k)
        cells = sl.get(d, 0) + sl.get(d + 1, 0)
        if cells > args.cap:
            skipped += 1
            continue
        key = "%d|%s|%d|%d" % (n, "-".join("%d.%d" % e for e in sorted(G.edges)), k, d)
        if key in res:
            t = res[key]["torsion"]
            took = "cached"
        else:
            t0 = time.time()
            t = lean.torsion_in_degree(G, k, d)
            took = "%.0fs" % (time.time() - t0)
            res[key] = {"torsion": t, "N": G.N, "M": G.M, "k": k, "degree": d,
                        "edgelist": sorted(tuple(e) for e in G.edges)}
            paths.save_json(DATA_FILE, res)
        if t:
            hits.append((n, G.M, sorted(tuple(e) for e in G.edges), t))
        print("%4d %4d %4d %7d %12s %-14s %s%s"
              % (i, G.N, G.M, nx.girth(Gx), format(cells, ","),
                 primary(t) if t else "none", took,
                 "   <== INDECOMPOSABLE" if t else ""), flush=True)
    print("-" * len(hdr))

    print("\n%d graphs tested, %d skipped as too large at this cap." % (len(pool) - skipped,
                                                                       skipped))
    if hits:
        print("\n%d planar graph(s) on at most 8 vertices with torsion in H_%d at k=%d:"
              % (len(hits), d, k))
        for n, m, edges, t in hits:
            print("   N=%d, M=%d, torsion %s" % (n, m, primary(t)))
            print("      edges %s" % (edges,))
        print("\nNo product can produce these, so this is INDECOMPOSABLE torsion above H_2")
        print("in a planar graph -- not a Kunneth class built out of H_2.")
    else:
        print("\nNo tested planar graph has torsion in H_%d at k = %d." % (d, k))
        if skipped:
            print("(%d graph(s) were skipped at this cap, so the search is not exhaustive;"
                  " raise --cap to include them.)" % skipped)
        print("Every planar example of higher torsion we have is still a product class.")
    print("\nwrote data/%s" % DATA_FILE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
