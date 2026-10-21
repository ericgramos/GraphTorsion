#!/usr/bin/env python3
"""
Indecomposable torsion above H_2 in a planar graph: a search where products are impossible.

THE QUESTION
------------
higher_degree_torsion.py produces 2-torsion in H_i of a planar graph for every i >= 2, but it
does so by taking products with disjoint cycles, so the classes are Kunneth classes -- a
torsion class in H_2 tensored with the fundamental classes of i-2 disjoint circles.  Those say
nothing that H_2 did not already say.  The question is whether a planar graph can carry torsion
in H_3 that is NOT of that form, i.e. whether degree two is the only degree in which a planar
graph CREATES torsion.

WHY THE VERTEX COUNT IS THE WHOLE POINT
---------------------------------------
The cleanest way to be sure a class is not a product is to look where no product can reach.
A product class in H_3(B_k Gamma) comes from two VERTEX-DISJOINT subgraphs G_1, G_2 with
k_1 + k_2 = k, and by Kunneth the only sources of torsion in total degree 3 are

    H_p(B_{k_1}G_1) (x) H_q(B_{k_2}G_2),  p + q = 3, with torsion in one factor, or
    Tor(H_p, H_q),                        p + q = 2.

Count vertices.  Torsion in H_1 needs a non-planar factor (Ko--Park), which a subgraph of a
planar graph cannot be; so Tor(H_1, H_1) and every term with torsion in an H_1 vanish.  Torsion
in H_2 needs at least 12 edges (true of every graph on at most 7 vertices, by the exhaustive
sweeps in sweep_small_graphs.py), so at least 6 vertices; and for H_2 (x) H_1 the other factor
must have nonzero H_1, so it needs a cycle: at least 3 more vertices.  (The remaining term,
H_3 (x) H_0, is torsion in H_3 of a smaller planar subgraph, which the search would itself find.)
Hence

    any planar graph on at most 8 vertices whose H_3 has torsion has H_3 torsion that is not a
    Kunneth class built from H_2.

THE POOL
--------
  * all connected planar graphs on 6 and 7 vertices with minimum degree >= 3, from networkx's
    graph atlas (complete up to 7 vertices): 9 on 6 vertices and 46 on 7;
  * a SAMPLE of 10 eight-vertex planar graphs: 3- and 4-regular ones found by random sampling
    up to isomorphism (seeded, so the sample is reproducible), the 4-antiprism, and maximal
    planar graphs grown by stacking a vertex into a random face.  The 8-vertex part is NOT
    exhaustive.

Minimum degree 3 loses only subdivisions (which do not change the configuration spaces) and
graphs with leaves; leaves are NOT covered by this search.

THE METHOD
----------
Only the torsion of H_d (default d = 3) is computed, with `lean.torsion_in_degree`, which uses
tors H_d = tors coker d_{d+1}; so only the degree-d and degree-(d+1) slices of the fully reduced
complex are built.  The slice sizes are known in advance (`reduced_cell_count_by_degree`), and a
graph whose two slices together exceed --cap is skipped.

WHAT THE SEARCH FOUND
---------------------
Nothing.  Of the 65 graphs in the pool, those whose two slices fit under the default cap were
computed: 59 at k = 5 (all 9 + 46 on 6 and 7 vertices, 4 of the 10 on 8 vertices) and 52 at
k = 6 (9 + 40 + 3), as stored in data/genuine_higher_torsion.json.  Every one has torsion-free
H_3.  The octahedron has torsion-free H_3 at k = 7 and H_4 at k = 6 and 7, and the
4-antiprism torsion-free H_3 at k = 6 and H_4 at k = 5 (`--deep`).  Every planar example of
torsion above H_2 that we have is still a Kunneth class.

There is a conceptual reason to expect this.  An--Knudsen's three generators of planar H_2 all
have homological dimension at most two: a cycle and a star have free fundamental group and no
homology above degree one, and the theta graph has two essential vertices, so Swiatkowski's
dimension bound gives H_i(B_k Theta) = 0 for i >= 3.  If their generation statement extends to
higher degrees with the same atoms, then every class in H_i for i >= 3 is built from atom
classes living in degrees <= 2, hence is a product -- and so is any torsion in it.  Degree two
would then be the only degree in which a planar graph creates torsion, which is the conjecture
this search supports.

WHAT IS PRINTED
---------------
A header giving the pool sizes, then one row per tested graph:

    #             position in the pool (0-8 six vertices, 9-54 seven, 55-64 eight)
    N, M          vertices and edges
    girth         length of the shortest cycle
    slice cells   cells in the degree-d and degree-(d+1) slices together
    H_d torsion   the torsion of H_d as a primary decomposition, or "none"
    time          seconds, or "cached" if the result was already stored

A row with torsion is flagged "<== INDECOMPOSABLE".  Graphs over the cap are skipped silently
and counted.  The run ends with "T graphs tested, S skipped as too large at this cap." and
either the list of hits with their edge lists or "No tested planar graph has torsion in H_d at
k = K."  Excerpt (default options, all cached, 3 seconds):

    65 planar graphs with minimum degree >= 3: 9 on 6, 46 on 7, 10 on 8 vertices
       #    N    M   girth  slice cells H_3 torsion    time
       0    6    9       3        9,360 none           cached
    ...
    59 graphs tested, 6 skipped as too large at this cap.
    No tested planar graph has torsion in H_3 at k = 5.

With --deep the table has columns graph, N, M, k, deg, slice cells, torsion, time, for the
fixed jobs octahedron (k, d) = (7, 3), (6, 4), (7, 4) and 4-antiprism (6, 3), (5, 4); a job
whose two slices exceed 2,200,000 cells is skipped.

USAGE
-----
    python3 experiments/genuine_higher_torsion.py [--k 5] [--degree 3] [--cap 400000]
    python3 experiments/genuine_higher_torsion.py --deep

    --k K        (default 5) number of particles.
    --degree D   (default 3) the homological degree whose torsion is computed.  The vertex
                 argument above is written for D = 3; for D = 4 a Kunneth class needs even more
                 vertices (an H_2 class and two disjoint cycles: 12), so it still applies.
    --cap C      (default 400000) skip a graph whose two slices have more than C cells.
    --deep       instead of the sweep, run the fixed deep jobs listed above.

Results are stored in data/genuine_higher_torsion.json.  Sweep records have key
"N|u.v-u.v-...|K|D" (vertex count, sorted edge list, particles, degree) and value
{"torsion": [...], "N", "M", "k", "degree", "edgelist"}; deep records have key
"deep|NAME|K|D" and no edge list.  Stored results are reused, so a rerun with the shipped
data takes seconds.  From scratch, k = 5 takes a few minutes and k = 6 an hour or so; --deep
reaches about a million cells per slice and takes a few minutes per job and a few GB.  The
closing line "wrote data/..." is printed even when every result was cached.

Exit code: always 0 -- the finding is the output, not the status.
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
    """All connected planar graphs on `nvert` vertices with minimum degree >= mindeg.

    Read from networkx's graph atlas, which lists every graph on at most 7 vertices exactly
    once up to isomorphism; hence the ValueError for nvert > 7.  Returns networkx graphs.
    """
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
        """Append Gx to the sample if connected, min degree >= 3, planar, and new up to iso."""
        if not nx.is_connected(Gx) or min(d for _, d in Gx.degree()) < 3:
            return
        if not nx.check_planarity(Gx)[0]:
            return
        if any(nx.is_isomorphic(Gx, H) for H in out):
            return
        out.append(Gx)

    rng = random.Random(0)          # fixed seed: the same sample on every run
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

    Both have at most 8 vertices, so nothing here can be a product.  Jobs are fixed (see the
    module docstring); each result is stored under "deep|NAME|K|D".  Returns 0.
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
    """Build the pool, compute the torsion of H_d for each graph under the cap, report hits."""
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
        # key: vertex count | sorted edge list as "u.v-u.v-..." | k | degree
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
