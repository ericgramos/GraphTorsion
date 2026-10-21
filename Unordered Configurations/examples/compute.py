#!/usr/bin/env python3
"""
compute.py -- the homology of B_k(Gamma) for a graph you type on the command line.

WHAT IT IS FOR
--------------
This is the shortest route from "I have a graph in mind" to "here is H_*(B_k Gamma; Z)".  It
needs no Python programming: you describe the graph either by name or by an edge list, choose
the number of particles, and the script picks the right pipeline, prints the cell counts
before it starts (so you know what you are in for), and prints the answer.

USAGE
-----
    python3 examples/compute.py GRAPH K [--degree D] [--invariant-factors]

GRAPH is one of

    a named graph      complete:7, bipartite:4,5, octahedron, dodecahedron, icosahedron,
                       tetrahedron, petersen, prism:5 (the 5-prism), antiprism:4, wheel:6,
                       theta:3,2 (3 paths of length 2)
    an edge list       "0-1 0-2 1-2 2-3"   (vertices are any integers; separate edges by
                       spaces or commas; the vertices are relabelled 0..N-1 in sorted order)

K is the number of particles (k >= 1).

    --degree D             compute only the torsion of H_D, from the two degree slices D and
                           D + 1 (uses gbt.lean.torsion_in_degree).  This reaches far larger
                           complexes than a full computation; see TUTORIAL.md, Lesson 5.
    --invariant-factors    print torsion as invariant factors (Z/2 + Z/6) instead of the
                           primary decomposition ((Z/2)^2 + Z/3), which is the default.

EXAMPLES (all from the top folder of the repository)
--------
    python3 examples/compute.py complete:5 2
        graph: 5 vertices, 10 edges;  k = 2
        cells of the fully reduced complex by degree: {0: 55, 1: 150, 2: 90}  (295 in all)
        pipeline: gbt.fast (full homology)
        H_0 = Z
        H_1 = Z^6 + Z/2
        H_2 = 0

    python3 examples/compute.py bipartite:4,5 3
        ...
        H_2 = Z^155 + (Z/2)^4 + Z/3

    python3 examples/compute.py "0-1 0-2 0-3 0-4 1-2 1-3 1-4 2-3 2-4 3-4" 3
        the same as complete:5 with three particles

    python3 examples/compute.py bipartite:5,7 4 --degree 3
        pipeline: gbt.lean.torsion_in_degree (torsion of H_3 only)
        tors H_3 = (Z/2)^14 + (Z/3)^6 + Z/4          (about 10 seconds, 0.5 GB)

HOW IT CHOOSES A PIPELINE
-------------------------
Full homology is computed with gbt.fast when the fully reduced complex has at most 600,000
cells, and with gbt.lean (a third of the memory, a little slower) up to 2,000,000 cells.  Above
that it refuses and suggests --degree.  The thresholds are the ones that fit in about 7 GB; see
README Section 5.  Every vertex of a graph built here has degree >= 1 (an edge list only
mentions endpoints), which is what the fast pipelines require.

Nothing is written to data/: this script only prints.
"""

import argparse
import os
import sys
import time

# Make "import gbt" work no matter which folder the script is started from: the repository root
# is the parent of the folder containing this file.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import gbt                                                                  # noqa: E402
from gbt import fast, fmt, lean, primary                                    # noqa: E402
from gbt.sizes import reduced_cell_count_by_degree                          # noqa: E402

FAST_LIMIT = 600_000        # reduced cells: above this gbt.fast may exhaust 7 GB
LEAN_LIMIT = 2_000_000      # above this even gbt.lean's full computation is impractical


def parse_graph(text):
    """Turn the GRAPH argument into a gbt.Graph.

    `text` is either NAME, NAME:ARGS (comma-separated integers), or an edge list such as
    "0-1 0-2 1-2".  Returns a Graph.  Raises SystemExit with a readable message if the text
    cannot be understood, so that a typo does not produce a Python traceback.
    """
    if "-" in text:                                   # an edge list: every token is u-v
        edges = []
        for token in text.replace(",", " ").split():
            try:
                u, v = (int(x) for x in token.split("-"))
            except ValueError:
                sys.exit("cannot read the edge %r: write edges as u-v with integer u, v" % token)
            if u == v:
                sys.exit("the edge %r is a loop; only simple graphs are supported" % token)
            edges.append((min(u, v), max(u, v)))
        if len(set(edges)) != len(edges):
            sys.exit("the edge list repeats an edge; only simple graphs are supported")
        return gbt.from_edges(edges)

    name, _, args = text.partition(":")
    nums = [int(x) for x in args.split(",")] if args else []
    makers = {
        "complete": lambda: gbt.complete(*nums),
        "bipartite": lambda: gbt.complete_bip(*nums),
        "prism": lambda: gbt.prism_graph(*nums),
        "antiprism": lambda: gbt.antiprism(*nums),
        "wheel": lambda: gbt.wheel(*nums),
        "theta": lambda: gbt.theta(*nums),
        "octahedron": gbt.octahedron,
        "tetrahedron": gbt.tetrahedron,
        "dodecahedron": gbt.dodecahedron,
        "icosahedron": gbt.icosahedron,
        "petersen": gbt.petersen,
    }
    if name not in makers:
        sys.exit("unknown graph %r; known names: %s, or give an edge list like \"0-1 1-2\""
                 % (name, ", ".join(sorted(makers))))
    try:
        return makers[name]()
    except TypeError:
        sys.exit("wrong number of parameters for %r (for example complete:7, bipartite:4,5)"
                 % name)


def show_torsion(tors, invariant_factors):
    """The torsion part of a group as text: '(Z/2)^4 + Z/3' (primary decomposition, the
    default) or 'Z/2 + Z/2 + Z/2 + Z/6' style invariant factors when asked; '0' if none."""
    if not tors:
        return "0"
    if invariant_factors:
        return fmt(0, tors)
    return primary(tors)


def main():
    """Parse the arguments, report the size of the complex, run the chosen pipeline, print."""
    ap = argparse.ArgumentParser(description="Homology of the unordered configuration space "
                                             "B_k of a graph.")
    ap.add_argument("graph", help='a name such as complete:7 or an edge list "0-1 1-2"')
    ap.add_argument("k", type=int, help="number of particles")
    ap.add_argument("--degree", type=int, default=None,
                    help="compute only the torsion of this homology degree")
    ap.add_argument("--invariant-factors", action="store_true",
                    help="print torsion as invariant factors instead of primary summands")
    args = ap.parse_args()

    G = parse_graph(args.graph)
    k = args.k
    if k < 1:
        sys.exit("k must be at least 1")
    print("graph: %d vertices, %d edges;  k = %d" % (G.N, G.M, k))

    by_degree = reduced_cell_count_by_degree(G, k)
    total = sum(by_degree.values())
    print("cells of the fully reduced complex by degree: %s  (%s in all)"
          % (by_degree, format(total, ",")))

    t0 = time.time()
    if args.degree is not None:
        d = args.degree
        slice_cells = by_degree.get(d, 0) + by_degree.get(d + 1, 0)
        print("pipeline: gbt.lean.torsion_in_degree (torsion of H_%d only; %s cells in the "
              "two slices)" % (d, format(slice_cells, ",")))
        tors = lean.torsion_in_degree(G, k, d)
        print("tors H_%d = %s" % (d, show_torsion(tors, args.invariant_factors)))
    else:
        if total <= FAST_LIMIT:
            print("pipeline: gbt.fast (full homology)")
            H, _, _ = fast.compute(G, k)
        elif total <= LEAN_LIMIT:
            print("pipeline: gbt.lean (full homology, memory-lean)")
            H, _, _ = lean.compute(G, k)
        else:
            sys.exit("the complex has %s cells, too many for a full computation; use "
                     "--degree D to get the torsion of a single H_D" % format(total, ","))
        for deg in sorted(H):
            free, tors = H[deg]
            if args.invariant_factors:
                text = fmt(free, tors)
            else:
                parts = (["Z^%d" % free if free > 1 else "Z"] if free else []) + \
                        ([primary(tors)] if tors else [])
                text = " + ".join(parts) if parts else "0"
            print("H_%d = %s" % (deg, text))
    print("(%.1f s)" % (time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
