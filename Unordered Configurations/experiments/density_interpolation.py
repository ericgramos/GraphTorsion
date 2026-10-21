#!/usr/bin/env python3
"""
Does density destroy the exotic odd torsion?  Interpolating from K_{4,5} to K_9.

THE QUESTION
------------
K_{4,5} at k=3 has H_2 = Z^155 + (Z/2)^4 + Z/3: the first odd torsion in a graph braid group.
K_9 at k=3 has H_2 = Z^1401 + (Z/2)^28: no odd torsion at all, and a much larger 2-torsion
multiplicity.  Both graphs have the same nine vertices, and the second contains the first.  So
somewhere along the way the Z/3 dies.  Where, and which edges kill it, says something about
where odd torsion comes from: is it a property of containing K_{4,5}, or of containing it with
nothing else on its large side?

THE METHOD
----------
Six graphs on the vertex set {0..8}, with A = {0,1,2,3} and B = {4,...,8}, are computed in full
at k particles (default 3), each by `gbt.fast` (or `gbt.lean` above 600,000 reduced cells):

    graph                     edges   added to K_{4,5}
    K_{4,5}                     20    nothing
    K_{4,5} + 3 in B            23    the first 3 pairs of B in lexicographic order:
                                      {4,5}, {4,6}, {4,7}  -- a star at vertex 4
    K_{4,5} + 6 in B            26    the first 6 pairs of B: all pairs containing 4 or {5,6},
                                      {5,7}
    K_{4,5} + all A             26    all 6 pairs inside A
    K_{4,5} + all B             30    all 10 pairs inside B
    K_9                         36    everything

RESULT (k = 3, stored in data/density_interpolation.json)
---------------------------------------------------------
    graph                 edges   rank H_1 torsion of H_2                 odd?     note
    K_{4,5}                  20         12 H_2: (Z/2)^4 + Z/3             [3]      ...
    K_{4,5} + 3 in B         23         15 H_2: (Z/2)^13 + Z/3            [3]      ...
    K_{4,5} + 6 in B         26         18 H_2: (Z/2)^18                  -        ...
    K_{4,5} + all A          26         18 H_2: (Z/2)^4 + Z/3             [3]      ...
    K_{4,5} + all B          30         22 H_2: (Z/2)^22                  -        ...
    K_9                      36         28 H_2: (Z/2)^28                  -        ...

Two conclusions.

  1.  Filling in the SMALL side changes nothing; filling in the LARGE side kills the Z/3
      somewhere between 3 and 6 added edges.  This fits the mechanism: the inclusion matrix
      W_{n-1,n}(b) governs the top differential of K_{n+1,b} only while the large side is
      INDEPENDENT.  Edges inside B add cells that are not matching configurations, and their
      boundaries supply the relations that divide out the Z/3.  It is the evidence behind two
      conjectures: odd torsion needs a K_{4,5} minor, and complete graphs carry only 2-primary
      torsion.

  2.  "Has odd torsion" is NOT closed under adding edges: K_9 contains K_{4,5} and has none.
      Equivalently, "has no odd torsion" is not closed under taking subgraphs (hence not under
      minors), so it cannot be characterised by a finite list of forbidden minors in the way
      Kuratowski's theorem characterises planarity.  Containing K_{4,5} is necessary for odd
      torsion (conjecturally), not sufficient.

WHAT IS PRINTED
---------------
For every graph not already stored, a progress block: "name m=EDGES k=K CELLS cells [time]",
the full homology line, and the torsion line.  Then the summary table above, whose columns are

    graph            the name used as the data key
    edges            number of edges
    rank H_1         the free rank of H_1(B_k)
    torsion of H_{k-1}   the torsion in degree k - 1 (degree 2 for k = 3), as primary summands
    odd?             the sorted list of odd parts > 1 of the invariant factors, or "-"
    note             a short description of the graph

and a verdict.  If the first graph has odd torsion and the last does not, the script prints
"The odd torsion present at 20 edges is GONE at 36 edges ..." followed by the conclusion that
"has no odd torsion" is not minor-closed; otherwise it prints a warning that the expected
pattern did not appear.

USAGE
-----
    python3 experiments/density_interpolation.py [--k 3]

    --k K   (default 3) the number of particles.  The table reports the torsion of H_{K-1},
            which is where the inclusion-matrix torsion of K_{4,5} lives only for K = 3; other
            values are allowed but are not what the experiment is about.

Results are stored in data/density_interpolation.json, key "NAME|K", record
{"H": {degree: [free rank, [invariant factors]]}, "edges": m, "note": "..."}.  Stored graphs
are not recomputed, so with the shipped data a run is instant; from scratch at k = 3 it takes
about ten seconds and a few hundred MB (K_9 is the largest, 142,710 reduced cells).

Exit code: always 0; the verdict is in the output.

PITFALL
-------
Every graph in the table HAS torsion (2-torsion in H_1 and H_2), so the experiment says nothing
about torsion in general.  Its conclusion is about ODD torsion only.
"""

import argparse
import itertools
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import Graph, reduced_cell_count, paths              # noqa: E402
from gbt import fast, lean                                    # noqa: E402
from gbt.chain import odd_part                                # noqa: E402
from gbt.report import homology_line, torsion_line            # noqa: E402

LEAN_THRESHOLD = 600_000
DATA_FILE = "density_interpolation.json"

A = list(range(4))
B = list(range(4, 9))
BASE = [(a, b) for a in A for b in B]                  # K_{4,5}
AA = list(itertools.combinations(A, 2))                # edges inside the small side
BB = list(itertools.combinations(B, 2))                # edges inside the large side


def cases():
    """The six graphs of the interpolation, as (name, edge list on {0..8}, note) triples.

    BASE is K_{4,5}; AA and BB are the pairs inside A and inside B in lexicographic order, so
    "+ 3 in B" adds {4,5}, {4,6}, {4,7} and "+ 6 in B" adds the first six pairs of BB.
    """
    return [
        ("K_{4,5}", BASE, "the sparse case governed by the inclusion matrix"),
        ("K_{4,5} + 3 in B", BASE + BB[:3], "large side partly filled"),
        ("K_{4,5} + 6 in B", BASE + BB[:6], "large side half filled"),
        ("K_{4,5} + all A", BASE + AA, "small side filled"),
        ("K_{4,5} + all B", BASE + BB, "large side filled"),
        ("K_9", BASE + AA + BB, "both sides filled: the complete graph"),
    ]


def main():
    """Compute (or load) the six graphs at k particles, print the table and the verdict."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=3)
    args = ap.parse_args()
    k = args.k

    res = paths.load_json(DATA_FILE)
    rows = []
    for name, E, note in cases():
        key = "%s|%d" % (name, k)
        G = Graph(9, E)
        if key not in res:
            sz = reduced_cell_count(G, k)
            mod = lean if sz > LEAN_THRESHOLD else fast
            t0 = time.time()
            H, raw, _ = mod.compute(G, k, verbose=(sz > LEAN_THRESHOLD))
            print("%-20s m=%-3d k=%d  %9s cells  [%.0fs]"
                  % (name, len(E), k, format(sum(raw.values()), ","), time.time() - t0))
            print("%-20s %s" % ("", homology_line(H)))
            print("%-20s torsion: %s" % ("", torsion_line(H)), flush=True)
            res[key] = {"H": paths.record(H), "edges": len(E), "note": note}
            paths.save_json(DATA_FILE, res)
        rec = res[key]
        rows.append((name, rec["edges"],
                     {int(d): tuple(v) for d, v in rec["H"].items()}, rec["note"]))

    print("\n" + "=" * 92)
    print("from K_{4,5} to K_9 on the same 9 vertices, k = %d" % k)
    print("=" * 92)
    print("%-20s %6s %10s %-30s %-8s %s"
          % ("graph", "edges", "rank H_1", "torsion of H_%d" % (k - 1), "odd?", "note"))
    oddness = []
    for name, m, H, note in rows:
        t = H.get(k - 1, (0, []))[1]
        # odd part of each invariant factor: e.g. 6 -> 3, so a Z/6 is reported as odd torsion
        odd = sorted({odd_part(x) for x in t if odd_part(x) > 1})
        oddness.append((name, m, bool(odd)))
        print("%-20s %6d %10d %-30s %-8s %s"
              % (name, m, H[1][0],
                 torsion_line({k - 1: H.get(k - 1, (0, []))}), odd or "-", note))

    print()
    first = oddness[0]
    last = oddness[-1]
    if first[2] and not last[2]:
        print("The odd torsion present at %d edges is GONE at %d edges, on the same"
              % (first[1], last[1]))
        print("vertex set, with the sparse graph a subgraph of the dense one.")
        print("So 'has odd torsion' is not closed under adding edges, and 'has no odd")
        print("torsion' is not minor-closed: odd torsion in H_2 cannot be characterised")
        print("as 'contains one of finitely many minors'.")
    else:
        print("The expected pattern (odd torsion at the sparse end only) did not appear;")
        print("check that k = 3 and that the K_{4,5} row really shows a Z/3.")
    print("\nwrote data/%s" % DATA_FILE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
