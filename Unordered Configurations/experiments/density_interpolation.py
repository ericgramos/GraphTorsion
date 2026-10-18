#!/usr/bin/env python3
"""
Does density destroy the exotic odd torsion?  Interpolating from K_{4,5} to K_9.

K_{4,5} at k=3 has H_2 = Z^155 + (Z/2)^4 + Z/3: the first odd torsion in a graph braid group.
K_9 at k=3 has H_2 = Z^1401 + (Z/2)^28: no odd torsion at all, and a much larger 2-torsion
multiplicity.  Both graphs have the same nine vertices, and the second contains the first.  So
somewhere along the way the Z/3 dies, and the 2-torsion multiplicity jumps from 4 to 28.

This script walks from one to the other by adding edges inside the parts:

    K_{4,5}                     20 edges    (Z/2)^4 + Z/3
    K_{4,5} + 3 edges in B      23 edges
    K_{4,5} + 6 edges in B      26 edges
    K_{4,5} + all edges in B    30 edges
    K_{4,5} + all edges in A    26 edges
    K_9                         36 edges    (Z/2)^28

Two conclusions, both of which matter for the conjectures.

  1.  The Z/3 disappears as soon as the sparse side is filled in.  This is the mechanism
      behind two conjectures: odd torsion needs a K_{4,5} minor, and complete graphs have
      only 2-primary torsion: the inclusion matrix W_{n-1,n}(b) only governs the top
      differential when the large side is INDEPENDENT.  Adding edges there adds cells that
      are not matching configurations, and the Wilson block stops being the whole story.

  2.  Having odd torsion is therefore NOT closed under adding edges.  Hence neither is
      having torsion, so the complement is not minor-closed and Robertson--Seymour cannot
      produce a finite forbidden list -- conjecturally there is no Kuratowski theorem above
      degree one.  Note how this cuts: K_9 CONTAINS K_{4,5} as a subgraph, and has no
      odd torsion.  Containment is necessary, not sufficient.

    python3 experiments/density_interpolation.py [--k 3]
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
    return [
        ("K_{4,5}", BASE, "the sparse case governed by the inclusion matrix"),
        ("K_{4,5} + 3 in B", BASE + BB[:3], "large side partly filled"),
        ("K_{4,5} + 6 in B", BASE + BB[:6], "large side half filled"),
        ("K_{4,5} + all A", BASE + AA, "small side filled"),
        ("K_{4,5} + all B", BASE + BB, "large side filled"),
        ("K_9", BASE + AA + BB, "both sides filled: the complete graph"),
    ]


def main():
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
        print("So 'has odd torsion' is not closed under adding edges, and")
        print("'is torsion-free' is not minor-closed: no finite forbidden-minor")
        print("characterisation of H_2 torsion can exist.")
    else:
        print("The expected pattern (odd torsion at the sparse end only) did not appear;")
        print("check that k = 3 and that the K_{4,5} row really shows a Z/3.")
    print("\nwrote data/%s" % DATA_FILE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
