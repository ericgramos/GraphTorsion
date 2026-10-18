#!/usr/bin/env python3
"""
Is there a Ko--Park theorem for intrinsic linking?  (No.)

Ko and Park (Discrete Comput. Geom. 48 (2012) 915-963) proved a clean statement in degree one:
H_1(B_k Gamma) has torsion if and only if Gamma is NOT planar, and then the torsion is a
single Z/2.  So in degree one the obstruction is a minor-closed property with a finite
forbidden list -- Kuratowski's K_5 and K_{3,3}.

The obvious guess one degree up is that H_2 torsion should detect the next topological
property in the sequence: intrinsic linking, whose forbidden minors are the seven graphs of
the PETERSEN FAMILY (Robertson--Seymour--Thomas).  This script tests that guess by computing
H_* for the whole family, at k = 3 and k = 4, and comparing with the linklessly embeddable
controls in `linkless_controls.py`.

What the data show (data/petersen_family.json):

    * at k = 3 all seven members have torsion-free H_2;
    * at k = 4 all seven members have exactly Z/2 in H_2;
    * but linklessly embeddable graphs -- the planar octahedron at k = 4, and the
      apex graph apex(planar8) already at k = 3 -- also carry H_2 torsion.

So "has H_2 torsion" does not characterise intrinsic linking at any fixed k: at k = 3 it
misses the whole family, and at k = 4 it also holds for linklessly embeddable graphs.

So there is no Ko--Park analogue in degree two.  The deeper reason: having
H_2 torsion is not closed under adding edges (see `density_interpolation.py`), so its
complement is not minor-closed and Robertson--Seymour does not apply at all.

The family is not hard-coded: it is generated as the closure of K_6 under Delta-Y and Y-Delta
exchanges, which is the construction that characterises it.

    python3 experiments/petersen_family.py [--k 3,4]

Requires networkx (for the isomorphism tests in the closure).  k=3 takes minutes; k=4 takes
a few hours, mostly on the Petersen graph itself.
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import petersen_family, reduced_cell_count, paths   # noqa: E402
from gbt import fast, lean                                  # noqa: E402
from gbt.report import has_torsion, homology_line, torsion_line   # noqa: E402

LEAN_THRESHOLD = 600_000
DATA_FILE = "petersen_family.json"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", default="3", help="comma-separated particle numbers")
    args = ap.parse_args()
    ks = [int(x) for x in args.k.split(",")]

    fam = petersen_family()
    print("the Petersen family: %d graphs, generated as the Delta-Y / Y-Delta closure of K_6\n"
          % len(fam))
    for name, G in fam:
        print("  %-12s %2d vertices %2d edges   degrees %s"
              % (name, G.N, G.M,
                 sorted((G.deg(v) for v in range(G.N)), reverse=True)))
    print()

    res = paths.load_json(DATA_FILE)
    for name, G in fam:
        for k in ks:
            key = "%s|%d" % (name, k)
            if key in res:
                continue
            sz = reduced_cell_count(G, k)
            mod = lean if sz > LEAN_THRESHOLD else fast
            t0 = time.time()
            H, raw, red = mod.compute(G, k, verbose=(sz > LEAN_THRESHOLD))
            print("%-12s n=%-2d k=%-2d  %9s cells -> %6s  [%.0fs]"
                  % (name, G.N, k, format(sum(raw.values()), ","),
                     format(sum(red.values()), ","), time.time() - t0))
            print("             %s" % homology_line(H))
            print("             torsion: %s" % torsion_line(H), flush=True)
            res[key] = {"H": paths.record(H), "vertices": G.N, "edges": G.M}
            paths.save_json(DATA_FILE, res)

    print("\n" + "=" * 72)
    print("VERDICT: does membership of the Petersen family predict H_2 torsion?")
    print("=" * 72)
    print("%-14s %4s %-30s %s" % ("graph", "k", "torsion", "torsion above H_1?"))
    any_free = []
    for name, G in fam:
        for k in ks:
            rec = res.get("%s|%d" % (name, k))
            if not rec:
                continue
            H = {int(d): tuple(v) for d, v in rec["H"].items()}
            above = has_torsion(H, min_degree=2)
            if not above:
                any_free.append((name, k))
            print("%-14s %4d %-30s %s" % (name, k, torsion_line(H),
                                          "yes" if above else "NO"))
    print()
    if any_free:
        print("Members of the Petersen family with NO torsion above H_1: %s" % any_free)
        print("So at these k, intrinsic linking does not imply H_2 torsion (at k = 4 every member")
        print("does acquire Z/2: run with --k 4).")
    else:
        print("Every member has torsion above H_1 at the k tried; but so do linklessly")
        print("embeddable graphs, so this does not characterise intrinsic linking.")
    print("For the converse -- linklessly embeddable graphs that DO have H_2 torsion --")
    print("run linkless_controls.py.  The octahedron is planar and has Z/2 in H_2 at k=4.")
    print("\nwrote data/%s" % DATA_FILE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
