#!/usr/bin/env python3
"""
H_*(B_k K_{a,b}) for every (a, b, k) that fits: the complete bipartite sweep.

This is where the exotic torsion lives.  For a = n+1 the top differential contains
the inclusion matrix W_{n-1,n}(b) (beside an identity block), whose cokernel supplies

        Z/3   in  H_2(B_3 K_{4,5}),           the first odd torsion in a graph braid group
        Z/4   in  H_3(B_4 K_{5,7}),           the first non-squarefree torsion

and more generally Z/n for K_{n+1,2n-1}.  The sweep confirms that and maps out where the
torsion stops being generic: in most cases the (k-1)-st homology has 2-torsion of
multiplicity exactly rank H_1, which is the "background" count, and the
interesting graphs are the ones that deviate from it -- K_{4,5} has (Z/2)^4 where the
background would give (Z/2)^12, plus a Z/3 on top.

Jobs are ordered by `gbt.reduced_cell_count` and results are written to
data/sweep_bipartite.json after each one, so the run is resumable.

    python3 experiments/sweep_bipartite.py [--cap CELLS] [--amax A] [--bmax B] [--kmax K]

    --cap   400000   a few minutes
    --cap   700000   an hour or so (the cap used for the stored data)
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import complete_bip, reduced_cell_count         # noqa: E402
from gbt import fast, lean, paths                        # noqa: E402
from gbt.chain import odd_part                           # noqa: E402
from gbt.report import (exponents, has_odd_torsion,       # noqa: E402
                        homology_line, torsion_line, torsion_summary)

LEAN_THRESHOLD = 600_000
DATA_FILE = "sweep_bipartite.json"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap", type=int, default=400_000)
    ap.add_argument("--amax", type=int, default=7)
    ap.add_argument("--bmax", type=int, default=8)
    ap.add_argument("--kmax", type=int, default=7)
    args = ap.parse_args()

    jobs = [(reduced_cell_count(complete_bip(a, b), k), a, b, k)
            for a in range(2, args.amax + 1)
            for b in range(a, args.bmax + 1)
            for k in range(2, args.kmax + 1)]
    jobs = sorted(j for j in jobs if j[0] <= args.cap)
    print("%d jobs with at most %s cells\n" % (len(jobs), format(args.cap, ",")),
          flush=True)

    res = paths.load_json(DATA_FILE)
    for sz, a, b, k in jobs:
        key = "K%d,%d|%d" % (a, b, k)
        if key in res:
            continue
        mod = lean if sz > LEAN_THRESHOLD else fast
        t0 = time.time()
        H, raw, red = mod.compute(complete_bip(a, b), k, verbose=(sz > LEAN_THRESHOLD))
        print("K_%d,%-3d k=%-2d %9s cells -> %6s  [%.0fs, %s]"
              % (a, b, k, format(sum(raw.values()), ","), format(sum(red.values()), ","),
                 time.time() - t0, mod.__name__.split(".")[-1]))
        print("           %s" % homology_line(H))
        print("           torsion: %s" % torsion_line(H), flush=True)
        res[key] = {"H": paths.record(H), "raw": sum(raw.values()),
                    "red": sum(red.values()), "cells_formula": sz}
        paths.save_json(DATA_FILE, res)

    # ---------------------------------------------------------------- summary
    print("\n" + "=" * 78)
    print("SUMMARY over %d stored cases" % len(res))
    print("=" * 78)

    allH = {key: {int(d): tuple(v) for d, v in rec["H"].items()}
            for key, rec in res.items()}
    print("invariant factors seen anywhere: %s"
          % sorted({x for H in allH.values() for x in exponents(H)}))

    print("\ncases with ODD torsion -- these are the ones coker W governs:")
    any_odd = False
    for key in sorted(allH, key=lambda s: (len(s), s)):
        H = allH[key]
        if not has_odd_torsion(H):
            continue
        any_odd = True
        ab, k = key.split("|")
        print("  K_%-10s k=%-2s  %s" % (ab[1:], k, torsion_line(H)))
    if not any_odd:
        print("  none at this cap -- raise it to at least 20,000 to reach K_{4,5} at k=3")

    print("\nthe a = n+1 family, where the theory predicts Z/n at b = 2n-1:")
    print("  %-14s %4s %-28s %s" % ("K_{a,b}", "k", "torsion of H_{k-1}", "b = 2k-1?"))
    for key in sorted(allH, key=lambda s: (len(s), s)):
        ab, k = key.split("|")
        a, b = (int(x) for x in ab[1:].split(","))
        k = int(k)
        if a != k + 1:
            continue
        H = allH[key]
        t = torsion_summary(H).get(k - 1, "none")
        print("  %-14s %4d %-28s %s"
              % ("K_{%d,%d}" % (a, b), k, t, "yes" if b == 2 * k - 1 else ""))

    print("\nthe background pattern: 2-torsion multiplicity in H_{k-1} vs rank H_1")
    print("  %-14s %4s %10s %10s %8s %s"
          % ("K_{a,b}", "k", "rank H_1", "m(Z/2)", "odd", "generic?"))
    for key in sorted(allH, key=lambda s: (len(s), s)):
        ab, k = key.split("|")
        k = int(k)
        H = allH[key]
        if k - 1 not in H:
            continue
        m = sum(1 for x in H[k - 1][1] if x % 2 == 0)
        odd = ",".join(str(odd_part(x)) for x in H[k - 1][1] if odd_part(x) > 1) or "-"
        b1 = H[1][0]
        print("  %-14s %4d %10d %10d %8s %s"
              % ("K_{%s}" % ab[1:], k, b1, m, odd,
                 "yes" if (m == b1 and b1) else ""))

    print("\nwrote data/%s" % DATA_FILE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
