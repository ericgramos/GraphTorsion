#!/usr/bin/env python3
"""
H_*(B_k K_n) for every (n, k) that fits: the complete-graph sweep.

The complete graphs are the densest graphs on a given vertex set, and the striking finding
is how TAME their configuration spaces are.  Across every case computed -- 3 <= n <= 12 and
2 <= k <= 12 -- the torsion is 2-torsion and nothing but 2-torsion, every summand is Z/2
with no Z/4, and for n >= 8 at k = 3 its multiplicity is exactly

        rank H_1(B_k K_n)  =  C(n-1, 2),

so the torsion looks like a shadow of the first homology rather than something new.  The conjecture
is that this holds for all n and k.  Contrast sweep_bipartite.py: complete BIPARTITE graphs with a
sparse large side carry Z/3 and Z/4.  Density appears to destroy odd torsion -- see
`density_interpolation.py` for the experiment that isolates this.

HOW THE SWEEP WORKS
-------------------
Jobs are ordered by `gbt.reduced_cell_count`, which is a closed form, so the cheap cases are
done first and an interrupted run loses at most one job: results are written to
data/sweep_complete.json after each one and reloaded on restart.  Pass a cell cap to control
how far it goes; the pipeline is chosen by size (`gbt.fast` below the threshold,
`gbt.lean` above it).

    python3 experiments/sweep_complete.py [--cap CELLS] [--nmax N] [--kmax K]

    --cap   400000   a few minutes
    --cap  1600000   a few hours, and needs a couple of GB
"""

import argparse
import os
import sys
import time
from math import comb

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import complete, reduced_cell_count            # noqa: E402
from gbt import fast, lean, paths                       # noqa: E402
from gbt.report import (exponents, has_odd_torsion,      # noqa: E402
                        homology_line, torsion_line, torsion_summary)

LEAN_THRESHOLD = 600_000
DATA_FILE = "sweep_complete.json"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap", type=int, default=400_000)
    ap.add_argument("--nmax", type=int, default=13)
    ap.add_argument("--kmax", type=int, default=12)
    args = ap.parse_args()

    jobs = [(reduced_cell_count(complete(n), k), n, k)
            for n in range(3, args.nmax + 1)
            for k in range(2, args.kmax + 1)]
    jobs = sorted(j for j in jobs if j[0] <= args.cap)
    print("%d jobs with at most %s cells\n" % (len(jobs), format(args.cap, ",")),
          flush=True)

    res = paths.load_json(DATA_FILE)
    for sz, n, k in jobs:
        key = "K%d|%d" % (n, k)
        if key in res:
            continue
        mod = lean if sz > LEAN_THRESHOLD else fast
        t0 = time.time()
        H, raw, red = mod.compute(complete(n), k, verbose=(sz > LEAN_THRESHOLD))
        print("K_%-3d k=%-2d  %9s cells -> %6s  [%.0fs, %s]"
              % (n, k, format(sum(raw.values()), ","), format(sum(red.values()), ","),
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
    odd = [key for key, H in allH.items() if has_odd_torsion(H)]
    print("cases with ODD torsion anywhere: %s" % (odd if odd else "none"))
    seen = sorted({x for H in allH.values() for x in exponents(H)})
    print("invariant factors seen anywhere: %s" % seen)
    print("  (a single value 2 means every torsion summand is Z/2: no Z/4, no odd torsion)")

    print("\nthe multiplicity pattern at k = 3:  (Z/2)^m in H_2  vs  rank H_1 = C(n-1,2)")
    print("  %-8s %10s %10s %10s %s" % ("graph", "rank H_1", "m", "C(n-1,2)", ""))
    for n in range(3, args.nmax + 1):
        rec = res.get("K%d|3" % n)
        if not rec:
            continue
        H = {int(d): tuple(v) for d, v in rec["H"].items()}
        m = sum(1 for x in H.get(2, (0, []))[1] if x % 2 == 0)
        print("  K_%-6d %10d %10d %10d %s"
              % (n, H[1][0], m, comb(n - 1, 2),
                 "match" if m == comb(n - 1, 2) else ""))

    print("\nfull table (torsion only, degrees >= 2):")
    ks = sorted({int(key.split("|")[1]) for key in res})
    ns = sorted({int(key.split("|")[0][1:]) for key in res})
    print("  %-7s%s" % ("", "".join("%15s" % ("k=%d" % k) for k in ks)))
    for n in ns:
        cells = []
        for k in ks:
            rec = res.get("K%d|%d" % (n, k))
            if not rec:
                cells.append(".")
                continue
            H = {int(d): tuple(v) for d, v in rec["H"].items()}
            t = {d: s for d, s in torsion_summary(H).items() if d >= 2}
            cells.append(";".join("H%d:%s" % (d, s) for d, s in t.items()) or "free")
        print("  K_%-5d%s" % (n, "".join("%15s" % c for c in cells)))

    print("\nwrote data/%s" % DATA_FILE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
