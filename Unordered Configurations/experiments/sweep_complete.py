#!/usr/bin/env python3
"""
H_*(B_k K_n) for every (n, k) that fits: the complete-graph sweep.

THE QUESTION
------------
The complete graphs are the densest graphs on a given vertex set, and the striking finding is
how TAME their configuration spaces are.  Across every case computed -- K_3 to K_12, with k up to
12 for the small ones and k = 2, 3 for K_9 to K_12 -- the torsion is 2-torsion and nothing but
2-torsion, every summand is Z/2 with no Z/4, and for n >= 8 at k = 3 its multiplicity in H_2 is
exactly

        rank H_1(B_k K_n)  =  C(n-1, 2)

(rank H_1 = C(n-1, 2) holds for n >= 5; K_4 has rank 4), so the torsion looks like a shadow of
the first homology rather than something new.  The conjecture is that complete graphs carry
only Z/2 for all n and k.  Contrast `sweep_bipartite.py`: complete BIPARTITE graphs with a
sparse large side carry Z/3 and Z/4.  Density appears to destroy odd torsion -- see
`density_interpolation.py` for the experiment that isolates this.

METHOD
------
All pairs 3 <= n <= NMAX, 2 <= k <= KMAX are sized with `gbt.reduced_cell_count` (a closed
formula), and those with at most CAP cells are run cheapest first.  Each job computes the full
homology with `gbt.fast`, or `gbt.lean` above LEAN_THRESHOLD = 600,000 cells, and is saved at
once, so an interrupted run loses at most one job and a restart skips everything stored.  K_8 at
k = 4 and K_12 at k = 3 are beyond the default cap; they are in `large_cases.py`.

WHAT IS PRINTED
---------------
1.  "<J> jobs with at most <CAP> cells".
2.  For each job not yet stored: "K_n  k=..  <raw> cells -> <reduced>  [time, fast|lean]",
    then the homology and its torsion.
3.  A SUMMARY over all stored cases:
      * the cases with odd torsion anywhere (expected: none), and the set of invariant factors
        seen (expected: [2]);
      * "the multiplicity pattern at k = 3": columns graph | rank H_1 | m | C(n-1,2) | match,
        where m is the number of even invariant factors of tors H_2(B_3 K_n);
      * "full table": rows K_n, columns k, each entry the torsion in degrees >= 2 ("free" if
        none, "." if not computed).

Excerpt (`--cap 20000`, everything stored, a second):

    cases with ODD torsion anywhere: none
    invariant factors seen anywhere: [2]
      graph      rank H_1          m   C(n-1,2)
      K_7               15          1         15
      K_8               21         21         21 match
      ...
      K_12              55         55         55 match

USAGE
-----
    python3 experiments/sweep_complete.py [--cap 400000] [--nmax 13] [--kmax 12]

    --cap N     skip jobs whose fully reduced complex has more than N cells.  400,000 (the
                default) takes a few minutes; 1,600,000 a few hours and a couple of GB.
    --nmax N    largest n (default 13; K_13 at k = 3 has 1.6 million cells and needs a cap
                above that, and more memory than 7 GB through the full pipeline).
    --kmax K    largest k (default 12).

DATA FILE
---------
data/sweep_complete.json, keyed "K<n>|<k>" (e.g. "K8|3"), each record
{"H": {degree: [free rank, [invariant factors > 1]]}, "raw": cells of the fully reduced
complex, "red": cells after the Morse reduction}, plus "cells_formula" in records written by
this version.

COST AND EXIT CODE
------------------
See --cap.  Always exits 0: the script reports what it finds.  `validate.py` and
`crosscheck_mod_p.py` test these numbers independently.
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

#: above this many reduced cells, use the memory-lean pipeline instead of gbt.fast
LEAN_THRESHOLD = 600_000
DATA_FILE = "sweep_complete.json"


def main():
    """Run every job under the cap that is not yet stored, then print the summary.  Returns 0."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap", type=int, default=400_000)
    ap.add_argument("--nmax", type=int, default=13)
    ap.add_argument("--kmax", type=int, default=12)
    args = ap.parse_args()

    jobs = [(reduced_cell_count(complete(n), k), n, k)
            for n in range(3, args.nmax + 1)
            for k in range(2, args.kmax + 1)]
    jobs = sorted(j for j in jobs if j[0] <= args.cap)     # cheapest first
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
        # number of even invariant factors of tors H_2 = number of Z/2 summands here,
        # since every factor seen is exactly 2
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
