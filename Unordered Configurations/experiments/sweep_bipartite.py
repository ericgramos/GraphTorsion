#!/usr/bin/env python3
"""
H_*(B_k K_{a,b}) for every (a, b, k) that fits: the complete bipartite sweep.

THE QUESTION
------------
This is where the exotic torsion lives.  For a = n+1 the top differential contains the
inclusion matrix W_{n-1,n}(b) (beside an identity block), whose cokernel supplies

        Z/3   in  H_2(B_3 K_{4,5}),           the first odd torsion found
        Z/4   in  H_3(B_4 K_{5,7}),           the first non-squarefree torsion found

and more generally Z/n for K_{n+1,2n-1}.  The sweep confirms the small cases and maps out where
the torsion stops being generic: in most cases the (k-1)-st homology has 2-torsion of
multiplicity exactly rank H_1, which is the "background" count, and the interesting graphs are
the ones that deviate from it -- K_{4,5} at k = 3 has (Z/2)^4 where the background would give
(Z/2)^12, plus a Z/3 on top.

METHOD
------
All triples 2 <= a <= AMAX, a <= b <= BMAX, 2 <= k <= KMAX are listed, their fully reduced
complexes are sized with `gbt.reduced_cell_count` (a closed formula, so nothing is built), and
those with at most CAP cells are run in order of size, cheapest first.  Each job computes the
full homology H_*(B_k K_{a,b}; Z) with `gbt.fast`, or `gbt.lean` above LEAN_THRESHOLD = 600,000
cells, and is saved immediately, so an interrupted run loses at most one job and a restart
skips everything stored.  Only a >= 2 and b >= a are run: K_{a,b} = K_{b,a}.

WHAT IS PRINTED
---------------
1.  "<J> jobs with at most <CAP> cells".
2.  For every job NOT already stored, a block
        K_a,b    k=<k>   <raw cells> cells -> <cells after Morse>  [time, fast|lean]
                   H_1=...  H_2=...
                   torsion: H_1: Z/2   H_2: ...
3.  A SUMMARY over ALL stored cases (not only this run's), in three tables:
      * "invariant factors seen anywhere" (at the stored data: [2, 6]);
      * "cases with ODD torsion": one line per case, K_a,b | k | torsion in every degree;
      * "the a = n+1 family": columns K_{a,b} | k | torsion of H_{k-1} | b = 2k-1?
        (only the cases with a = k+1, the family governed by coker W);
      * "the background pattern": columns
            K_{a,b} | k | rank H_1 | m(Z/2) | odd | generic?
        where m(Z/2) is the number of EVEN invariant factors of tors H_{k-1} (which equals the
        number of 2-primary summands, since every factor seen is 2 or 6), odd lists the odd
        parts > 1 of the factors, and generic? is "yes" when m(Z/2) = rank H_1 > 0.

Excerpt (`--cap 20000`, everything stored, a few seconds):

    cases with ODD torsion -- these are the ones coker W governs:
      K_4,5        k=3   H_1: Z/2   H_2: (Z/2)^4 + Z/3
      ...
      K_5,6        k=4   H_1: Z/2   H_2: (Z/2)^21   H_3: (Z/2)^5 + Z/3
    ...
      K_{a,b}           k   rank H_1     m(Z/2)      odd generic?
      K_{5,5}           3         16         16        - yes

USAGE
-----
    python3 experiments/sweep_bipartite.py [--cap 400000] [--amax 7] [--bmax 8] [--kmax 7]

    --cap N     skip jobs whose fully reduced complex has more than N cells (default 400,000:
                a few minutes; 700,000, the cap used for the stored data, takes an hour or so).
    --amax A    largest small side a (default 7).
    --bmax B    largest large side b (default 8).
    --kmax K    largest particle number k (default 7).

DATA FILE
---------
data/sweep_bipartite.json, keyed "K<a>,<b>|<k>" (e.g. "K4,5|3"), each record
{"H": {degree: [free rank, [invariant factors > 1]]}, "raw": total cells of the fully reduced
complex, "red": total cells after the Morse reduction}.  Records written by this version also
carry "cells_formula" (the predicted size); two older records hold only "H".

COST AND EXIT CODE
------------------
See --cap.  The largest jobs at the default cap need about a GB.  Always exits 0: the sweep
reports, it does not test.  `wilson_cokernel.py` and `top_degree_vanishing.py` test the
predictions against these numbers.

PITFALL
-------
To recompute a case delete its key from the data file; a stored case is never redone, even
when the cap changes.
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

#: above this many reduced cells, use the memory-lean pipeline instead of gbt.fast
LEAN_THRESHOLD = 600_000
DATA_FILE = "sweep_bipartite.json"


def main():
    """Run every job under the cap that is not yet stored, then print the summary.  Returns 0."""
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
    jobs = sorted(j for j in jobs if j[0] <= args.cap)     # cheapest first
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

    # every stored case, with degrees back to int and groups to (free, factors) tuples
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
        # number of even invariant factors = number of 2-primary cyclic summands
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
