#!/usr/bin/env python3
"""
The expensive regression checks: the million-cell cases, one at a time.

WHY
---
These are headline numbers that cannot be recomputed in minutes, so they are kept out of
`validate.py`.  Each is a single full homology computation with `gbt.fast`; run them
individually, because together they are hours of CPU and the larger ones want several GB.

THE CASES
---------
    #  case                                  cells (reduced)  expected
    0  K_{5,5}, k = 4                        617,235          H_4 = 0;  H_3 = Z^1260 + Z/2
    1  K_{5,6}, k = 4                        1,355,675        H_4 = 0;  H_3 = Z^3010 + (Z/2)^4 + Z/6
    2  K_8, k = 4                            1,065,001        H_2 = Z^1434 + (Z/2)^22
    3  K_12, k = 3                           971,036          H_2 = Z^6071 + (Z/2)^55

Cases 0 and 1 test the top-degree vanishing (b <= 2n-1) and the torsion one degree down, which
is the torsion of coker W_{3,4}(b); cases 2 and 3 test the complete-graph pattern.

WHAT IS PRINTED
---------------
With --list, or with no arguments: the table above (# | case | cells | expected).  Otherwise,
for each selected case: "=== <case> ===", the progress messages of gbt.fast, one line
"[ok  ] H_d = <group>" (or "[FAIL] ... expected ...") per checked degree, and the time; then
"ALL SELECTED CHECKS PASSED" or "<m> FAILURE(S): [...]".  Case 0 took about two and a half
minutes on the machine used here.

USAGE
-----
    python3 experiments/validate_slow.py --list
    python3 experiments/validate_slow.py 0
    python3 experiments/validate_slow.py 0 1 2

    which      case numbers to run (0-3).
    --list     print the table and exit (also the behaviour with no case numbers).

EXIT CODE AND OUTPUT FILES
--------------------------
Exits with 1 if a selected check fails, else 0.  Writes nothing to data/ (the same numbers are
stored by `large_cases.py`).

ON THE RESOURCE CEILING
-----------------------
The Python unit-pivot peel is the bottleneck.  Above roughly two million cells it dominates the
runtime, and above four million it exhausts 7 GB -- which is why K_{5,7} at k = 4 and K_13 at
k = 3 sit beyond what the full-homology pipelines reach (the torsion of a single degree goes
further: see `prime_census.py`).  A compiled implementation of the peel, or a sparse Smith form
over Z with better pivoting, is the obvious next step for anyone wanting to go further.
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import complete, complete_bip, fmt, reduced_cell_count   # noqa: E402
from gbt import fast, lean                                       # noqa: E402

# (description, builder, k, [(degree, expected (free rank, invariant factors))], pipeline)
CASES = [
    ("K_{5,5} k=4  [vanishing: H_4 = 0 since b = 5 <= 2n-1 = 7]",
     lambda: complete_bip(5, 5), 4,
     [(4, (0, ())), (3, (1260, (2,)))], "fast"),
    ("K_{5,6} k=4  [vanishing again, and the torsion one degree down]",
     lambda: complete_bip(5, 6), 4,
     [(4, (0, ())), (3, (3010, (2, 2, 2, 2, 6)))], "fast"),
    ("K_8 k=4      [does the multiplicity keep tracking rank H_1?]",
     lambda: complete(8), 4,
     [(2, (1434, tuple([2] * 22)))], "fast"),
    ("K_12 k=3     [the far end of the k=3 family: C(11,2) = 55]",
     lambda: complete(12), 3,
     [(2, (6071, tuple([2] * 55)))], "fast"),
]


def main():
    """List the cases, or run the selected ones and compare with the expected groups.

    Returns 1 if a selected check fails, else 0."""
    ap = argparse.ArgumentParser()
    ap.add_argument("which", nargs="*", type=int, help="case indices; see --list")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    if args.list or not args.which:
        print("%3s %-56s %12s %s" % ("#", "case", "cells", "expected"))
        for i, (desc, build, k, checks, _) in enumerate(CASES):
            print("%3d %-56s %12s %s"
                  % (i, desc, format(reduced_cell_count(build(), k), ","),
                     ";  ".join("H_%d = %s" % (d, fmt(w[0], list(w[1])))
                                for d, w in checks)))
        print("\nrun one with:  python3 experiments/validate_slow.py 0")
        return 0

    failures = []
    for i in args.which:
        desc, build, k, checks, pipe = CASES[i]
        print("=== %s ===" % desc, flush=True)
        t0 = time.time()
        mod = fast if pipe == "fast" else lean
        H, _, _ = mod.compute(build(), k, verbose=True)
        for d, want in checks:
            # a missing degree means the complex has no cells there, i.e. H_d = 0
            got = (H.get(d, (0, ()))[0], tuple(H.get(d, (0, ()))[1]))
            ok = (got == want)
            if not ok:
                failures.append("%s, H_%d" % (desc, d))
            print("  [%s] H_%d = %s%s"
                  % ("ok  " if ok else "FAIL", d, fmt(got[0], list(got[1])),
                     "" if ok else "   expected %s" % fmt(want[0], list(want[1]))),
                  flush=True)
        print("  [%.0fs]" % (time.time() - t0), flush=True)
    print()
    if failures:
        print("%d FAILURE(S): %s" % (len(failures), failures))
        return 1
    print("ALL SELECTED CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
