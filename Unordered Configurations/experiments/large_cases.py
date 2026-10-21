#!/usr/bin/env python3
"""
The individually expensive cases, run one at a time with the memory-lean pipeline.

The sweeps stop at a cell cap.  These are the cases beyond it that are worth the cost
because the theory makes a specific prediction, or because they decide a pattern:

    K55    K_{5,5} at k=4   the a = n+1, b <= 2n-1 case: the top group H_4 must VANISH
                            (top-degree vanishing), and coker W predicts Z/2 in H_3.
                            COMPLETED: H_4 = 0, H_3 = Z^1260 + Z/2.
    K56    K_{5,6} at k=4   the intermediate case: coker W predicts (Z/2)^5 + Z/3
                            in H_3.  COMPLETED: invariant factors 2,2,2,2,6.
    K8-4   K_8 at k=4       does the 2-torsion multiplicity keep tracking rank H_1?
                            COMPLETED: (Z/2)^22 in H_2.
    K12-3  K_12 at k=3      the same question at the far end of the k=3 family.
                            COMPLETED: (Z/2)^55 = (Z/2)^{rank H_1}.
    K57    K_{5,7} at k=4   where the flip-chain class predicts Z/4 in H_3 -- the first
                            non-squarefree torsion.  2.6M cells.  NOT COMPLETED HERE: the
                            full complex exhausted 7 GB during the peel.  The torsion of
                            H_3 alone has since been computed from two degree slices by
                            prime_census.py, and agrees with coker W_{3,4}(7).
    K13-3  K_13 at k=3      1.6M cells.  NOT COMPLETED, same reason.

Each case writes its own entry into data/large_cases.json, so they can be run separately,
in any order, on different machines, and the results pooled (by merging the JSON entries by
hand).  The two unfinished cases are included because they are the obvious targets for a
faster implementation of the peel -- see the note at the foot of validate_slow.py.

THE METHOD
----------
The full homology is computed by `gbt.lean.compute`: the fully reduced Swiatkowski complex is
built one degree at a time from the top down, unit pivots are cancelled (discrete Morse
reduction), and the small residual complex goes to the exact Smith normal form.  This is the
same answer as `gbt.fast.compute` in roughly a third of the memory.

WHAT IS PRINTED
---------------
With --list, or with no case names, a table with columns

    name             the case name to pass on the command line
    case             the graph and k
    cells            number of cells of the fully reduced complex (`reduced_cell_count`)
    what it decides  the prediction or question the case settles
    fits 7GB?        whether the full computation finished on a 7 GB machine
    stored?          "yes" if data/large_cases.json already has it

Real output (abridged):

    name     case                    cells  what it decides                        fits 7GB?
    K55      K_{5,5}, k=4          617,235  vanishing: H_4 = 0 (b = 5 <= 2n-1 = 7); ...  yes
    K57      K_{5,7}, k=4        2,615,835  flip-chain class: Z/4 in H_3           NO

When running a case: a header "=== CASE (N cells) -- what it decides ===", the progress lines
of `lean.compute` (verbose), then the full homology line, the torsion line with the run time,
and "wrote data/large_cases.json".  A stored case is not recomputed; its stored torsion
summary is printed instead ("... already stored: ...").

USAGE
-----
    python3 experiments/large_cases.py --list
    python3 experiments/large_cases.py K57
    python3 experiments/large_cases.py K55 K56
    python3 experiments/large_cases.py K8-4 --force

    cases     one or more of K55, K56, K8-4, K12-3, K57, K13-3, run in the order given
    --list    print the table above and exit (also the behaviour when no case is named)
    --force   recompute a case even if it is stored, overwriting the stored entry

Data: data/large_cases.json, key = the case name, record {"case": description, "k": k,
"H": {degree: [free rank, [invariant factors]]}, "raw": total cells of the reduced complex
before the Morse reduction, "red": total cells after it, "seconds": run time, "summary":
the torsion line}.  The stored K55 and K56 entries were written by an earlier version of the
script and have "raw" and "red" set to null and no "seconds".

Exit code: 1 if a case name is unknown, 0 otherwise (also after a successful run).  A run
killed by the operating system for lack of memory ends with exit code 137 from the shell and
writes nothing for that case.

WARNING on resources.  These are hours of CPU, and the two largest want several GB of RAM.
The Python unit-pivot peel is what breaks first: above roughly two million cells it becomes
the bottleneck, and above four million it exhausts 7 GB.  Run them detached, e.g.

    nohup python3 -u experiments/large_cases.py K57 > k57.log 2>&1 &
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import complete, complete_bip, reduced_cell_count    # noqa: E402
from gbt import lean, paths                                   # noqa: E402
from gbt.report import homology_line, torsion_line            # noqa: E402

DATA_FILE = "large_cases.json"

# name -> (description, graph builder, k, what it decides, known to finish in 7 GB?)
CASES = {
    "K55":   ("K_{5,5}, k=4", lambda: complete_bip(5, 5), 4,
              "vanishing: H_4 = 0 (b = 5 <= 2n-1 = 7); coker W: Z/2 in H_3", True),
    "K56":   ("K_{5,6}, k=4", lambda: complete_bip(5, 6), 4,
              "coker W: (Z/2)^5 + Z/3 in H_3", True),
    "K8-4":  ("K_8, k=4",     lambda: complete(8), 4,
              "2-torsion multiplicity vs rank H_1", True),
    "K12-3": ("K_12, k=3",    lambda: complete(12), 3,
              "multiplicity vs rank H_1 = C(11,2) = 55", True),
    "K57":   ("K_{5,7}, k=4", lambda: complete_bip(5, 7), 4,
              "flip-chain class: Z/4 in H_3", False),
    "K13-3": ("K_13, k=3",    lambda: complete(13), 3,
              "the largest complete graph attempted", False),
}


def main():
    """List the cases, or run the named ones in order and store each result as it finishes."""
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="*", help="case names; see --list")
    ap.add_argument("--list", action="store_true", help="show the cases and their sizes")
    ap.add_argument("--force", action="store_true", help="recompute a stored case")
    args = ap.parse_args()

    if args.list or not args.cases:
        stored = paths.load_json(DATA_FILE)
        print("%-8s %-16s %12s  %-38s %-10s %s"
              % ("name", "case", "cells", "what it decides", "fits 7GB?", "stored?"))
        for name, (desc, build, k, why, feasible) in CASES.items():
            print("%-8s %-16s %12s  %-38s %-10s %s"
                  % (name, desc, format(reduced_cell_count(build(), k), ","), why,
                     "yes" if feasible else "NO", "yes" if name in stored else ""))
        print("\nrun one with:  python3 experiments/large_cases.py <name>")
        return 0

    res = paths.load_json(DATA_FILE)
    for name in args.cases:
        if name not in CASES:
            print("unknown case %r; see --list" % name)
            return 1
        desc, build, k, why, feasible = CASES[name]
        if name in res and not args.force:
            print("%s already stored: %s" % (desc, res[name]["summary"]))
            continue
        G = build()
        print("=== %s  (%s cells)  --  %s ==="
              % (desc, format(reduced_cell_count(G, k), ","), why), flush=True)
        if not feasible:
            print("    NOTE: this case exhausted 7 GB in our runs.  Proceeding anyway;")
            print("    expect it to need more memory than that, or to be killed.",
                  flush=True)
        t0 = time.time()
        H, raw, red = lean.compute(G, k, verbose=True)
        dt = time.time() - t0
        print("   %s" % homology_line(H))
        print("   torsion: %s   [%.0fs]" % (torsion_line(H), dt), flush=True)
        res[name] = {"case": desc, "k": k, "H": paths.record(H),
                     "raw": sum(raw.values()), "red": sum(red.values()),
                     "seconds": round(dt), "summary": torsion_line(H)}
        paths.save_json(DATA_FILE, res)
        print("   wrote data/%s" % DATA_FILE, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
