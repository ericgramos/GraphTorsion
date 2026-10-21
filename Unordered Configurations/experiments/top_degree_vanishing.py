#!/usr/bin/env python3
"""
Top-degree vanishing and the trichotomy conjecture: when does H_n(B_n K_{a,b}) vanish?

THE QUESTION
------------
The top interesting degree of B_n(Gamma) for a bipartite Gamma is n, and the expected answer is

    (1)  a <= n                      =>  H_n(B_n K_{a,b}) = 0        (too few A-vertices)
    (2)  a = n+1 and b <= 2n-1       =>  H_n(B_n K_{a,b}) = 0        (inclusion-matrix rank)

together with the exact dimension in case a = n+1 for every b,

         dim H_n(B_n K_{n+1,b})  =  max( 0, C(b,n) - C(b,n-1) ) ,

which holds because the top differential is the unsigned inclusion matrix W_{n-1,n}(b) beside
an identity block, so the top cycles are exactly ker W.  The vanishing in (2) is then the
numerical fact

         b <= 2n-1   <=>   C(b,n-1) >= C(b,n) ,

a GLOBAL count: no local propagation argument sees it, which is why the proof has to pass
through Gottlieb's rank theorem for inclusion matrices.

Conjecturally the list is complete -- H_n(B_n K_{a,b}) != 0 whenever a >= n+2, and whenever
a = n+1 and b >= 2n.  That converse is open, and this script is the evidence for it: every
triple (a,b,n) small enough to compute is checked against the prediction, and any disagreement
is reported as a counterexample.

METHOD
------
All triples 2 <= a <= AMAX, a <= b <= BMAX, 2 <= n <= NMAX whose fully reduced complex has at
most CAP cells are run, cheapest first.  For each, the full homology of B_n K_{a,b} is
computed (gbt.fast, or gbt.lean above 600,000 cells) and dim H_n -- the free rank in degree n;
H_n is the top degree, so it has no torsion -- is compared with predicted_dim(a, b, n):
an exact value for a <= n+1, and "nonzero" for a >= n+2.  Since only a <= b is run, the
smaller side is always a.

WHAT IS PRINTED
---------------
1.  "<J> triples (a,b,n) with at most <CAP> reduced cells", then one row per triple, columns
        K_{a,b} | n | cells | dim H_n | predicted ("-" when only the conjecture applies) |
        verdict | H_n
    where verdict is "ok" (exact prediction met), "MISMATCH", "conj !=0" (a >= n+2 and H_n is
    nonzero, as conjectured) or "CONJ FAILS".
2.  The number of cases, how many vanishing and nonzero exact predictions were confirmed, how
    many conjectural cases were checked, and either "no disagreements" or the list of them.
3.  "the a = n+1 family": columns K_{n+1,b} | n | C(b,n-1) | C(b,n) | dim H_n | b <= 2n-1?,
    for the cases of THIS run (sorted numerically by a, b, n).
4.  "the torsion it forces in degree n-1": for every a = n+1 case of this run below 600,000
    cells, the torsion of H_{n-1} (recomputed), which `wilson_cokernel.py` compares with
    coker W_{n-1,n}(b).

Excerpt (`--cap 5000`, under a second):

    K_{a,b}         n      cells    dim H_n  predicted verdict   H_n
    K_{2,2}         2         32          0          0 ok        0
    ...
    39 cases in 1 s
      vanishing cases confirmed       : 22
      nonzero cases with exact formula: 7
      cases with a >= n+2 (conjecture): 10
    no disagreements: every computable triple matches the prediction.

USAGE
-----
    python3 experiments/top_degree_vanishing.py [--cap 400000] [--amax 7] [--bmax 10]
                                                [--nmax 5]

    --cap N     skip triples whose fully reduced complex has more than N cells (default
                400,000: a few minutes; a few million reaches further, with gbt.lean).
    --amax A    largest a (default 7).     --bmax B    largest b (default 10).
    --nmax N    largest particle number n (default 5).

DATA FILE
---------
data/top_degree_vanishing.json, keyed "K<a>,<b>|<n>" (e.g. "K4,5|3"), each record
{"cells": reduced cells, "H_n": [free rank, [invariant factors]] (H_n only),
 "predicted": the predicted dimension, or null for a >= n+2}.  New results are MERGED into the
stored file, so a run with a smaller cap never discards stored cases.

COST AND EXIT CODE
------------------
See --cap; memory as for the bipartite sweep.  Exits with 1 if any triple disagrees with the
prediction (an exact mismatch, or a vanishing H_n for a >= n+2), otherwise 0.

PITFALL
-------
Nothing is cached on input: every run recomputes every triple under the cap (the stored file
is only written to), and step 4 recomputes the a = n+1 cases once more.
"""

import argparse
import os
import sys
import time
from math import comb

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import complete_bip, fmt, primary, reduced_cell_count   # noqa: E402
from gbt import fast, lean, paths                               # noqa: E402

LEAN_THRESHOLD = 600_000       # above this many cells, use the memory-lean pipeline


def predicted_dim(a, b, n):
    """The predicted dim H_n(B_n K_{a,b}) (an int), or None where only the conjecture applies.

    a <= n gives 0; a = n+1 gives max(0, C(b,n) - C(b,n-1)), the rank of ker W_{n-1,n}(b);
    a >= n+2 gives None (conjecturally nonzero)."""
    if a <= n:
        return 0                                   # case (1)
    if a == n + 1:
        return max(0, comb(b, n) - comb(b, n - 1))  # exact: ker of the inclusion matrix
    return None                                   # a >= n+2: conjecturally nonzero


def main():
    """Check every triple under the cap, print the tables, merge the results into the data
    file.  Returns 1 if some triple disagrees with the prediction, else 0."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap", type=int, default=400_000,
                    help="skip any case with more than this many reduced cells")
    ap.add_argument("--amax", type=int, default=7)
    ap.add_argument("--bmax", type=int, default=10)
    ap.add_argument("--nmax", type=int, default=5)
    args = ap.parse_args()

    jobs = []
    for a in range(2, args.amax + 1):
        for b in range(a, args.bmax + 1):
            for n in range(2, args.nmax + 1):
                G = complete_bip(a, b)
                sz = reduced_cell_count(G, n)
                if sz <= args.cap:
                    jobs.append((sz, a, b, n))
    jobs.sort()
    print("%d triples (a,b,n) with at most %s reduced cells\n"
          % (len(jobs), format(args.cap, ",")), flush=True)

    hdr = ("%-12s %4s %10s %10s %10s %-9s %s"
           % ("K_{a,b}", "n", "cells", "dim H_n", "predicted", "verdict", "H_n"))
    print(hdr)
    print("-" * len(hdr))

    results = {}
    bad = []
    counts = {"vanishes": 0, "nonzero": 0, "exact dim matched": 0, "no prediction": 0}
    t_start = time.time()
    for sz, a, b, n in jobs:
        G = complete_bip(a, b)
        mod = lean if sz > LEAN_THRESHOLD else fast
        H, _, _ = mod.compute(G, n)
        # H has no degree-n entry when n exceeds the number of vertices: then there are
        # no n-cells at all and H_n = 0.
        free, tors = H.get(n, (0, []))
        pred = predicted_dim(a, b, n)

        if pred is None:
            verdict = "conj !=0" if free or tors else "CONJ FAILS"
            counts["no prediction"] += 1
            if not (free or tors):
                bad.append((a, b, n, "conjecture predicts nonzero, got 0"))
        elif (free, bool(tors)) == (pred, False):
            verdict = "ok"
            counts["exact dim matched"] += 1
            counts["vanishes" if pred == 0 else "nonzero"] += 1
        else:
            verdict = "MISMATCH"
            bad.append((a, b, n, "got dim %d%s, expected %d"
                        % (free, " + torsion" if tors else "", pred)))

        print("%-12s %4d %10s %10d %10s %-9s %s"
              % ("K_{%d,%d}" % (a, b), n, format(sz, ","), free,
                 "-" if pred is None else pred, verdict, fmt(free, tors)), flush=True)
        results["K%d,%d|%d" % (a, b, n)] = {
            "cells": sz, "H_n": [free, list(tors)], "predicted": pred,
        }

    print("-" * len(hdr))
    print("\n%d cases in %.0f s" % (len(jobs), time.time() - t_start))
    print("  vanishing cases confirmed       : %d" % counts["vanishes"])
    print("  nonzero cases with exact formula: %d" % counts["nonzero"])
    print("  cases with a >= n+2 (conjecture): %d" % counts["no prediction"])
    if bad:
        print("\n%d DISAGREEMENT(S):" % len(bad))
        for a, b, n, why in bad:
            print("   K_{%d,%d}, n=%d: %s" % (a, b, n, why))
    else:
        print("\nno disagreements: every computable triple matches the prediction.")

    # The a = n+1 row is where the vanishing theorem lives; print it on its own, with the
    # numerical equivalence b <= 2n-1  <=>  C(b,n-1) >= C(b,n) made explicit.
    print("\nthe a = n+1 family (vanishing, and torsion one degree down):")
    print("  %-12s %3s %9s %9s %12s %s"
          % ("K_{n+1,b}", "n", "C(b,n-1)", "C(b,n)", "dim H_n", "b <= 2n-1?"))
    for key, rec in sorted(results.items(),
                           key=lambda kv: [int(x) for x in kv[0][1:].replace('|', ',').split(',')]):
        ab, n = key.split("|")
        a, b = (int(x) for x in ab[1:].split(","))
        n = int(n)
        if a != n + 1:
            continue
        print("  %-12s %3d %9d %9d %12d %s"
              % ("K_{%d,%d}" % (a, b), n, comb(b, n - 1), comb(b, n),
                 rec["H_n"][0], "yes" if b <= 2 * n - 1 else "no"))

    # Torsion one degree down, where coker W_{n-1,n}(b) predicts it exactly.
    print("\nand the torsion it forces in degree n-1 (see wilson_cokernel.py):")
    for sz, a, b, n in jobs:
        if a != n + 1 or sz > LEAN_THRESHOLD:
            continue
        H, _, _ = fast.compute(complete_bip(a, b), n)
        if H.get(n - 1, (0, []))[1]:
            print("  K_{%d,%d}, n=%d:  H_%d torsion = %s"
                  % (a, b, n, n - 1, primary(H[n - 1][1])))

    # Merge rather than overwrite, so a run with a smaller --cap never discards stored cases.
    stored = paths.load_json("top_degree_vanishing.json")
    stored.update(results)
    paths.save_json("top_degree_vanishing.json", stored)
    print("\nwrote data/top_degree_vanishing.json (%d new or updated, %d stored in all)"
          % (len(results), len(stored)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
