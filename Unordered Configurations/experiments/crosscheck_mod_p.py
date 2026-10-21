#!/usr/bin/env python3
"""
An independent check of every torsion claim: Betti numbers over F_p from the RAW complex.

THE QUESTION
------------
Every integral answer in this repository comes out of the same two steps: a unit-pivot
(discrete Morse) reduction, then an exact Smith normal form.  Both are places where a bug would
be invisible.  A reduction that forgets to delete a cancelled generator produces a complex whose
homology is simply wrong, and that did happen once during this project.  So it is worth having
a second route to the same information that shares none of that machinery.  This script is
that route.

THE METHOD
----------
For a prime p, dim H_d(X; F_p) is computed straight from the boundary matrices of the RAW
(unreduced) Swiatkowski complex by sparse elimination over F_p (`gbt.linalg.betti_mod_p`): no
Morse reduction, no Smith normal form, and over F_2 not even division, just xor on bitsets.
The universal coefficient theorem then turns the integral answer into a prediction: if

        H_d(X; Z)  =  Z^{f_d}  (+)  (+)_i Z/q_{d,i}

then

        dim H_d(X; F_p)  =  f_d  +  #{i : p | q_{d,i}}  +  #{i : p | q_{d-1,i}} .

The last term is the Tor contribution from one degree down.  So the integral answer PREDICTS
the mod-p Betti numbers, and this script tests the prediction (`predicted_betti`).  Comparing
p = 2 with a large prime (1000003, standing in for Q) also isolates the 2-torsion on its own:
the difference between the two Betti vectors in degree d counts the 2-primary summands of H_d
plus those of H_{d-1}.  The raw complex is used deliberately -- the one that has had nothing
done to it.

Two passes are made:

  (1) ten named cases (K_5, K_{3,3}, the triangular prism, the octahedron, K_6, K_7, K_{4,5} at
      various k), each checked over F_2 AND over F_1000003;
  (2) every record in data/sweep_complete.json and data/sweep_bipartite.json that HAS torsion
      somewhere (torsion-free records are skipped: there is nothing to cross-check), checked
      over F_2 only, since that is the prime that detects the 2-torsion.  Odd torsion is not
      tested in this pass; it is covered by pass (1) (K_{4,5}) and by `wilson_cokernel.py`.

WHAT IS PRINTED
---------------
Pass (1) is a table with columns

    case             the graph and the number of particles k
    raw cells        number of cells of the UNREDUCED complex (`gbt.cell_count`)
    F_2 Betti        [dim H_0, dim H_1, ...] over F_2, computed from the raw complex
    F_1000003 Betti  the same over F_1000003
    match            "ok" if BOTH vectors equal the universal-coefficient prediction from the
                     integral homology, otherwise "MISMATCH" (followed by two lines giving the
                     predicted and computed vectors); then the run time of that row.

Under a row, "2-torsion witness" lists {degree: (F_2 Betti) - (F_q Betti)} wherever the two
differ.  Cases whose raw complex exceeds --cap are listed as "(skipped: above the cap)".
Pass (2) prints one line per stored record checked, then "N verified, M skipped as too large".
The final verdict is either

    ALL CROSS-CHECKS PASSED: the mod-p Betti numbers of the raw complexes agree
    with the integral answers from the reduce-then-SNF pipeline in every case.

or "N MISMATCH(ES): [...]" naming the failing cases.  Excerpt (with --cap 3000):

    case              raw cells F_2 Betti                  F_1000003 Betti            match
    K_5, k=2                555 [1 7 1]                    [1 6 0]                    ok   [0s]
          2-torsion witness (F_2 minus F_q, by degree): {1: 1, 2: 1}
    K_{3,3}, k=3          4,685  (skipped: above the cap)
    ...
    checking stored results in data/sweep_bipartite.json
       K4,5|2         raw     2,256 cells  F_2 Betti [1, 13, 28]  ok
       7 verified, 48 skipped as too large

Read the witness for K_5 as: H_1(B_2 K_5) has one Z/2 (degree 1), and that same Z/2 shows up
again in degree 2 through the Tor term, although H_2 itself is zero.

USAGE
-----
    python3 experiments/crosscheck_mod_p.py [--cap CELLS]

    --cap CELLS   (default 200000) skip any case, named or stored, whose RAW complex has more
                  than CELLS cells.  The raw complex is much larger than the reduced one (K_7 at
                  k = 3: 46,704 raw cells against 25,256 reduced), and the reference pipeline
                  `gbt.swiatkowski` that builds it is slow, so the cap matters.

Nothing is written to data/; the stored sweep files are only read.  With the default cap the
run takes about two minutes and well under 1 GB; `--cap 3000` takes a few seconds.

Exit code: 0 if every prediction holds, 1 if any case mismatches.

PITFALLS
--------
* A torsion summand whose order is prime to both test primes is invisible to this check.  In
  pass (1) the large prime is 1000003, so in practice only 2-torsion and the odd torsion of the
  named cases are tested.
* `betti_mod_p` with the large prime is exact only if 1000003 divides no invariant factor; for
  the complexes here that is certain (all invariant factors are tiny), but it is an assumption.
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import (complete, complete_bip, octahedron, prism,            # noqa: E402
                 swiatkowski, cell_count, paths)
from gbt import fast                                                  # noqa: E402
from gbt.linalg import betti_mod_p                                     # noqa: E402

BIG_PRIME = 1000003


def predicted_betti(H, p):
    """dim H_d(X; F_p) predicted from the integral homology, by universal coefficients.

    H is a homology dict {d: (free rank, list of invariant factors > 1)}.  For each degree d the
    prediction is  free_d + #{q in tors_d : p | q} + #{q in tors_{d-1} : p | q}:  an invariant
    factor divisible by p contributes one Z/p to H_d(-) (x) F_p, and one to Tor(H_{d-1}, F_p).
    Returns {d: predicted dimension} over the same degrees as H.
    """
    out = {}
    for d in sorted(H):
        free, tors = H[d]
        below = H.get(d - 1, (0, []))[1]
        out[d] = (free
                  + sum(1 for q in tors if q % p == 0)
                  + sum(1 for q in below if q % p == 0))
    return out


def named_cases():
    """The fixed list of small cases checked in pass (1), as (label, Graph, k) triples.

    Chosen to cover: Ko--Park torsion in H_1 (K_5, K_{3,3}), a torsion-free planar control (the
    triangular prism at k = 4), the planar Z/2 in H_2 (octahedron at k = 4), complete graphs on
    either side of the first H_2 torsion (K_6, K_7), and odd torsion (K_{4,5} at k = 3, which
    is the only case here where the large prime sees a difference from the integral free rank).
    """
    return [
        ("K_5, k=2", complete(5), 2),
        ("K_{3,3}, k=2", complete_bip(3, 3), 2),
        ("K_{3,3}, k=3", complete_bip(3, 3), 3),
        ("prism, k=4", prism(), 4),
        ("octahedron, k=3", octahedron(), 3),
        ("octahedron, k=4", octahedron(), 4),
        ("K_6, k=3", complete(6), 3),
        ("K_6, k=4", complete(6), 4),
        ("K_7, k=3", complete(7), 3),
        ("K_{4,5}, k=3", complete_bip(4, 5), 3),
    ]


def main():
    """Run pass (1) on the named cases and pass (2) on the stored sweeps; return the exit code."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap", type=int, default=200_000,
                    help="skip cases whose RAW complex exceeds this many cells")
    args = ap.parse_args()
    failures = []

    print("integral homology (Morse + SNF)  vs  mod-p Betti numbers (raw complex)")
    print("p = 2 and p = %d; the difference between them counts the 2-torsion.\n"
          % BIG_PRIME)
    hdr = ("%-16s %10s %-26s %-26s %s"
           % ("case", "raw cells", "F_2 Betti", "F_%d Betti" % BIG_PRIME, "match"))
    print(hdr)
    print("-" * len(hdr))

    for name, G, k in named_cases():
        raw = cell_count(G, k)
        if raw > args.cap:
            print("%-16s %10s  (skipped: above the cap)" % (name, format(raw, ",")))
            continue
        t0 = time.time()
        # the integral answer, through the usual pipeline (fully reduced complex)
        H, _, _ = fast.compute(G, k)
        # and the independent one, from the raw unreduced complex
        C = swiatkowski(G, k)
        b2 = betti_mod_p(C, 2)
        bq = betti_mod_p(C, BIG_PRIME)

        p2 = predicted_betti(H, 2)
        pq = predicted_betti(H, BIG_PRIME)
        ok2 = all(b2.get(d, 0) == p2.get(d, 0) for d in set(b2) | set(p2))
        okq = all(bq.get(d, 0) == pq.get(d, 0) for d in set(bq) | set(pq))
        if not (ok2 and okq):
            failures.append(name)

        # print a Betti dict {d: dim} as "[b_0 b_1 ...]" in increasing degree
        fmt_b = lambda b: "[" + " ".join(str(b.get(d, 0)) for d in sorted(b)) + "]"
        print("%-16s %10s %-26s %-26s %s   [%.0fs]"
              % (name, format(raw, ","), fmt_b(b2), fmt_b(bq),
                 "ok" if (ok2 and okq) else "MISMATCH", time.time() - t0), flush=True)
        if not (ok2 and okq):
            print("      predicted F_2 %s, got %s" % (fmt_b(p2), fmt_b(b2)))
            print("      predicted F_q %s, got %s" % (fmt_b(pq), fmt_b(bq)))
        # (F_2 Betti) - (F_q Betti) in degree d = 2-primary summands of H_d plus those of H_{d-1}
        diff = {d: b2[d] - bq[d] for d in b2 if b2[d] != bq.get(d, 0)}
        if diff:
            print("      2-torsion witness (F_2 minus F_q, by degree): %s" % diff)
    print("-" * len(hdr))

    # ------------------------------------------- stored sweep results, if any
    for fname in ("sweep_complete.json", "sweep_bipartite.json"):
        stored = paths.load_json(fname)
        if not stored:
            continue
        print("\nchecking stored results in data/%s" % fname)
        n_ok = n_skip = 0
        for key, rec in sorted(stored.items()):
            H = {int(d): (v[0], v[1]) for d, v in rec["H"].items()}
            if not any(H[d][1] for d in H):
                continue                       # nothing to cross-check
            # keys are "K{n}|{k}" (complete) or "K{a},{b}|{k}" (complete bipartite)
            spec, k = key.split("|")
            k = int(k)
            if fname.startswith("sweep_complete"):
                G = complete(int(spec[1:]))
            else:
                a, b = (int(x) for x in spec[1:].split(","))
                G = complete_bip(a, b)
            if cell_count(G, k) > args.cap:
                n_skip += 1
                continue
            C = swiatkowski(G, k)
            b2 = betti_mod_p(C, 2)
            p2 = predicted_betti(H, 2)
            ok = all(b2.get(d, 0) == p2.get(d, 0) for d in set(b2) | set(p2))
            if ok:
                n_ok += 1
            else:
                failures.append(key)
                print("   %-14s MISMATCH: predicted %s, got %s" % (key, p2, b2))
            print("   %-14s raw %9s cells  F_2 Betti %s  %s"
                  % (key, format(cell_count(G, k), ","),
                     [b2.get(d, 0) for d in sorted(b2)], "ok" if ok else "MISMATCH"),
                  flush=True)
        print("   %d verified, %d skipped as too large" % (n_ok, n_skip))

    print()
    if failures:
        print("%d MISMATCH(ES): %s" % (len(failures), failures))
        return 1
    print("ALL CROSS-CHECKS PASSED: the mod-p Betti numbers of the raw complexes agree")
    print("with the integral answers from the reduce-then-SNF pipeline in every case.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
