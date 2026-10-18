#!/usr/bin/env python3
"""
An independent check of every torsion claim: Betti numbers over F_p from the RAW complex.

Everything else in this repository goes through the same two steps -- unit-pivot (discrete
Morse) reduction, then exact Smith normal form.  Both are places where a bug would be
invisible: a reduction that forgets to delete a cancelled generator produces a complex whose
homology is simply wrong, and that happened during this project.  So it is worth having a
route to the same information that shares none of that machinery.

THE CHECK
---------
For a prime p, dim H_d(X; F_p) is computed straight from the boundary matrices by sparse
elimination -- no Morse reduction, no SNF, and over F_2 not even division, just xor on bitsets.
Then universal coefficients gives the comparison: if

        H_d(X; Z)  =  Z^{f_d}  (+)  (+)_i Z/q_{d,i}

then

        dim H_d(X; F_p)  =  f_d  +  #{i : p | q_{d,i}}  +  #{i : p | q_{d-1,i}} .

The last term is the Tor contribution from one degree down.  So the integral answer PREDICTS
the mod-p Betti numbers, and this script tests the prediction.  Comparing p = 2 against a large
prime (1000003, standing in for Q) also isolates the 2-torsion on its own: the difference
between the two Betti vectors counts exactly the 2-primary summands in that degree and the one
below.

The raw complex is used deliberately -- the one that has had nothing done to it.

    python3 experiments/crosscheck_mod_p.py [--cap CELLS]

It checks the small named cases below over F_2 and a large prime, and then the stored sweep
results in data/ under the cap over F_2 (which is what detects the 2-torsion).  Exits non-zero
if a prediction fails.
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import (complete, complete_bip, cube, octahedron, prism,      # noqa: E402
                 swiatkowski, cell_count, paths)
from gbt import fast                                                  # noqa: E402
from gbt.linalg import betti_mod_p                                     # noqa: E402

BIG_PRIME = 1000003


def predicted_betti(H, p):
    """dim H_d(X; F_p) predicted from the integral homology, by universal coefficients."""
    out = {}
    for d in sorted(H):
        free, tors = H[d]
        below = H.get(d - 1, (0, []))[1]
        out[d] = (free
                  + sum(1 for q in tors if q % p == 0)
                  + sum(1 for q in below if q % p == 0))
    return out


def named_cases():
    return [
        ("K_5, k=2", complete(5), 2),
        ("K_{3,3}, k=2", complete_bip(3, 3), 2),
        ("K_{3,3}, k=3", complete_bip(3, 3), 3),
        ("prism, k=4", prism(), 4),
        ("octahedron, k=3", octahedron(), 3),
        ("octahedron, k=4", octahedron(), 4),
        ("cube Q_3, k=4", cube(), 4),
        ("K_6, k=3", complete(6), 3),
        ("K_6, k=4", complete(6), 4),
        ("K_7, k=3", complete(7), 3),
        ("K_{4,5}, k=3", complete_bip(4, 5), 3),
    ]


def main():
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

        fmt_b = lambda b: "[" + " ".join(str(b.get(d, 0)) for d in sorted(b)) + "]"
        print("%-16s %10s %-26s %-26s %s   [%.0fs]"
              % (name, format(raw, ","), fmt_b(b2), fmt_b(bq),
                 "ok" if (ok2 and okq) else "MISMATCH", time.time() - t0), flush=True)
        if not (ok2 and okq):
            print("      predicted F_2 %s, got %s" % (fmt_b(p2), fmt_b(b2)))
            print("      predicted F_q %s, got %s" % (fmt_b(pq), fmt_b(bq)))
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
