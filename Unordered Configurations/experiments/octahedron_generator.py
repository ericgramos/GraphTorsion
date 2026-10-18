#!/usr/bin/env python3
"""
An explicit generator of the octahedron's Z/2 in H_2(B_4), and where it goes in K_6.

THE HALVING CONSTRUCTION
------------------------
To produce a torsion class by hand, take a chain x with d(x) = 0 mod 2 but d(x) != 0 over Z.
Then every coefficient of d(x) is even, so

        z := d(x) / 2

is defined, satisfies d(z) = 0 (because d(d(x)) = 0 and we only divided), and has 2z = d(x) a
boundary.  So [z] is a class of order dividing 2.  Finding such an x means finding an element
of ker(d_3 mod 2) that is NOT in the image of ker(d_3 over Z), and the map

        lambda : ker_{F_2}(d_3)  ->  H_2[2] = Z/2,       x |-> [d(x)/2]

is linear and onto, so half of the mod-2 kernel works.

Whether a particular z generates is then decided by the adjoin-and-test idiom: compare
coker(d_3) with coker(d_3 | z).  If the Z/2 dies, z generated it.

A TRAP WORTH RECORDING
----------------------
"Half of the kernel works" does not mean half of a BASIS works.  The mod-2 kernel basis
produced by echelon elimination is highly structured, and in this case the first fourteen
basis vectors and the first twelve prefix sums all landed in the wrong coset -- lambda = 0 for
every one of them.  A random element of the kernel works with probability 1/2, so a few random
tries suffice (with the default seed the first one already works).  Hence the three strategies below; `random` is the default because it is the
only one that is not fighting the structure of the basis.

    random    xor a uniformly random subset of the kernel basis; retry
    prefix    try the prefix sums in order (included to document that it fails here)
    bisect    bisect on prefix sums to locate a single good basis vector, which is
              informative when one exists but costs a cokernel computation per step

PUSHING THE CLASS INTO K_6
--------------------------
The octahedron is a subgraph of K_6, and K_6 at k=4 also has exactly Z/2 in H_2.  Are they the
same class?  The inclusion induces a chain map (relabel edge indices; the vertex labels are
already shared), so push z forward and test it in K_6 by the same idiom.  The answer is yes:
i_*(z) generates the K_6 torsion.  That is the second instance (after theta_class.py) of
torsion created by relations rather than carried by a generator.

    python3 experiments/octahedron_generator.py [--strategy random|prefix|bisect] [--seed S]

Takes a few minutes with the default strategy.
"""

import argparse
import os
import random
import sys
import time
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import complete, octahedron, swiatkowski, primary          # noqa: E402
from gbt.linalg import coker_torsion, gf2_kernel, indexed_columns   # noqa: E402

K = 4


def halve(cols, mask):
    """z = d(x)/2 for x the column subset encoded by `mask`, or None if d(x) = 0 over Z."""
    acc = defaultdict(int)
    j = 0
    m = mask
    while m:
        if m & 1:
            for r, v in cols[j].items():
                acc[r] += v
        m >>= 1
        j += 1
    acc = {r: v for r, v in acc.items() if v}
    if not acc:
        return None                      # x was an integral cycle: lambda is not defined
    assert all(v % 2 == 0 for v in acc.values()), "d(x) is not even: x is not a mod-2 cycle"
    return {r: v // 2 for r, v in acc.items()}


def is_cycle(z, cols_down):
    """Check d(z) = 0, as a guard against an indexing slip."""
    acc = defaultdict(int)
    for r, v in z.items():
        for rr, vv in cols_down[r].items():
            acc[rr] += v * vv
    return not any(acc.values())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--strategy", choices=("random", "prefix", "bisect"),
                    default="random")
    ap.add_argument("--seed", type=int, default=4)
    ap.add_argument("--tries", type=int, default=12)
    args = ap.parse_args()

    OCT = octahedron()
    CO = swiatkowski(OCT, K)
    print("octahedron, weight-%d complex %s" % (K, CO.sizes()), flush=True)
    # indexed_columns returns (columns, {row generator: row index}, column generators).
    # The rows of d_3 are the degree-2 generators, which is where z lives.
    d3, index2, _ = indexed_columns(CO, 3)
    d2, _, gens2 = indexed_columns(CO, 2)       # gens2[r] is the r-th degree-2 generator

    baseline = coker_torsion(d3, "octahedron d_3")
    print("baseline torsion of coker(d_3) = %s" % primary(baseline), flush=True)
    if not baseline:
        print("no torsion to find; aborting")
        return 1

    t0 = time.time()
    kern = gf2_kernel(d3)
    print("dim ker(d_3) over F_2 = %d   [%.0fs]" % (len(kern), time.time() - t0),
          flush=True)

    # ---------------------------------------------------------------- search
    found = None
    if args.strategy == "random":
        rng = random.Random(args.seed)
        for trial in range(args.tries):
            mask = 0
            for t in kern:
                if rng.random() < 0.5:
                    mask ^= t
            z = halve(d3, mask)
            if z is None:
                continue
            res = coker_torsion(d3 + [z], "random trial %d" % trial)
            print("  trial %d: |supp z| = %d -> torsion %s"
                  % (trial, len(z), primary(res)), flush=True)
            if res != baseline:
                found = z
                print("  *** generator found on trial %d ***" % trial, flush=True)
                break
    elif args.strategy == "prefix":
        mask = 0
        for j, t in enumerate(kern, start=1):
            mask ^= t
            z = halve(d3, mask)
            if z is None:
                continue
            res = coker_torsion(d3 + [z], "prefix %d" % j)
            print("  prefix %d: torsion %s" % (j, primary(res)), flush=True)
            if res != baseline:
                found = z
                print("  *** generator found at prefix %d ***" % j, flush=True)
                break
            if j >= args.tries:
                print("  (stopping after %d prefixes -- see the docstring: the prefix"
                      % args.tries)
                print("   sums are exactly where this search does NOT work)", flush=True)
                break
    else:   # bisect
        def lam(j):
            mask = 0
            for t in kern[:j]:
                mask ^= t
            z = halve(d3, mask)
            if z is None:
                return 0, None
            res = coker_torsion(d3 + [z], "prefix %d" % j, verbose=False)
            return (1 if res != baseline else 0), z

        hi = None
        probes = sorted({len(kern), len(kern) // 2, len(kern) // 4, len(kern) // 8,
                         len(kern) // 16, 30, 15, 7, 3, 1}, reverse=True)
        for j in probes:
            v, _ = lam(j)
            print("  lambda(prefix %d) = %d" % (j, v), flush=True)
            if v == 1:
                hi = j
                break
        if hi is None:
            print("  no prefix with lambda = 1 among those probed")
        else:
            lo = 0
            while hi - lo > 1:
                mid = (lo + hi) // 2
                v, _ = lam(mid)
                print("  bisect: lambda(prefix %d) = %d" % (mid, v), flush=True)
                if v == 1:
                    hi = mid
                else:
                    lo = mid
            print("  *** basis vector #%d is a generator ***" % (hi - 1), flush=True)
            found = halve(d3, kern[hi - 1])

    if found is None:
        print("\nno generator located with strategy %r." % args.strategy)
        if args.strategy == "prefix":
            print("That is the expected outcome, and the point of including it:")
            print("run with --strategy random instead.")
            return 0
        return 1

    z = found
    print("\ngenerator z = d(x)/2, supported on %d cells of C_2(octahedron)" % len(z))
    print("  d(z) = 0 : %s" % is_cycle(z, d2))
    print("  adjoining z to d_3 gives torsion %s (was %s)"
          % (primary(coker_torsion(d3 + [z], "oct + z")), primary(baseline)), flush=True)

    # ------------------------------------------------- push forward into K_6
    print("\n=== pushing i_*(z) into K_6 ===", flush=True)
    K6 = complete(6)
    eK6 = {frozenset(e): i for i, e in enumerate(K6.edges)}
    eOCT = {i: frozenset(e) for i, e in enumerate(OCT.edges)}
    CK = swiatkowski(K6, K)
    print("K_6 weight-%d complex %s" % (K, CK.sizes()), flush=True)
    d3k, index2k, _ = indexed_columns(CK, 3)
    base_k = coker_torsion(d3k, "K_6 d_3")
    print("K_6 baseline torsion = %s" % primary(base_k), flush=True)

    def push(gen):
        """An octahedron basis element -> the K_6 element with the same vertex states."""
        mono, st = gen
        nst = [0] * 6
        for v, x in enumerate(st):
            if x == 1:
                nst[v] = 1
            elif x >= 2:
                nst[v] = 2 + eK6[eOCT[x - 2]]
        return (tuple(sorted(eK6[eOCT[e]] for e in mono)), tuple(nst))

    zk = {index2k[push(gens2[r])]: v for r, v in z.items()}
    print("i_*(z) supported on %d cells of C_2(K_6)" % len(zk), flush=True)
    res_k = coker_torsion(d3k + [zk], "K_6 + i_*(z)")
    print("K_6 with i_*(z) adjoined: %s" % primary(res_k), flush=True)
    print()
    if res_k != base_k:
        print("=> i_*(z) GENERATES the K_6 torsion: the Z/2 of K_6 at k=4 is the image")
        print("   of the octahedron's class under the inclusion of a planar subgraph.")
    else:
        print("=> i_*(z) is a BOUNDARY in K_6: the K_6 class is not pushed forward from")
        print("   the octahedron, and its origin is still open.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
