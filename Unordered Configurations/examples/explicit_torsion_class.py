#!/usr/bin/env python3
"""
explicit_torsion_class.py -- write down a torsion class by hand, prove its order, and push it
into a bigger graph.

WHAT IT IS FOR
--------------
Computing H_1(B_2 K_5) = Z^6 + Z/2 tells you that a class of order 2 exists; it does not tell
you what it is.  This script carries out, on the smallest interesting example, the three steps
that the experiments use on larger ones (octahedron_generator.py, theta_class.py,
two_torsion_anatomy.py):

  1. HALVING.  Find an integral chain x of degree d+1 whose boundary is divisible by 2 but not
     zero.  Then z = (d x)/2 is an integral cycle (d z = d d x / 2 = 0) and 2 z = d x is a
     boundary, so [z] has order 1 or 2.  Such an x is an element of ker(d_{d+1} mod 2) that does
     not lift to ker(d_{d+1}) over Z.
  2. ADJOIN AND TEST.  [z] is nonzero exactly when z is not a boundary.  Compare the torsion of
     coker(d_{d+1}) with the torsion of coker(d_{d+1} | z), i.e. with z added as an extra
     column: adjoining z divides out the subgroup it generates, so the Z/2 disappears exactly
     when [z] generates it.  (By the identity tors coker d_{d+1} = tors H_d, the cokernel is
     all we need.)
  3. PUSH FORWARD.  K_5 is a subgraph of K_6 on the first five vertices.  The inclusion of
     graphs induces a chain map of Swiatkowski complexes that simply renames edges, so z gives
     a cycle i_*(z) of K_6, which we test the same way in K_6.

HOW TO READ THE OUTPUT
----------------------
    complex of K_5 at k = 2: {0: 115, 1: 280, 2: 160} generators
    torsion of H_1 (from coker d_2): Z/2
    ker(d_2 mod 2) has dimension 1; trying random elements of it
      try 1: d x has 100 nonzero entries, all even -> z = d x / 2 is a 1-cycle: True
      torsion after adjoining z: 0   => [z] generates the Z/2
    pushed into K_6: z uses 100 cells of the complex of K_6 at k = 2; is a cycle: True
      torsion of H_1(B_2 K_6): Z/2;  after adjoining i_*(z): 0  => i_*[z] generates it

and it prints the first few terms of z in the notation of look_at_the_complex.py.

WHY RANDOM ELEMENTS
-------------------
The map x -> [d x / 2] from ker(d mod 2) onto H_d[2] is linear and onto, so half of all
kernel elements work -- but the elements of a basis produced by elimination can all land in
the wrong coset (this happened for the octahedron).  Random sums of basis vectors avoid that.
Here the kernel is one-dimensional, so the first try works.

    python3 examples/explicit_torsion_class.py          (well under a second; prints only)
"""

import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import complete, primary, swiatkowski                              # noqa: E402
from gbt.linalg import coker_torsion, gf2_kernel, indexed_columns           # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from look_at_the_complex import describe_cell                               # noqa: E402

K = 2           # particles
DEGREE = 1      # we look for torsion in H_1, so we work with d_2 (degree-2 chains x)


def boundary_of(cols, x_bits):
    """d x for x given as a bitmask over the columns of d (bit j set = column j with
    coefficient 1).  Returns {row index: integer coefficient}, zeros removed."""
    out = {}
    for j, col in enumerate(cols):
        if x_bits >> j & 1:
            for r, c in col.items():
                out[r] = out.get(r, 0) + c
    return {r: c for r, c in out.items() if c}


def is_cycle(C, d, chain):
    """True if the degree-d chain {generator: coefficient} has zero boundary in C."""
    total = {}
    for g, a in chain.items():
        for f, c in C.B[d].get(g, {}).items():
            total[f] = total.get(f, 0) + a * c
    return not any(total.values())


def main():
    """Find the Z/2 of H_1(B_2 K_5) explicitly, check it, and push it into K_6."""
    rng = random.Random(1)
    G = complete(5)
    C = swiatkowski(G, K)                 # unreduced: every generator has a readable meaning
    print("complex of K_5 at k = %d: %s generators" % (K, C.sizes()))

    # Columns of d_2 indexed by the degree-2 generators; rows by the degree-1 generators.
    d2, row_index, _ = indexed_columns(C, DEGREE + 1)
    rows = {i: g for g, i in row_index.items()}            # row number -> generator
    base_torsion = coker_torsion(d2, verbose=False)
    print("torsion of H_1 (from coker d_2): %s" % (primary(base_torsion)))

    kernel = gf2_kernel(d2)                # basis of ker(d_2 mod 2), as bitmasks of columns
    print("ker(d_2 mod 2) has dimension %d; trying random elements of it" % len(kernel))
    z = None
    for attempt in range(1, 21):
        x = 0
        for b in kernel:                   # a uniformly random element of the F_2-span
            if rng.random() < 0.5:
                x ^= b
        dx = boundary_of(d2, x)
        if not dx:
            continue                       # x lifts to an integral cycle: no information
        assert all(c % 2 == 0 for c in dx.values()), "x was not in the kernel mod 2"
        candidate = {r: c // 2 for r, c in dx.items()}
        chain = {rows[r]: c for r, c in candidate.items()}
        print("  try %d: d x has %d nonzero entries, all even -> z = d x / 2 is a 1-cycle: %s"
              % (attempt, len(dx), is_cycle(C, DEGREE, chain)))
        after = coker_torsion(d2 + [candidate], verbose=False)
        verdict = "[z] generates the Z/2" if not after else "[z] is zero; try again"
        print("  torsion after adjoining z: %s   => %s" % (primary(after), verdict))
        if not after:
            z = chain
            break
    if z is None:
        print("no generator found in 20 tries")
        return 1

    print("\nthe first terms of z (coefficient, generator):")
    for g, c in sorted(z.items(), key=lambda kv: repr(kv[0]))[:6]:
        print("  %+d  %s" % (c, describe_cell(G, g)))

    # Push forward along K_5 -> K_6 (vertices 0..4 are shared).  A generator (mono, st) of K_5
    # becomes (mono', st') of K_6: rename each edge index to the index of the same edge in K_6,
    # and give the new vertex 5 the state 0 (empty).  No signs change, because the vertex
    # order -- which fixes the Koszul signs -- is the same.
    H = complete(6)
    emap = [H.edge_index(u, v) for (u, v) in G.edges]

    def push(cell):
        """The image of a K_5 generator in the complex of K_6."""
        mono, st = cell
        new_mono = tuple(sorted(emap[e] for e in mono))
        new_st = tuple(2 + emap[s - 2] if s >= 2 else s for s in st) + (0,)
        return (new_mono, new_st)

    pushed = {push(g): c for g, c in z.items()}
    D = swiatkowski(H, K)
    d2H, row_index_H, _ = indexed_columns(D, DEGREE + 1)
    column = {row_index_H[g]: c for g, c in pushed.items()}
    print("\npushed into K_6: z uses %d cells of the complex of K_6 at k = %d; is a cycle: %s"
          % (len(pushed), K, is_cycle(D, DEGREE, pushed)))
    before = coker_torsion(d2H, verbose=False)
    after = coker_torsion(d2H + [column], verbose=False)
    print("  torsion of H_1(B_2 K_6): %s;  after adjoining i_*(z): %s  => %s"
          % (primary(before), primary(after),
             "i_*[z] generates it" if before and not after else "i_*[z] does not generate it"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
