#!/usr/bin/env python3
"""
Is the H_2 2-torsion of B_3(K_n) a product of smaller classes?  For K_8, yes -- entirely.

THE QUESTION
------------
A product class in H_2(B_3 Gamma) comes from a pair of VERTEX-DISJOINT subgraphs
G_1, G_2 with k_1 + k_2 = 3 particles, through

        B_{k_1}(G_1) x B_{k_2}(G_2)  ->  B_3(Gamma),

and by Kunneth the only ways to get torsion in degree 2 are
H_i(B_{k_1} G_1) (x) H_j(B_{k_2} G_2) with i+j = 2 and torsion in one factor, or
Tor(H_1, H_1).  The term H_2 (x) H_0 would need H_2 torsion with at most 2 particles on a
proper subgraph, which never occurs in the small-graph sweeps, so the useful term is
H_1 (x) H_1.  Torsion in H_1(B_k G) needs G non-planar, so at least five vertices
(Ko--Park); a nonzero H_1(B_1 G) = H_1(G) needs a cycle, so at least three.  A product
therefore needs at least EIGHT vertices, and K_8 is the first complete graph where one exists:
G_1 = K_5 on a 5-subset with two particles, G_2 = the complementary triangle with one.

THE TEST
--------
Build every such product explicitly as a 2-cycle of the Swiatkowski complex of K_8, adjoin
all of them to the image of d_3, and watch what happens to the torsion:

        coker(d_3)                 = (Z/2)^21        the torsion of H_2(B_3 K_8)
        coker(d_3 | all products)  = ?

If the second is trivial, the products generate the whole (Z/2)^21 and every torsion class in
H_2(B_3 K_8) is decomposable.  That is what happens, with 9,296 product columns.

This is NOT the general picture: `kunneth_K6.py` shows K_6 at k=4 has no product classes
available at all, and K_7 at k=3 and the whole of the seven-vertex sweep are indecomposable
for the same counting reason.  The conjectured threshold is 2k+2 vertices (so 8 for k=3), and it
matches exactly where the complete-graph 2-torsion multiplicity starts tracking rank H_1.

THE METHOD, STEP BY STEP
------------------------
  1. Build the UNREDUCED weight-3 Swiatkowski complex of K_n (`gbt.swiatkowski`) and the
     columns of d_3, indexed by the degree-2 generators (`indexed_columns`).  The torsion of
     coker(d_3) is the torsion of H_2 (the identity "torsion is a cokernel").
  2. Compute a Z-basis of the 1-cycles Z_1 of the weight-2 complex of K_5 (`kernel_basis`).
     These are ALL 1-cycles, free classes and the Ko--Park class alike; using the whole cycle
     lattice is harmless, since adjoining boundaries changes nothing.
  3. For every 5-subset S of the vertices of K_n, every triangle in the complement, and every
     basis cycle of Z_1(K_5): copy the K_5 cycle onto S (`translate`), multiply it with the
     fundamental 1-cycle of the triangle carrying one particle (`triangle_cycle`, `merge`), and
     record the result as a column indexed by the degree-2 generators of K_n.
  4. Check on 50 sampled products that d(product) = 0, and abort if not.
  5. Compare the torsion of coker(d_3) with that of coker(d_3 | all products).

The product construction needs care with signs: merging a basis element of G_1 with one of
G_2 reorders the moving slots, so the merge carries the Koszul sign (-1)^(number of
inversions between the two slot sets).  Omitting it produces chains that are not cycles, which
is why step 4 exists.

WHAT IS PRINTED
---------------
    === K_n, k = 3 ===
       weight-3 complex {degree: number of generators}
          d_3 alone: C cols -> peeled P, residual R x S     (progress of `coker_torsion`)
       baseline: torsion of coker(d_3) = ...                the torsion of H_2(B_3 K_n)
       Z_1(K_5, weight 2): G generators  [time]
       built N product columns
       sanity: of 50 sampled products, B fail d(product) = 0
       torsion of coker(d_3 | products) = ...

then a summary line "baseline X -> with products adjoined Y" and one of three verdicts:
"The products kill the torsion completely ... DECOMPOSABLE", "kill some but not all", or
"kill none ... INDECOMPOSABLE".  For n < 8 no products exist and the script stops after
"built 0 product columns".  Real output with --n 7 (half a second):

    === K_7, k = 3 ===
       weight-3 complex {0: 3864, 1: 15624, 2: 19656, 3: 7560}
          d_3 alone: 7560 cols -> peeled 7559, residual 572x1
       baseline: torsion of coker(d_3) = Z/2
       Z_1(K_5, weight 2): 166 generators  [0s]
       built 0 product columns
       no products exist here, so the torsion is indecomposable by default.

USAGE
-----
    python3 experiments/decomposability_K8.py [--n 8]

    --n N   (default 8) work in K_N at k = 3.  N < 8 prints a note and finds no products;
            N = 8 is the interesting case; N >= 9 is much more expensive (the number of
            products grows like C(N,5) C(N-5,3)).

Takes about twenty minutes for K_8 and several GB; building the unreduced complex of K_8 and
the final cokernel are the expensive parts.  --n 7 takes under a second.  Nothing is written
to data/.

WHAT THE TEST PROVES, AND WHY ADJOINING FREE PRODUCTS IS HARMLESS
-----------------------------------------------------------------
Adjoining columns computes the torsion of A / P, where A = coker(d_3) and P is the subgroup
generated by the products.  The adjoined products include products of FREE classes (step 2
uses all of Z_1), and one might worry that a free column could make torsion disappear without
the torsion lying in P.  It cannot: if t in A has finite order m, its image in A / P also has
order dividing m, so it is a torsion element of A / P.  Hence

    tors(A / P) = 0   implies   every torsion element of A lies in P.

(For instance, in Z (+) Z/2 the quotient by the free element (1, 1) is Z/2, not 0: the
relation e + t = 0 forces 2e = 0.)  So "the torsion disappears" proves that every torsion class
of H_2(B_3 K_n) is a Z-linear combination of product classes, modulo boundaries.  It does not say
that each torsion class is a single product, only that the products span the torsion.

Exit code: 1 if the sign sanity check fails (the test would be meaningless), 0 otherwise --
including when the torsion turns out indecomposable, which is a result, not an error.
"""

import argparse
import itertools
import os
import sys
import time
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import complete, swiatkowski, primary                      # noqa: E402
from gbt.linalg import coker_torsion, indexed_columns, kernel_basis  # noqa: E402


def edge_index_map(G):
    """{frozenset({u, v}): edge index} for the graph G, so an edge can be looked up unordered."""
    return {frozenset(e): i for i, e in enumerate(G.edges)}


def translate(gen, vmap, edges5, eidx, N):
    """A K_5 basis element -> the corresponding K_N element under the vertex map `vmap`.

    gen       a basis element (mono, st) of the Swiatkowski complex of K_5 (see gbt.swiatkowski)
    vmap      a sequence of 5 distinct vertices of K_N: vertex v of K_5 goes to vmap[v]
    edges5    the edge list of K_5, so that edge index e of K_5 is the pair edges5[e]
    eidx      {frozenset({u, v}): edge index in K_N}
    N         number of vertices of K_N
    Returns the basis element of K_N obtained by relabelling vertices and re-indexing edges.
    The parked monomial is re-sorted because edge indices change order under the relabelling.
    No sign arises: the moving slots of a single factor keep their relative order only if vmap
    is increasing, which holds because vmap is a sorted 5-subset (itertools.combinations).
    """
    mono, st = gen
    nst = [0] * N
    for v5, x in enumerate(st):
        if x == 0:
            continue
        if x == 1:
            nst[vmap[v5]] = 1
        else:
            nst[vmap[v5]] = 2 + eidx[frozenset(vmap[w] for w in edges5[x - 2])]
    nmono = tuple(sorted(eidx[frozenset(vmap[w] for w in edges5[e])] for e in mono))
    return (nmono, tuple(nst))


def moving_slots(st):
    """The vertices whose state is a moving generator h_{v,e} (state code >= 2), increasing."""
    return [v for v, x in enumerate(st) if x >= 2]


def merge(b1, b2):
    """Tensor two basis elements on disjoint vertex sets, with the Koszul sign.

    b1, b2 are basis elements (mono, st) of the same complex, supported on disjoint vertex sets
    (so their state tuples can simply be added).  The product lists the moving generators of
    b1 first and then those of b2; the basis element lists them in increasing vertex order, and
    reordering costs (-1)^(number of pairs x in b1, y in b2 with x > y).
    Returns (basis element, sign).
    """
    m1, s1 = b1
    m2, s2 = b2
    st = tuple(a + b for a, b in zip(s1, s2))
    M1, M2 = moving_slots(s1), moving_slots(s2)
    inv = sum(1 for x in M1 for y in M2 if x > y)
    return (tuple(sorted(m1 + m2)), st), (-1) ** inv


def triangle_cycle(tri, eidx, N):
    """The fundamental 1-cycle of the weight-1 complex on the triangle `tri`.

    With one particle the configuration space of a triangle is the triangle itself, and the
    cycle is  sum over its three edges {u, w} of (h_{u,e} - h_{w,e}): each term has boundary
    (e - u) - (e - w) = w - u, and the three boundaries cancel around the triangle.
    Returns {basis element of K_N: coefficient}, with the basis elements of K_N's complex.
    """
    a, b, c = tri
    out = {}
    for (u, w) in ((a, b), (b, c), (c, a)):
        e = eidx[frozenset((u, w))]
        for (x, s) in ((u, 1), (w, -1)):
            st = [0] * N
            st[x] = 2 + e
            key = ((), tuple(st))
            out[key] = out.get(key, 0) + s
            if out[key] == 0:
                del out[key]
    return out


def main():
    """Run the five steps of the module docstring for K_n; return the exit code."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=8, help="work in K_n")
    args = ap.parse_args()
    n = args.n
    if n < 8:
        print("NOTE: products need at least 8 vertices (see the docstring); with n = %d"
              % n)
        print("there are none, and the baseline torsion is automatically indecomposable.")

    print("=== K_%d, k = 3 ===" % n, flush=True)
    G = complete(n)
    eN = edge_index_map(G)
    C = swiatkowski(G, 3)
    print("   weight-3 complex %s" % C.sizes(), flush=True)

    # the rows of d_3 are the degree-2 generators; index2 maps a generator to its row
    cols3, index2, _ = indexed_columns(C, 3)
    base = cols3
    baseline = coker_torsion(base, "d_3 alone")
    print("   baseline: torsion of coker(d_3) = %s" % primary(baseline), flush=True)

    # --- the cycle lattice of K_5 at weight 2 -----------------------------
    G5 = complete(5)
    edges5 = list(G5.edges)
    C5 = swiatkowski(G5, 2)
    cols1_5, _, gens1_K5 = indexed_columns(C5, 1)
    t0 = time.time()
    ker5 = kernel_basis(cols1_5, len(C5.gens[0]))
    print("   Z_1(K_5, weight 2): %d generators  [%.0fs]"
          % (len(ker5), time.time() - t0), flush=True)

    # --- every product over every 5-subset and complementary triangle -----
    # A product column is (K_5 cycle on S) x (triangle cycle on tri), written in the degree-2
    # generators of K_n.  Products whose terms are not degree-2 generators of K_n cannot occur
    # (S and tri are disjoint), but `index2.get` guards against it anyway.
    prods = []
    for S in itertools.combinations(range(n), 5):
        rest = [v for v in range(n) if v not in S]
        if len(rest) < 3:
            continue
        for tri in itertools.combinations(rest, 3):
            z2 = triangle_cycle(tri, eN, n)
            for kv in ker5:
                col = {}
                for j, cf in kv.items():
                    b1 = translate(gens1_K5[j], S, edges5, eN, n)
                    for b2, cf2 in z2.items():
                        mg, sg = merge(b1, b2)
                        r = index2.get(mg)
                        if r is None:
                            continue
                        col[r] = col.get(r, 0) + cf * cf2 * sg
                        if col[r] == 0:
                            del col[r]
                if col:
                    prods.append(col)
    print("   built %d product columns" % len(prods), flush=True)
    if not prods:
        print("   no products exist here, so the torsion is indecomposable by default.")
        return 0

    # --- sanity: a product must be a cycle --------------------------------
    cols2, _, _ = indexed_columns(C, 2)   # cols2[r] is d_2 of the r-th 2-generator
    bad = 0
    for col in prods[:50]:
        acc = defaultdict(int)
        for r, x in col.items():
            for rr, xx in cols2[r].items():
                acc[rr] += x * xx
        if any(acc.values()):
            bad += 1
    print("   sanity: of 50 sampled products, %d fail d(product) = 0" % bad, flush=True)
    if bad:
        print("   ABORTING: the Koszul merge sign is wrong, so the test is meaningless.")
        return 1

    # --- adjoin and see what survives -------------------------------------
    after = coker_torsion(base + prods, "d_3 + products")
    print("   torsion of coker(d_3 | products) = %s" % (primary(after)), flush=True)

    print()
    print("baseline %s   ->   with products adjoined %s"
          % (primary(baseline), primary(after)))
    if baseline and not after:
        print("The products kill the torsion completely: every torsion class in")
        print("H_2(B_3 K_%d) is a Kunneth product of a Ko--Park class on a K_5 with a" % n)
        print("cycle on a disjoint triangle.  The 2-torsion is DECOMPOSABLE.")
    elif after != baseline:
        print("The products kill some but not all of the torsion: %s survives."
              % primary(after))
    else:
        print("The products kill none of the torsion: it is INDECOMPOSABLE.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
