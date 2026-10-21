Explicit generators of the large torsion, and where it lives (October 19 to 20, 2026)
===================================================================================

Gamma_16 is K_4 with the edges 01, 02, 03, 12, 13 doubled and 23 single, and n = 5. Its torsion is
    T = Tors H_2(Conf_5 Gamma_16; Z) = (Z/2)^7 + (Z/8)^4 + Z/16,   |T| = 2^23.
All programs named below are in certificates/ or analysis/ (see the main README, Sections 7 and 8).

graphs/                       input files ("V E n" + edges)
  K4dbl5_n5.txt               Gamma_16 in the edge order (0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(0,1),(0,2),(0,3),(1,2),(1,3)
  K4dbl5_Hfirst_n5.txt        Gamma_16 with the triangle-with-leaf H listed first:
                              (0,1),(1,2),(1,2),(1,3),(1,3),(2,3) | (0,1),(0,2),(0,2),(0,3),(0,3)
  triangle_leaf_n5.txt        H alone, engine format;  triangle_leaf_analysis_format.txt  H for analysis/decomp
  doubledK4_n5.txt            the doubled K_4 (Gamma_16 plus a second copy of 23, added last)

1. Generators
-------------
generators_K4dbl5_n5.txt.gz   "certificates/certq 2 2 4 -1 gens.txt < graphs/K4dbl5_n5.txt": twelve integral 2-cycles
                              z_0..z_11 generating T (seven of order 2, four of order 8, one of order 16), blocks
                              "Z id order" + "coef cell".  83,000 to 549,000 cells each.  (Written before certq also
                              exported the chains w_t and cocycles; rerunning the command adds them.)
vertex_deletion_none.json     orderq on these generators with nothing deleted, as a self-test.  Each z_t has exactly its
                              stated order and together they generate T.

2. Images of subgraphs
----------------------
For a subgraph H (with the base half-edges of its vertices inside H), [z] is in the image of H_2(Conf_5 H) exactly when
z lies in C_2(H) + B_2(Gamma), i.e. when z vanishes in M = C_2 / (C_2(H) + B_2).  certificates/orderq computes the
coordinates of the z_t in M; kernel.py enumerates T and prints the kernel of T -> M, which is the torsion that comes
from H.  Vanishing mod 2^K is exact once M has no 2-torsion of order >= 2^K, which is checked against the rank of the
same matrix over F_1000003 (minimal_*_rankcheck.log).
vertex_deletion_v2.json, _v3  deleting vertex 2 or 3 leaves only (Z/2)^6 of T in the image.
support_search_n5.jsonl       the downward search (support_search.py): 90 sub-multigraphs up to the symmetries (01), (23)
                              of Gamma_16, with T cap image for each; 40 contain an element of order 16.  Deleting one copy
                              of 01 or of 02 keeps all of T.  Deleting the single edge 23 keeps (Z/2)^7 + (Z/4)^4 + Z/16,
                              although that subgraph (the doubled K_4 minus an edge) has torsion-free H_2.
minimal_A_K10.json, minimal_E_K10.json, minimal_*_rankcheck.log
                              two minimal cases redone mod 2^10, with the rank check.
image_triangle_leaf_Hfirst_K10.json
                              the same test for the triangle with a leaf in the H-first edge order (generators made in that
                              order); T cap image = (Z/2)^4 + (Z/8)^4 + Z/16, containing z_2 + z_11 of order 16.
image_in_doubledK4.json       the generators pushed into the doubled K_4.  z_11 (order 16) goes to an element of order 4
                              with a unit coordinate on the Z/4 summand, so the Z/4 of the doubled K_4 is the image of the
                              Z/16.

The minimal subgraphs whose image contains an element of order 16 (all with six edges and torsion-free H_2(Conf_5)):
  type  edges (multiplicities)       contains 23   T cap image                      b_2(Conf_5 H)
  (a)   01, 12^2, 13^2, 23           yes           (Z/2)^4 + (Z/8)^4 + Z/16         4374
  (a)   01^2, 03, 12, 13^2           no            (Z/2)^4 + (Z/4)^4 + Z/16         4374
  (b)   01^2, 03^2, 12, 13           no            (Z/2)^3 + (Z/4)^4 + Z/16         5734
  (c)   01, 03, 12, 13^2, 23         yes           (Z/2)^3 + (Z/8)^4 + Z/16         3204
  (c)   01^2, 02, 03, 12, 13         no            (Z/2)^4 + (Z/4)^4 + Z/16         3204
  (d)   01^2, 03, 12, 13, 23         yes           (Z/2)^3 + (Z/8)^4 + Z/16         1974
  (a) triangle with two doubled edges, leaf at the vertex where they meet; (b) the same triangle, leaf at a vertex where
  a doubled edge meets the single edge; (c) K_4 minus an edge, opposite edge doubled; (d) K_4 minus an edge, adjacent
  edge doubled.  Every edge of Gamma_16 is needed for order 16 in Gamma_16 itself.  Deleting a copy of 01, a copy of 02,
  or 23 lowers the largest order to 4, to 8, or to no torsion at all.

3. The explicit cycle on the triangle with a leaf
-------------------------------------------------
H = edges (0,1),(1,2),(1,2),(1,3),(1,3),(2,3) of graphs/K4dbl5_Hfirst_n5.txt.
target_x_z2_plus_z11.txt.gz    x = z_2 + z_11 (generators made in the H-first order), an element of order 16 of the image.
triangle_leaf_cycle_y.txt.gz   y with x = y + d u exactly over Z, from
                               "SOLVER=certificates/solve_residual.py certificates/exactsolve 2 6 x.txt sol < graphs/K4dbl5_Hfirst_n5.txt"
                               (about 10 minutes, 6 GB; residual block 1274 x 1055 of rank 3).  40,839 cells, coefficients
                               up to 6.8e8, every cell on H.  The solution depends on pivot choices: the current exactsolve
                               (which drops cells using an elimination mod a large prime) returns a different valid y
                               with 39,986 cells; validation/check_generator_tools.py builds and verifies its certificate.
triangle_leaf_cycle_y_homology_coords.txt.gz
                               y reduced modulo boundaries of Conf_5(H) by certificates/reducecyc: 2323 nonzero coordinates
                               in a basis of H_2(Conf_5 H) = Z^4374.
certificates/cert_triangle_leaf_into_K4dbl5_order16.txt.gz
                               the certificate (SUBGRAPH 6): W = 8 w_2 + w_11 - 16 u with d W = 16 y, the mod-16 cocycle of
                               z_11, and f(y) = 1 mod 16.  So [y] has infinite order in Conf_5(H) and its image in
                               Conf_5(Gamma_16) has order exactly 16.

4. The class is a sum of tori
-----------------------------
tori/decomp_p2.json, decomp_p1000003.json
                               analysis/decomp on H: the 54,540 tori (a star class at an essential vertex x times a 1-cycle of
                               H minus x; with three essential vertices this covers every product of disjoint 1-classes) span
                               rank 4369 of H_2(Conf_5 H) = Z^4374, the same mod 2 and mod 1000003.
tori_triangle_leaf_n5.txt.gz   the tori as integral chains ("DUMP=file analysis/decomp 5 1000003 0 < graphs/triangle_leaf_analysis_format.txt").
tori/y_block_six_edges.txt.gz  y with its cells written for the six edges of H.
tori/y_mod_tori_*.json, .log   "certificates/orderq 2 p K tori/y_block_six_edges.txt - - tori_triangle_leaf_n5.txt < graphs/triangle_leaf_n5.txt"
                               for (p, K) = (2, 8), (3, 1), (5, 1), (1000003, 1), (998244353, 1).  M = H_2 / tori is free of
                               rank 5 with no 2-torsion, and y vanishes in it in every case.  So y is in the rational span of
                               the tori (modular evidence) and then in their span over Z_(2).  So the class of order 16 is a sum of
                               images of tori.  An explicit integral combination has not been computed.

Scripts
-------
support_search.py   the downward search (usage in its docstring).
kernel.py           T cap image from an orderq output, by enumerating T.

5. The class of order 16 as a sum of two tori (October 20, 2026)
----------------------------------------------------------------
tori/ (scripts and results; graphs/K4dbl5_Hfirst_n5.txt is the edge order used)
  minimal_tori.py           greedy search for a small set of tori of H whose span (modulo the boundaries of Gamma_16)
                            contains an element of order 16, with certificates/reducemod and certificates/spanq (mod 2^24).
  minimal_tori_search.*     its result: 14 tori, all with the star at vertex 1, in two families.
  exact_check.py            the exact test: rational relations among the images (two large primes, rational
                            reconstruction) and the mod-16 cocycle on them.
  two_tori_exact_check.json the unique relation among the 14 images, coefficients +-1, cocycle value 7 (odd): order 16.
  two_tori_cycle_ystar.txt  the cycle y* = s_24 x c_24 + s_34 x c_34 (653 cells, coefficients +-1, +-2), certified by
                            certificates/cert_two_tori_into_K4dbl5_order16.txt.gz (16 y* = d W exactly, from
                            "certificates/exactsolve 2 0 ...").
  describe_tori.py, walks.py  read a torus back as two loops (closed walks of single moves through vertices).
  shortest_loops.py         searches the simple loops of the 1-skeleton of Conf_3(H minus vertex 1) for short loops
                            homologous to the second factors (normal forms modulo boundaries, two primes).
  second_factor_loops_L12.txt, _L16.txt   c_24 = two loops of 12 moves; c_34 is not a single loop or a sum of two loops of
                            length <= 16, and is written as its three closed walks (12, 30, 12 moves).
  simplify_piece.py, second_factor_exchange_ranks.txt
                            H_1(Conf_3(H minus vertex 1)) has rank 40; two-particle exchanges at vertices 2 and 3 with
                            the third particle parked span only rank 24, and neither second factor lies in their span.
The tori:  s_24 and s_34 are 8-move exchanges of particles {2,4} and {3,4} at vertex 1 (valence 5) using the leaf edge 01,
both copies of 12 and one copy of 13.  The second factors are loops of the other three particles in the H-shaped tree
H minus vertex 1 (vertices 2 and 3 with two hairs each, joined by 23); the particles cross the bridge 23 and return in a
different order.
