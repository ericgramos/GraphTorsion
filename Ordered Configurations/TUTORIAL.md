# Tutorial: computing, proving and locating torsion in graph configuration spaces

This tutorial is a guided path through the package for a reader who knows the mathematics
(configuration spaces, cellular homology, Smith normal forms) and is comfortable typing
commands into a terminal and editing a short Python script, but is not a programmer. Each
lesson says what question it answers, gives commands that work as they stand from the top
folder of the package, shows what the output looks like, and explains how to read it. The
README is the reference manual; this file is the walk-through. Times and memory are for one
core of an ordinary computer with 8 GB of memory.

Contents

1. Setting up
2. Graph files, cells and the other file formats
3. Homology of Conf_n(Γ) from Python
4. Large cases: the local Smith form
5. Proving that torsion exists: certificates
6. Explicit generators of the torsion
7. Where does a class live? Images of subgraphs
8. Writing a class explicitly on a subgraph
9. Is a class a sum of tori?
10. Checking the software yourself
11. Practical advice and common mistakes

---

## 1. Setting up

    make                      # builds everything; takes about a minute
    make check                # the fast validation suite; about half an hour

`make` needs a C++17 compiler, and the Python parts need `python-flint` (and `networkx`
for a few scripts). Section 1 of the README explains installation on Linux and macOS.

Every program prints its results on standard output (usually one line of JSON, a format
that Python reads with `json.loads`) and progress messages on standard error. A JSON object
looks like `{"betti": [1, 18, 161], "torsion": [[], [], []]}`: names in quotes, lists in
square brackets.

## 2. Graph files, cells and the other file formats

**Graphs.** A graph is a text file with `V E n` on the first line (vertices, edges,
particles) and then one line `u v` per edge. Vertices are numbered 0, ..., V−1, particles
0, ..., n−1. Parallel edges are repeated lines; a loop is `v v`. The doubled K_4 with five
particles is `examples/graphs/dk4_5.in`. The programs in `analysis/` take the same file
without the particle number (first line `V E`); they take n on the command line.

**Half-edges and bases.** Edge e = (u, v) has two half-edges: 2e at its tail u and 2e+1 at
its head v. The base half-edge of a vertex is its lowest-numbered half-edge. So **the order
of the edges matters**: it decides the bases, hence which cells exist. This becomes
important in Lessons 7 and 8, where a subgraph must be listed first.

**Cells.** A cell of the reduced Świątkowski complex places every particle either at a
vertex, together with a non-base half-edge h there (standing for the generator h − h0), or
inside an edge, in an order from tail to head. Its degree is the number of occupied
vertices. In files a cell is written

    V[2/7,-,0/5,-] E[1 3||4||||]

meaning: particle 2 sits at vertex 0 with half-edge 7, vertex 1 is empty, particle 0 sits at
vertex 2 with half-edge 5, vertex 3 is empty; edge 0 carries particles 1 then 3 (from tail
to head), edge 2 carries particle 4, the other edges are empty.

**Chains and cycles.** A chain is a list of lines `coefficient cell`. A file of several
chains uses blocks, each starting with a header line `Z ...`:

    Z 11 16
    3 V[...] E[...]
    -1 V[...] E[...]

**Certificates** (Lesson 5) are text files, often gzipped, with lines `GRAPH V E n`,
`EDGE u v`, `DEGREE i`, `PRIME P`, optionally `SUBGRAPH h`, then `W coef cell` (a chain of
degree i+1) and `F coef cell` (a cochain of degree i).

## 3. Homology of Conf_n(Γ) from Python

The simplest way in is the function `homology` of `python/gconf.py`:

    cd python
    python3 -c "from gconf import homology; print(homology([(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)], 4))"

The result is a dictionary. The two entries that matter most are

    'betti':   [1, 18, 161, 0, 0]        # ranks of H_0, H_1, H_2, ...
    'torsion': [[], [], [], [], []]      # elementary divisors > 1 of each torsion group

so H_*(Conf_4 K_4; Z) = Z, Z^18, Z^161 with no torsion. The entries `chi` and `gal` are the
Euler characteristic of the complex and Gal's formula; they must agree. The first torsion
example:

    python3 -c "from gconf import homology; print(homology([(0,1),(0,1),(0,2),(0,2),(1,2),(1,2)], 6)['torsion'])"
    # [[], [], [2, 2], []]      H_2(Conf_6 T(2,2,2)) has torsion (Z/2)^2; two minutes, 1 GB

What happens underneath: the engine `engine/oconf` builds every boundary map, eliminates it
over Z using only pivots ±1, and hands the small leftover block to FLINT for an exact Smith
form (README Section 4.5). Useful variants:

    homology(edges, n, prime=2)              # dimensions over F_2 instead of Z
    homology(edges, n, prime=1000003)        # over a large prime field: the rational Betti numbers
    homology(edges, n, colors="aabcd")       # particles 0 and 1 indistinguishable (Conf_n / S_2)

By the universal coefficient theorem, dim H_i(F_p) − b_i counts the cyclic p-primary summands
of H_i and H_{i−1} together. Comparing `prime=2` with `prime=1000003` is the cheapest torsion
test there is.

## 4. Large cases: the local Smith form

Exact elimination over Z can stall when the leftover entries grow. For anything large, ask
for the local Smith form over Z/p^K instead:

    engine/oconf -q 2 3 < examples/graphs/dk4_5.in        # doubled K_4, n = 5; 5 minutes, 1.5 GB

Each boundary map in the output has a list `stages`: s_0 pivots that are units, s_1
elementary divisors exactly 2, s_2 exactly 4, and so on up to p^(K−1). For d_3 of the doubled
K_4 you see one 4 and eight 2's, which is the torsion (Z/2)^8 ⊕ Z/4 of H_2.

**The one thing to remember.** Elementary divisors divisible by p^K are invisible mod p^K:
they are simply zero. To rule them out, compare the total number of pivots of each map with
its rank over a large prime (`engine/oconf -p 1000003`). They agree exactly when nothing is
hidden. Every torsion statement in this package was checked this way (README Section 15).

`python/gconf.py` has a function `local_torsion(edges, n, p, k)` that runs this and reports
the orders found. Group quotients (README Section 5) make much larger cases feasible, but they
see only part of the homology, so read that section before using them.

## 5. Proving that torsion exists: certificates

A computation says that H_i has an element of order P. A **certificate** proves it in a form
that anyone can check with a short program that shares no code with the engine:

* an integral chain w of degree i+1 whose boundary is divisible by P, so that z = ∂w/P is an
  integral cycle with Pz a boundary;
* a cocycle f modulo P with f(z) prime to p.

If (P/p)·z were a boundary ∂u, then (P/p)·f(z) = f(∂u) = 0 mod P, so p would divide f(z).
Hence [z] has order exactly P. (For P = p this just says that z is not a boundary.)

To check one:

    python3 certificates/verify_cert.py certificates/cert_T222_n6.txt.gz

The verifier prints each check as it passes: (a) every cell named exists; (b) ∂w is divisible
by P; (c) ∂z = 0; (d) f vanishes mod P on the boundary of every cell of degree i+1 (all of
them are enumerated, and their number is compared with the closed formula); (e) f(z) is
prime to p. The last line begins with `VERIFIED`. A wrong or corrupted certificate fails at
one of these steps.

To make one:

    # prime order p (here 3, for T(1,3,3) at n = 6): certp degree prime seed outfile
    certificates/certp 2 3 5 my_cert.txt < examples/graphs/t133_6.in

    # prime-power order p^k (here 4 = 2^2, the doubled K_4 at n = 5): certq degree p k index outfile
    certificates/certq 2 2 2 0 cert4.txt < examples/graphs/dk4_5.in      # 6 minutes, 3 GB
    python3 certificates/verify_cert.py cert4.txt                          # 1 minute

`certq` eliminates modulo p^(k+1), computes the local Smith form of the leftover block with its
transformations, picks a diagonal entry p^k·unit (the `index`-th one, counting from 0), and
turns the transformations back into w and f. It prints the local Smith form it found, which
is a useful report in itself (`1:8 2:1` means eight entries of valuation 1 and one of
valuation 2).

## 6. Explicit generators of the torsion

With index −1, `certq` writes one integral cycle for **every** cyclic summand of order below
p^(k+1):

    G=data/survey/k4_family/support/graphs
    certificates/certq 2 2 4 -1 gens.txt < $G/K4dbl5_n5.txt          # 4 minutes, 3 GB

For K_4 with five doubled edges (torsion (Z/2)^7 ⊕ (Z/8)^4 ⊕ Z/16) this writes twelve blocks
`Z t order`. Each holds the cycle z_t (`coef cell` lines), a chain w_t with ∂w_t = order·z_t
(`W` lines), and, for the summands of the top order, a cocycle mod p^k detecting it (`F`
lines). These z_t generate the 2-primary torsion: every torsion class is Σ c_t [z_t] for
unique c_t mod order_t. They are large (up to half a million cells each); the next lessons
are about finding out what they are.

## 7. Where does a class live? Images of subgraphs

Let H ⊂ Γ be a subgraph (possibly with fewer copies of a multiple edge). If the base
half-edge of every vertex of H lies in H, the cells of Γ that use only edges and half-edges
of H form a subcomplex, which is the complex of H, and its inclusion models
Conf_n H → Conf_n Γ. A cycle z of Γ has its class in the image of H_i(Conf_n H) exactly when
z = y + ∂u with y supported on H. The program `orderq` tests this for a list of cycles:

    certificates/orderq 2 2 10 gens.txt 5 - < $G/K4dbl5_n5.txt > img.json

Arguments: degree, prime, K (work mod p^K), the cycle file, the edges to delete (comma
separated, here edge 5 = the single edge 23), and base overrides (`-` for none; needed only
when the lowest half-edge of a vertex of H lies on a deleted edge, written `v:h`). The output
lists, for every input cycle, its coordinates in the torsion of
M = C_i / (B_i + C_i(H)) and the order of its image there (`order_in_M`). A cycle of the
torsion comes from H exactly when its image is 0. To see which torsion classes come from H:

    python3 data/survey/k4_family/support/kernel.py img.json

This enumerates the whole torsion group (2^23 elements here, a few seconds) and prints the
subgroup T ∩ image, for example `{2: 7, 4: 4, 16: 1}` for (Z/2)^7 ⊕ (Z/4)^4 ⊕ Z/16.

**Exactness.** orderq works mod p^K. Its conclusions are exact once M has no torsion of order
p^K or more. Check it by running the same command with a large prime and K = 1 (for example
`orderq 2 1000003 1 ...`) and comparing: unit pivots plus the number of diagonal entries mod
p^K must equal the rank over the large prime (both are printed on standard error).

The search of `data/survey/k4_family/support/support_search.py` runs this over many
subgraphs; its README describes the results.

## 8. Writing a class explicitly on a subgraph

Once orderq says that a class comes from H, `exactsolve` finds the cycle on H. List the edges
of H first, so that its bases lie inside it (the file `$G/K4dbl5_Hfirst_n5.txt` does this for
the triangle with a leaf, H = the first 6 edges), make generators in that edge order
(Lesson 6), put the chosen combination of generators in a file of `coef cell` lines, and run

    SOLVER=$PWD/certificates/solve_residual.py certificates/exactsolve 2 6 x.txt sol < $G/K4dbl5_Hfirst_n5.txt

It writes `sol.y` (the cycle y on H, exactly over Z) and `sol.u` (the chain u with x = y + ∂u).
About 10 minutes and 6 GB here. To make a certificate for the image of y, put
W = (the chains w_t of the combination, scaled to order 16) − 16·u, the cocycle F of the
order-16 generator, `PRIME 16` and `SUBGRAPH 6` into a file in the format of Lesson 2; the
verifier then also checks that z = ∂W/16 lives on the first six edges. The certificate
`certificates/cert_triangle_leaf_into_K4dbl5_order16.txt.gz` was made this way, and
`validation/check_generator_tools.py` (step `exact`) contains a few lines of Python that
assemble such a certificate. With h = 0 (empty subgraph) `exactsolve` simply decides whether
x is a boundary and writes x = ∂u; that is how the certificate of Lesson 9 was made.

`certificates/reducecyc` reduces a cycle of a graph whose maps reduce with ±1 pivots modulo
boundaries, giving its coordinates in a basis of H_i.

## 9. Is a class a sum of tori?

A torus is a product of disjoint 1-classes, for example a star class at a vertex x (particles
exchanging at x) times a 1-cycle of Γ minus x. `analysis/decomp` lists them and measures their
span; with `DUMP=file` it also writes each torus as a chain, with a header saying which vertex,
which particles and which basis cycles it is made of:

    DUMP=tori.txt analysis/decomp 5 1000003 0 < $G/triangle_leaf_analysis_format.txt

To test whether a cycle lies in the span of the tori, give them to orderq as relations (an
optional seventh argument): `orderq 2 2 8 y.txt - - tori.txt < graph`. A cycle is in
(boundaries + span of the tori) exactly when its image vanishes.

For many such questions (which tori are really needed?) use the fast path. `reducemod`
reduces all cycles modulo the boundaries once (15 minutes), and `spanq` then answers each
query in seconds:

    certificates/reducemod 2 2 24 blocks.txt blocks.red < $G/K4dbl5_Hfirst_n5.txt
    certificates/spanq blocks.red        # then type lines "span 3 17 40 ..." (relation numbers)

The search script `data/survey/k4_family/support/tori/minimal_tori.py` drives spanq. **spanq
works mod p^K and can be fooled** when the images of the tori are divisible by high powers of
p, so every final answer must be confirmed exactly with `tori/exact_check.py`. That script
finds the rational relations among the images of the chosen tori (modulo two large primes,
with rational reconstruction) and evaluates the order-16 cocycle on them: an odd value means
order 16. Then make a certificate as in Lesson 8. Finally, `tori/describe_tori.py` and
`tori/walks.py` print each torus as two closed walks of single moves ("particle 2 passes
through vertex 1 from edge 01 to edge 12"), and `tori/shortest_loops.py` looks for short
loops homologous to a given 1-cycle. The README of `data/survey/k4_family/support/` records
what came out: the class of order 16 is a sum of two tori.

## 10. Checking the software yourself

| command | what it checks | time |
|---|---|---|
| `make check` | the engine against FLINT, Abrams' model and the vertex explosion; quotients; dense phases; local Smith forms; certificates | 30 min |
| `python3 validation/check_generator_tools.py` | every computation behind the explicit generators: certificates, generators, orderq, exactsolve, reducecyc, the torus span, and spanq against orderq on random sets | 2 hours |
| `g++ -O3 -fopenmp -std=c++17 -o oconf26 history/engine_versions/oconf26.cpp` and then `python3 validation/regression_vs_reference.py engine/oconf ./oconf26` | the engine against its reference version on 33 tests (every option and mode) | 2 min |

Logs of our runs are in `validation/logs/`.

## 11. Practical advice and common mistakes

* **Memory is the limit.** The number of cells is n!·e_k(d_v − 1)·C(n − k + E − 1, E − 1)
  (`engine/oconf -c` prints it). Around 10 million cells per degree needs a few GB. Exit code
  137 means the system ran out of memory.
* **Mod p^K hides large divisors.** Always compare with a large prime (Lesson 4).
* **Edge order decides the bases.** For subgraph questions list the subgraph first; otherwise
  use base overrides in orderq.
* **Fast answers are not proofs.** spanq and the local modes find candidates; certificates
  and the independent verifier prove.
* **Do not kill programs with `pkill -f name`** from a shell whose own command line contains
  the same word: it kills the shell too. Use the process id (`kill 12345`).
* **Big intermediate files.** Generator files (Lesson 6) can be hundreds of MB; compress them
  with gzip when you keep them.
