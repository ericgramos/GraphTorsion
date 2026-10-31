# gconf-torsion

Software for computing the homology of ordered configuration spaces of graphs, with
particular attention to torsion, together with every script, data file and certificate
from a study of the question whether these homology groups can have torsion.

The short version of what was found is the following. The homology of the ordered
configuration space Conf_n(Γ) of a graph Γ does have torsion, and it already happens for
planar graphs. The smallest examples we know are the triangle multigraphs T(a,b,c) at six
particles, for instance

    H_2(Conf_6 T(2,2,2); Z) = Z^71848 ⊕ (Z/2)^2
    H_2(Conf_6 T(1,3,3); Z) = Z^224725 ⊕ (Z/3)^6

These triangle examples were the first torsion found. We later found torsion in Conf_5 of the doubled K_4, with H_2 = Z^83615 ⊕ (Z/2)^8 ⊕ Z/4.
Related doublings of K_4 have torsion of order 8 and 16, and for those classes the package
contains certificates of exact order, explicit generators, and a cycle on a six-edge subgraph
with torsion-free configuration space (a triangle with two doubled edges and a leaf) whose
image has order exactly 16 (Section 8, and `data/survey/k4_family/support/`). A long
search for torsion in configuration spaces of non-planar graphs (K_5, K_{3,3} and many
graphs built from them) found none in any case that fitted in 8 GB of memory.

If you are new to the package, start with `TUTORIAL.md`, a guided path through every tool with
worked examples. This README explains how to install the package, how to use the engine and the Python
library, and what every file in the package does. It is long because it is meant to be
usable by someone who has never seen the code. Sections 1 to 3 are enough to start
computing. Section 4 is the reference for the engine's options. Section 5 explains group
quotients, which need some care. Sections 6 to 13 go through the remaining folders one
file at a time.

## Contents

1. Installation
2. A first computation
3. The mathematical model in brief
4. The engine `engine/oconf`: every option, with examples
5. Group quotients and characters
6. The Python library `python/`
7. The analysis tools `analysis/`
8. Certificates `certificates/`
9. Validation `validation/`
10. Surveys and campaign scripts `surveys/`
11. Data `data/`
12. Example inputs `examples/graphs/`
13. History `history/`
14. Reproducing the main results
15. Limits, pitfalls and practical advice
16. References

---

## 1. Installation

The package needs

* a C++17 compiler with OpenMP support. On Linux, `g++` version 9 or later works. On macOS,
  install `g++` with Homebrew (`brew install gcc`) and build with `make CXX=g++-14` (use the
  version number Homebrew installed), or use Apple clang together with `libomp`;
* Python 3.8 or later;
* the Python package `python-flint`, which provides exact Smith normal forms over the
  integers. Install it with `pip install python-flint` (on some systems
  `pip install --user python-flint` or `pip install --break-system-packages python-flint`);
* optionally `networkx` (`pip install networkx`), used only by `python/toric.py`;

To build everything, run `make` in the top folder:

    cd gconf-torsion
    make            # builds engine/oconf, certificates/cert, certp, certq, orderq, exactsolve, reducecyc, analysis/*
    make check      # runs the fast validation scripts (about 25 minutes)

If you only want the engine, `make engine` is enough. The Makefile passes
`-O3 -march=native -fopenmp -std=c++17`. The flag `-march=native` makes the binary specific to
the processor it was built on; remove it from `CXXFLAGS` in the Makefile if you build on one
machine and run on another. Without OpenMP the code still compiles if you remove `-fopenmp`,
but the dense phase and the parallel matrix construction then run on one core.

The Python code finds the engine at `engine/oconf` relative to its own location. To use a
binary somewhere else, set the environment variable `OCONF_BIN` to its path. Temporary files
(the residual matrices written by the engine in exact mode) go to the folder named by the
environment variable `GCONF_SCRATCH`, and to the system temporary folder if that is not set.
Some residual files can be large, so point `GCONF_SCRATCH` at a disk with space if you run big
exact computations.

## 2. A first computation

The engine reads a graph from standard input. The first line is `V E n` (number of vertices,
number of edges, number of particles), followed by one line `u v` per edge. Vertices are
numbered from 0. Edge e = (u, v) is oriented from u to v; the orientation only fixes a
convention for the cells and does not change the space.

The complete graph K_4 with four particles:

    printf "4 6 4\n0 1\n0 2\n0 3\n1 2\n1 3\n2 3\n" | engine/oconf -p 2

The engine prints progress to standard error and a single JSON object to standard output:

    {"V":4,"E":6,"n":4,"P":2,"dim":4,"cells":[3024,10752,12096,4608,384],
     "maps":[{"k":1,"rows":10752,"cols":3024,"pivots":3023,...}, ...]}

`cells` lists the number of cells in each degree, and each entry of `maps` describes one
boundary map d_k, with `pivots` its rank over the chosen coefficients. The Betti number in
degree k is then cells[k] minus the rank of d_k minus the rank of d_{k+1}.

Doing that bookkeeping by hand is tedious, so normally you call the engine from Python:

    cd python
    python3 -c "from gconf import homology; print(homology([(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)], 4))"

which prints a dictionary containing

    'betti': [1, 18, 161, 0, 0], 'torsion': [[], [], [], [], []], 'chi': 144, 'gal': 144, ...

So H_0 = Z, H_1 = Z^18, H_2 = Z^161, no torsion, and the Euler characteristic of the complex
(144) agrees with Gal's formula (144). The `torsion` entry lists, for each degree i, the
elementary divisors larger than 1 of the torsion subgroup of H_i. For the first torsion
example in the triangle family:

    python3 -c "from gconf import homology; print(homology([(0,1),(0,1),(0,2),(0,2),(1,2),(1,2)], 6)['torsion'])"
    # [[], [], [2, 2], []]      that is, H_2 has torsion (Z/2)^2; about two minutes, 1 GB of memory

## 3. The mathematical model in brief

This section gives the definitions needed to use the software.

**Configuration spaces.** For a finite graph Γ (loops and multiple edges allowed), Conf_n(Γ)
is the space of n-tuples of pairwise distinct points of Γ, and UConf_n(Γ) = Conf_n(Γ)/S_n is
the unordered version. Both are aspherical, and their fundamental groups are the pure and
the full graph braid groups (Abrams 2000; Ghrist 2001).

**The reduced ordered Świątkowski complex.** Świątkowski (2001) built a small cubical complex
for the unordered space. The ordered version that the engine uses is written out by Hainaut,
Knudsen and Wawrykow (2025, Definition 3.12 and Theorem 3.14). Its cells are described as
follows.

* Each edge e = (u, v) has two half-edges. Half-edge 2e is the end of e at its tail u, and
  half-edge 2e+1 the end at its head v.
* Each vertex has a chosen base half-edge h0 (by default the lowest-numbered half-edge at
  the vertex).
* A cell places each of the n labelled particles either at a vertex v, where it carries one
  of the reduced generators h − h0 for a half-edge h ≠ h0 at v, or in the interior of an edge.
  At most one particle sits at each vertex. The particles inside an edge are ordered from
  tail to head.
* The degree of a cell is the number of particles sitting at vertices.
* The boundary of a cell is the alternating sum over its occupied vertices (in increasing
  vertex order) of "push the particle into the edge of h at the h end" minus "push it into
  the edge of h0 at the h0 end".

The number of cells of degree k is n! · e_k(d_v − 1) · C(n − k + E − 1, E − 1), where e_k is
the k-th elementary symmetric function of the numbers d_v − 1 over all vertices (d_v is the
valence) and C is a binomial coefficient. This grows fast. For six particles on a graph with
three vertices and seven edges there are already about 13 million cells.

**What "torsion" means computationally.** The engine computes the rank of each boundary map
and, when asked, its elementary divisors (exactly over Z), its rank over a prime field F_p,
or its local Smith form over Z/p^K, which records how many elementary divisors are exactly
p, p², and so on. The torsion of H_{k−1} is read off from the elementary divisors of d_k.

## 4. The engine `engine/oconf`: every option, with examples

`engine/oconf.cpp` is one C++ file of about 2400 lines; `engine/gr42.h` holds the Galois
ring arithmetic. The top of `oconf.cpp` has a long comment that explains the mathematics,
the data structures and the elimination strategy. Every function has a comment saying what
it does, what it expects and what it returns.

### 4.1 Input

    V E n
    u_0 v_0
    ...
    u_{E-1} v_{E-1}

Limits: at most 10 particles, and 3E < 255 (each half-edge and each edge must be addressable
with 8 bits). Loops `v v` and repeated lines (parallel edges) are allowed. Isolated vertices
are allowed and do nothing.

### 4.2 Output

One JSON object on standard output:

| key | meaning |
|---|---|
| `V`, `E`, `n`, `P` | the input sizes and the prime (0 for exact) |
| `dim` | the top degree of the complex |
| `cells` | list of cell counts by degree |
| `maps` | one object per boundary map d_k |

Each object in `maps` has

| key | meaning |
|---|---|
| `k` | the degree of the map d_k : C_k → C_{k−1} |
| `rows`, `cols` | number of k-cells and number of (k−1)-cells that remained as columns (fewer than cells[k−1] when column dropping applies, see 4.5) |
| `pivots` | number of pivots found by elimination; over a field this is the rank; in exact mode it is the number of unit pivots |
| `res_rows`, `res_cols`, `res_nnz`, `maxabs` | (exact mode) size, number of nonzero entries and largest absolute entry of the residual block that could not be pivoted on with units |
| `res_file` | (exact mode) the file that contains the residual block; the Python driver computes its Smith form with FLINT and deletes it |
| `stages` | (local mode `-q`) list s_0, s_1, s_2, ... where s_j is the number of elementary divisors of d_k exactly equal to p^j (s_0 counts units) |
| `leftover_nnz` | (local mode) number of nonzero entries that remained after the last stage. This is normally 0. It does **not** detect elementary divisors divisible by p^K: those entries are 0 mod p^K and simply vanish. See Section 15 for how to detect them |

With a group (Section 5), the JSON also has a field `group` with the group order and the
chosen base half-edges.

Progress messages go to standard error. They include the number of cells, the time of each
step and, during elimination, a line every 200,000 pivots with the number of live rows and
the average fill per row.

### 4.3 Options

| option | meaning |
|---|---|
| `-p P` | coefficients. `P = 0` (the default) works exactly over Z. `P = 2` uses a fast routine over F_2. For 3 ≤ P ≤ 13 a fast routine over F_P. Any larger prime uses the general sparse routine. P must be prime. |
| `-q p K` | local Smith form over Z/p^K, with p prime and p^K < 32 (so for p = 2, K ≤ 4; p = 3, K ≤ 3; p = 5, K ≤ 2; p ≤ 31, K = 1). Reports `stages` for every map. This is the most reliable way to find p-primary torsion in large complexes, because entries never grow. |
| `-o PREFIX` | prefix for the residual files `PREFIX_k.txt` (default `res`, so the files `res_1.txt`, `res_2.txt`, ... appear in the current folder). They are written in exact mode and with large primes; the Python driver chooses its own prefix in `GCONF_SCRATCH` and deletes the files afterwards. When you call the engine by hand you can delete them. |
| `-d m` | build only d_1, ..., d_m. Then H_0, ..., H_{m−1} are computed completely, including the torsion of H_{m−1}. Useful when the top maps are too big. |
| `-c` | only count cells and exit. |
| `-z` | check that d∘d = 0 on every cell of degree ≥ 2, and exit. A sanity check of the complex (not available together with `-W`). |
| `-T` | with `-p 2`, compute the rank of the top map through its transpose. Sometimes much faster. |
| `-k COLORS` | quotient by a Young subgroup. COLORS is a word with one letter per particle; particles with the same letter become indistinguishable. `-k aaaa` with n = 4 gives UConf_4. |
| `-g GENS` | a group of label permutations, given by generators separated by `;`. Each generator is a comma-separated permutation of 0..n−1; a shorter list is extended by the identity. Example `-g "1,2,3,4,0"` is the 5-cycle on the first five labels. |
| `-a GENS` | a group of graph automorphisms, generators separated by `;`. A generator is `v_0,...,v_{V−1}` (the image of each vertex), optionally followed by `|e_0,...,e_{E−1}` (the image of each edge). An edge image written `~e` means the edge is mapped onto edge e with its orientation reversed. Without the edge part, parallel edges are matched in input order. |
| `-X TW` | a character twist. TW contains `l` for the sign of the label permutation, `e` for the sign of the edge permutation, both (`le`), or `n` for no twist. |
| `-Y s_1,s_2,...` | a sign ±1 for each automorphism generator. The engine checks that these define a character of the group (exit code 12 otherwise). |
| `-W t_1,t_2,...` | Galois ring mode. The coefficients become GR(4,2) = Z/4[ω], ω² + ω + 1 = 0, and automorphism generator i acts by ω^{t_i}. Used for the part of the homology on which an automorphism of order 3 acts by a primitive cube root of unity, at the prime 2. The group must have odd order (exit code 13). |
| `-F` | fast orbit enumeration. Always use it with `-g` or `-a`; it gives the same complex much faster. |
| `-D m` | allow the dense finishing phase in local mode once at most m columns are alive. `-D 0` disables it, `-D -1` starts dense immediately (only for small matrices). |
| `-M a b` | switch to the dense phase as soon as at most a columns are alive, or else once the sparse part holds more than b nonzero entries. Defaults a = 30000, b = 120000000. |
| `-C c` | consider the dense phase only once the cheapest remaining column has at least c entries (default 48). |
| `-G m` | the dense threshold for the Galois ring engine (default 70000). |
| `-m` | generic matrix mode. Read a sparse integer matrix (`rows cols`, then `r c v` triples) from standard input and eliminate it. The Python driver uses this for auxiliary matrices. |

Unknown options are an error.

### 4.4 Examples

    # Betti numbers of Conf_5(K_5) over F_3, only H_0 and H_1 (maps d_1, d_2)
    printf "5 10 5\n0 1\n0 2\n0 3\n0 4\n1 2\n1 3\n1 4\n2 3\n2 4\n3 4\n" | engine/oconf -p 3 -d 2

    # exact 3-primary torsion of Conf_6(T(1,3,3)) up to exponent 2 (stages s_1 count Z/3 summands)
    engine/oconf -q 3 2 < examples/graphs/t133_6.in

    # the 2-torsion of the doubled K_4 at n = 5, in the quotient by the label rotation Z/5
    engine/oconf -q 2 3 -g 1,2,3,4,0 -F < examples/graphs/dk4_5.in

    # unordered configuration space UConf_2(K_5): H_1 = Z^6 ⊕ Z/2 (Abrams)
    printf "5 10 2\n0 1\n0 2\n0 3\n0 4\n1 2\n1 3\n1 4\n2 3\n2 4\n3 4\n" | engine/oconf -k aa -q 2 2

    # check d∘d = 0 (on Conf_4(K_4); for the large examples this takes a few minutes)
    printf "4 6 4\n0 1\n0 2\n0 3\n1 2\n1 3\n2 3\n" | engine/oconf -z

### 4.5 How the engine works, in a few paragraphs

The cells of each degree are enumerated, encoded as 128-bit integers (12 bits per particle:
8 bits for its location, 4 bits for its position on an edge) and stored in sorted arrays, so
that a face of a cell can be found by binary search. The matrix of d_k is built with one row
per k-cell, and reduced by sparse Gaussian elimination with the Markowitz rule (pivot in a
column with the fewest entries, using the shortest row in that column). Over Z only pivots
±1 are used. Whatever is left is a small "residual" block, whose Smith normal form is
computed exactly by FLINT. In local mode only entries prime to p are used as pivots; when
none remains, all remaining entries are divisible by p, they are divided by p, and the
next stage begins. The number of pivots in stage s is the number of elementary divisors
exactly p^s.

Column dropping is the key memory saving. Once d_k is reduced with unit pivots, the k-cells
used as pivot rows can be deleted from the columns of d_{k+1}, because the kernel of d_k
projects isomorphically onto the remaining coordinates and the image of d_{k+1} lies in that
kernel. This roughly halves the second matrix.

When elimination has consumed the easy pivots, fill-in grows quickly (in the logs, the
average row length jumps from about 5 to over 20 within a few hundred thousand pivots once
fewer than about 40,000 columns are alive). At that point the local engine switches to a
dense phase: the remaining rows are streamed through a dense row-reduced basis. For Z/4 and
Z/9 the dense rows are bit-sliced (two or four bit planes, so 64 entries are processed per
machine word), and for GR(4,2) a similar bit-sliced arithmetic is used. Rows that are
divisible by p after reduction ("leftovers") feed the next stage; they are compressed on the
fly against the span of the previous leftovers mod p.

### 4.6 Exit codes

| code | meaning |
|---|---|
| 0 | success |
| 1 | bad input, or graph or n too large |
| 2 | a face of a cell was not found (an internal consistency check) |
| 3 | integer overflow in exact mode (use `-q` instead) |
| 4 | local mode found a non-unit entry not divisible by p (internal check) |
| 5 | unsupported modulus in a dense phase |
| 6, 7 | internal checks of the dense phase |
| 8 | invalid automorphism, or no base half-edges invariant under the automorphisms exist |
| 9 | out of memory in a dense phase |
| 10 | the group order is divisible by p; coinvariants would not compute the right thing |
| 11 | the group does not split as labels × automorphisms (internal check) |
| 12 | the values given by `-Y` or `-W` do not define a character |
| 13 | Galois ring mode needs a group of odd order |
| 137 | (from the shell) the operating system killed the process for using too much memory |

## 5. Group quotients and characters

Large complexes become manageable by passing to a quotient by a group of symmetries. This
section explains what such a quotient computes, because it is easy to misread the output.

Let G be a finite group acting on the complex C = C_*(Conf_n Γ) by permuting particle labels
(the group L given by `-g`) and by graph automorphisms (the group A given by `-a`), so
G = L × A. Let χ : G → {±1} (or, in Galois ring mode, G → ⟨ω⟩) be a character, chosen by
`-X`, `-Y` or `-W`. The engine builds the complex of twisted coinvariants

    C_G^χ = C / ⟨ g·c − χ(g) c ⟩.

Its cells are orbit representatives. A cell whose stabiliser acts on it by a scalar other
than χ becomes zero. The action of a group element on a cell is a signed permutation: the
sign comes from re-sorting the occupied vertices after applying the automorphism (the
boundary formula depends on their order) and from χ.

If |G| is invertible in the coefficient ring R (the engine checks this and exits with code
10 otherwise), averaging over G shows that

    H_*(C_G^χ; R) ≅ H_*(C; R)^χ,

the part of H_*(Conf_n Γ; R) on which G acts through χ. When R is a field or a local ring
Z/p^K with p not dividing |G|, the p-primary torsion of H_*(Conf_n Γ) splits into isotypic
parts in the same way.

The consequence to keep in mind is that a quotient computation sees only one isotypic part.
If it finds torsion, then Conf_n Γ has torsion. If it finds none, the torsion could sit in
another part. For abelian G every irreducible representation over a field containing the
right roots of unity is a character, and running all characters is a complete check. For
non-abelian G (for instance the dihedral or Frobenius groups used in the campaign) the
characters see only the one-dimensional representations; the higher-dimensional ones are
invisible. `validation/check_quotients.py` checks the engine's quotients against an independent computation of the invariant
subcomplex.

The base half-edges must be preserved by the automorphisms. The engine chooses invariant
base half-edges automatically when they exist (it reports them in the `group` field), and
exits with code 8 when they do not.

Examples:

    # tripled C_4 at n = 5, quotient by the label rotation Z/5: H_2 has 3-torsion (Z/3)^3
    engine/oconf -q 3 2 -g 1,2,3,4,0 -F < examples/graphs/tc4_5.in

    # T(1,3,3) with an apex at n = 6: Frobenius group F_20 on labels, two parallel-edge swaps,
    # sign twist on labels. H_2 has (Z/3)^2 in this isotypic part.
    engine/oconf -q 3 2 -g "1,2,3,4,0;0,2,4,1,3" -a "0,1,2,3|0,2,1,3,4,5,6,7,8,9;0,1,2,3|0,1,2,3,5,4,6,7,8,9" -X l -F -D 70000 < examples/graphs/t133a_6.in

## 6. The Python library `python/`

All modules have a docstring at the top that explains the mathematics and the usage, and
every function has its own docstring.

### `gconf.py`

The main driver. Functions:

| function | what it does |
|---|---|
| `homology(edges, n, prime=0, tag=None, keep=False, timeout=None, colors=None)` | Runs the engine and returns a dictionary with `cells`, `betti`, `torsion` (a list per degree of the elementary divisors > 1), `rank`, `chi` (Euler characteristic of the complex), `gal` (Gal's formula) and `residual` (sizes of the residual blocks). With `prime = p` the `betti` entry holds dimensions over F_p and no torsion is reported. `colors` gives a Young quotient as in `-k`. |
| `local_torsion(edges, n, p, k=None, tag=None, colors=None, extra=())` | Exact p-primary torsion through the local Smith form. `k` defaults to the largest exponent with p^k < 32. `extra` passes further engine options (for example a group). |
| `gal_euler(V, edges, n)` | The Euler characteristic of Conf_n from Gal's formula (Gal 2001), used as an independent check of the cell counts. |
| `abrams_homology(edges, n)` | Homology of Abrams' discretised configuration space of the graph (it should be subdivided first with `subdivide`). An independent model written in pure Python, usable for small cases. |
| `subdivide(edges, m)` | Subdivides every edge into m pieces. |
| `normalize(edges)` | Renumbers vertices to 0..V−1 and returns (V, edges). |
| `smooth(edges)` | Removes vertices of valence 2 (homeomorphism class). |
| `snf_divisors(fn)`, `snf_divisors_modp(fn, p)` | Rank and elementary divisors (or rank mod p) of a residual file written by the engine, computed with FLINT. |
| `sparse_snf(nr, nc, triplets)` | Smith form of an arbitrary sparse integer matrix through the engine's `-m` mode. |
| `summarize(json)` | A one-line readable summary of an engine run. |

### `unordered.py`

An independent implementation, in Python, of the reduced unordered Świątkowski complex
(cells are sets of occupied vertices with half-edge generators and numbers of particles per
edge). `uhomology(edges, n)` returns the homology of UConf_n. It reproduces the known
2-torsion of H_1(UConf_n) for non-planar graphs (Ko and Park 2012), which shows that the
elimination machinery does detect torsion when it exists. Run as a script it prints UConf_n of
K_4, K_5, K_{3,3} for n = 2, 3, 4.

### `graphs.py`

Enumerates connected multigraphs with all degrees at least 3, by genus, up to isomorphism.
These are exactly the homeomorphism classes of connected graphs without leaves or vertices
of valence 2. It was used to produce `data/survey/multigraphs.json`, the list of all 113
classes of genus 3 and 4 with at least three essential vertices.

### `explode.py`

Rational bookkeeping for the ordered vertex explosion (the complex of Γ is the mapping cone of a
map Δ between complexes of the graph Γ_v in which v is exploded; see `validation/check_explosion.py`
and `analysis/vexp.cpp`). From the Betti
numbers of Conf_n Γ, Conf_n Γ_v and Conf_{n−1} Γ_v it computes the ranks of the connecting maps
of the long exact sequence.

### `toric.py`

Computes the "toric bound", the largest degree that products of disjoint basic classes can
reach (a basic class is a particle moving around a cycle, or two particles exchanging at a
vertex of valence at least 3). `python3 toric.py` compares it with the top nonzero degree for every
run of the exact survey and lists the exceptions (16 runs, 14 distinct cases, all one degree above the
bound; none below it). Needs `networkx`.

### `fisharp.py`

The binomial inversion test. If the homology were a sum of induced representations of FI#
type (Church, Ellenberg and Farb 2015), the numbers w_k = Σ_j (−1)^{k−j} C(k,j) b_i(Conf_j)
would be non-negative integers. The script prints them for several graphs.

## 7. The analysis tools `analysis/`

These are standalone C++ programs that ask finer questions than the engine. They share a
block of helper code (graph structure, cell encoding, enumeration, boundary, an F_p
eliminator) which is documented line by line in `prod2.cpp`; the other files say so in their
header. **Their graph input has no particle number on the first line**: it is `V E`
followed by the edges, and n is a command-line argument.

| program | usage | what it computes |
|---|---|---|
| `prod2` | `prod2 n prime [k] < graph` | For a graph with exactly two essential vertices u, w, whether H_2(Conf_n) (the top homology) is spanned by star×star tori, the products of a star class at u and a star class at w. Prints the dimension of H_2 mod the prime and the rank of the span of the tori. |
| `g2tors` | `g2tors n x prime < graph` | The "triangle shortcut". For a graph Γ such that deleting vertex x leaves two essential vertices, it builds the relation lattice M from the vertex explosion and prints its rank mod the prime. The p-torsion of H_2(Conf_n Γ) has rank_Q(M) − rank_p(M) cyclic summands. |
| `g2split` | `g2split n x p [samples] < graph` | Which split types of tori (k particles at u, n − k at w) produce the torsion found by `g2tors`. |
| `phi` | `phi n prime < graph` | Rank of an auxiliary "double push" map C_3 → C_1 used in the analysis of the free part of H_2 for graphs with three essential vertices. |
| `decomp` | `decomp n prime [mode] < graph` | Dimension of the span of all tori (a star class at an essential vertex x times a 1-cycle of Γ minus x) inside H_2, mod a prime. With the environment variable `DUMP=file` (and a large prime) it also writes every torus as an integral chain, for use as relations in `certificates/orderq`. Mode 1 keeps only the basis cycles supported at one vertex; it is only a lower bound for the star×star span (the exact star×star span is Σ_x b_2(Conf_n Γ_x) − (rank ∂_3 − rank Φ), from `prod2` and `phi`; see the header of `decomp.cpp`). |
| `cutimg3` | `cutimg3 n prime mode [edges ...] < graph` | Adds to the tori the images of H_2 of graphs obtained by cutting given edges. Used to show that the classes of T(2,2,2) not spanned by tori are supported on Θ_3 subgraphs. |
| `vexp` | `vexp n v k p K < graph` | p-primary torsion of H_k(Conf_n Γ) through the ordered vertex explosion at vertex v. It checks its own hypothesis (no p-torsion in H_{k−1}(Conf_{n−1} Γ_v)) and exits with code 3 if it fails. |

Example, with T(2,2,2) as the graph (vertex 2 is x):

    printf "3 6\n0 1\n0 1\n0 2\n0 2\n1 2\n1 2\n" | analysis/g2tors 6 2 2         # rankM = 71548
    printf "3 6\n0 1\n0 1\n0 2\n0 2\n1 2\n1 2\n" | analysis/g2tors 6 2 1000003   # rankM = 71550
    # the difference 2 is the number of Z/2 summands in H_2(Conf_6 T(2,2,2))

## 8. Certificates `certificates/`

### 8.1 What a certificate proves

A torsion certificate is a short, independently checkable proof that a homology group has an
element of a given order P, where P is a prime p or a prime power p^k. It consists of

* an integral chain w of degree i+1 with ∂w ≡ 0 mod P, so that z = ∂w/P is an integral cycle
  and P·z is a boundary, and
* a cocycle f of degree i with values in Z/P such that f(z) is prime to p.

If (P/p)·z were a boundary ∂u, then (P/p)·f(z) = f(∂u) = 0 mod P, so p would divide f(z).
Hence [z] has order exactly P. For P = p this is the statement that z is not a boundary.

A certificate may also contain the line `SUBGRAPH h`. Then the verifier also checks that z
lives on the subgraph H formed by the first h edges (no particle on a later edge, and at every
occupied vertex both the generator half-edge and the base half-edge lie in H). Since the
edges of H come first, the base half-edges of H lie in H and these cells form the reduced
complex of H. So z is a cycle of Conf_n(H), and the certificate proves that the image of its
class in H_i(Conf_n Γ) has order exactly P.

### 8.2 Programs

| file | what it does |
|---|---|
| `verify_cert.py` | The independent verifier (a copy of `validation/verify_cert.py`). It shares no code with the engine or the finders. It rebuilds the complex from the definition, checks that every cell named in the certificate exists, recomputes ∂w, checks divisibility by P and ∂z = 0, checks the subgraph condition if one is given, enumerates every cell of degree i+1 (comparing their number with the closed formula) and checks that f vanishes mod P on all their boundaries, and checks that f(z) is prime to p. It rejects a corrupted certificate. Coefficients may be negative and arbitrarily large. |
| `cert.cpp` | Finder for mod-2 certificates. Usage `cert deg seed outfile < graph_with_n` (the input has the engine's format, with n on the first line). |
| `certp.cpp` | Finder for mod-p certificates, any prime. Usage `certp deg prime seed outfile < graph_with_n`. It reduces the maps with unit pivots, draws random kernel vectors and stops at the first pair (w, f) with nonzero pairing. |
| `certq.cpp` | Finder for certificates of prime-power order p^k, and generator export. Usage `certq deg p k index outfile [K] < graph_with_n`. It eliminates ∂_{i+1} mod p^(k+1) with unit pivots while logging every row operation, takes the local Smith form of the residual block with both transformations, and turns a diagonal entry p^k·unit back into a chain w and a cocycle f. With index −1 it writes integral cycles generating every cyclic summand of order below p^K, each with a chain w_t (∂w_t = order · z_t) and, for the summands of order p^k, a detecting cocycle. |
| `orderq.cpp` | Where does a cycle live? Coordinates of given cycles in the torsion of M = C_i / (B_i + C_i(H) + relations), localised at p, where H is Γ minus some edges and the relations are optional extra cycles. Usage `orderq deg p K cycles deleted_edges base_overrides [relations] < graph_with_n`. A torsion cycle is in the image of H_i(Conf_n H) (plus the relations) exactly when its coordinates vanish. The header explains how to make the conclusion exact (compare with the rank over a large prime). |
| `exactsolve.cpp` | Writes a cycle x as y + ∂u exactly over Z, with y supported on the subgraph of the first h edges. Usage `exactsolve deg h target out_prefix < graph_with_n`. Pivots ±1 only, a dropping step that keeps the residual small, and `solve_residual.py` (Hermite normal form) for the residual block. |
| `solve_residual.py` | The exact integral solver for the residual block of `exactsolve`. |
| `reducemod.cpp` | Eliminates the boundaries of Γ once (mod p^K) and writes many cycles as short reduced vectors, so that span questions become small. Usage `reducemod deg p K blocks out < graph_with_n`. |
| `spanq.cpp` | Answers queries "which torsion lies in the span of these cycles modulo boundaries?" on the output of `reducemod`, one per line of stdin. Its answers are mod p^K and must be confirmed exactly (see its header and `data/survey/k4_family/support/tori/exact_check.py`). |
| `reducecyc.cpp` | Reduces an integral cycle modulo boundaries with pivots ±1, giving its coordinates in a basis of H_i when the maps reduce completely. Usage `reducecyc deg cycle out < graph_with_n`. |
| `swiat_cells.h`, `elim_modq.h`, `elim_z.h` | Code shared by the four programs above. They contain the cells of the complex (encoding, enumeration, boundary, text form), sparse elimination mod p^K with an operation log and a dense local Smith form, and exact elimination over Z with pivots ±1. Each header documents its conventions and the facts the programs rely on. |

### 8.3 Certificate files

| file | what it proves | verification |
|---|---|---|
| `cert_T222_n6.txt.gz` | an element of order 2 in H_2(Conf_6 T(2,2,2)) | 1 min |
| `cert_T123_n6.txt.gz` | an element of order 2 in H_2(Conf_6 T(1,2,3)) | 1 min |
| `cert_T133_n6_p3.txt.gz` | an element of order 3 in H_2(Conf_6 T(1,3,3)) | 1 min |
| `cert_doubledK4_n5_p2.txt.gz` | an element of order 2 in H_2(Conf_5 of the doubled K_4) | 2.5 min, 3 GB |
| `cert_tripledC4_n5_p3.txt.gz` | an element of order 3 in H_2(Conf_5 of the tripled C_4) | 3 min, 3 GB |
| `cert_doubledK4_n5_order4.txt.gz` | an element of order exactly 4 in H_2(Conf_5 of the doubled K_4) | 1 min |
| `cert_K4dbl5_n5_order8.txt.gz` | an element of order exactly 8 in H_2(Conf_5) of K_4 with 01, 02, 03, 12, 13 doubled | 1 min |
| `cert_K4dbl5_n5_order16.txt.gz` | an element of order exactly 16 in the same group | 1 min |
| `cert_two_tori_into_K4dbl5_order16.txt.gz` | the class of order 16 as a SUM OF TWO TORI on the triangle with a leaf: a 653-cell cycle with coefficients ±1, ±2 whose image has order exactly 16 (`SUBGRAPH 6`) | 1 min |
| `cert_triangle_leaf_into_K4dbl5_order16.txt.gz` | a cycle on the triangle with two doubled edges and a leaf (the first six edges) whose image in the same group has order exactly 16 (`SUBGRAPH 6`) | 1 min |

To check a certificate: `python3 certificates/verify_cert.py certificates/cert_T222_n6.txt.gz`.
It prints each check (a) to (e) as it passes and a final line beginning with `VERIFIED`.

To make new ones:

    printf "3 7 6\n0 1\n0 2\n0 2\n0 2\n1 2\n1 2\n1 2\n" | certificates/certp 2 3 5 my_cert.txt
    G=data/survey/k4_family/support/graphs
    certificates/certq 2 2 4 0 cert16.txt < $G/K4dbl5_n5.txt        # order 16, about 4 min and 3 GB

### 8.4 Explicit generators and subgraphs, step by step

This is how the cycle on the triangle with a leaf was obtained (all files are in
`data/survey/k4_family/support/`, see its README for the results).

1. `certq 2 2 4 -1 gens.txt < graphs/K4dbl5_Hfirst_n5.txt` writes twelve generators of the
   torsion (Z/2)^7 ⊕ (Z/8)^4 ⊕ Z/16, in an edge order that lists the triangle with a leaf first.
2. `orderq 2 2 10 gens.txt 6,7,8,9,10 - < graphs/K4dbl5_Hfirst_n5.txt > img.json` gives their
   coordinates modulo the image of the subgraph (edges 6 to 10 deleted), and `kernel.py img.json`
   finds the torsion that comes from it, (Z/2)^4 ⊕ (Z/8)^4 ⊕ Z/16, and an element of order 16 in it
   (z_2 + z_11).
3. `exactsolve 2 6 x.txt sol < graphs/K4dbl5_Hfirst_n5.txt` writes that element as y + ∂u with
   y on the subgraph, exactly over Z.
4. The certificate is W = w_2·8 + w_11 − 16u, the cocycle of z_11, `PRIME 16` and `SUBGRAPH 6`.

## 9. Validation `validation/`

Each script checks one part of the software against something independent. All of them
print one line per test and a final verdict.

| script | what it checks | time |
|---|---|---|
| `check_explosion.py` | The ordered vertex explosion theorem at the level of chains: for several graphs and vertices, the complex of Γ is compared cell by cell and sign by sign with the mapping cone of Δ. Uses its own implementation of the complex. Prints `ALL OK`. | 1 min |
| `check_quotients.py` | The engine's quotient complexes. For 22 combinations of graph, group and character it builds the complex independently, applies each group element as a signed permutation, checks that every group element commutes with the boundary, and computes the homology of the invariant subcomplex with FLINT. Prints `ALL OK`. | 10 min |
| `verify_cert.py` | The certificate verifier (Section 8). | 1 min |
| `abrams_vs_swiatkowski.py` | The engine against Abrams' discretised model, an entirely different cell complex for the same space, on small cases including the triangle graphs at n = 3, 4. Every line must end with `agree=True`. | a few min |
| `elimination_vs_flint.py` | The engine's exact sparse elimination against FLINT's Smith form on 300 random integer matrices. Prints `mismatches = 0`. | 1 min |
| `local_smith_suite.py` | The local mode `-q` on every known torsion example, against the known answers. `--quick` skips the two half-hour T(2,2,3) runs. | 10 min (quick) |
| `dense_phase_suite.py` | Every dense finishing phase (byte, Z/4, Z/9) against pure sparse elimination, forcing early switches. All lines must say MATCH. | 3 min |
| `vexp_validate.py` | `analysis/vexp` against the direct computation, on T(2,2,2), T(1,3,3), T(1,2,2), K_4, K_5, a theta graph with leaves, the tripled C_4 and the doubled K_4. Needs about 4 GB. | about 1 hour |
| `check_generator_tools.py` | Every computation behind the explicit generators (Section 8): all certificates, the generator export, the subgraph image test, exactsolve and reducecyc against the saved results, the torus span, a rebuilt two-tori certificate, and the fast span path (reducemod + spanq) against the independent slow path (orderq) on random sets of tori. Steps can be run separately (`certs gens exact tori span`). | 2 hours, 6 GB |
| `regression_vs_reference.py` | Usage `regression_vs_reference.py NEW REF`. Compares two engine binaries on 33 small tests that cover every mode; used to confirm that the commented release engine behaves exactly like version 26 (`history/engine_versions/oconf26.cpp`). | 2 min |

The outputs of these scripts from the final run of the package are in `validation/logs/`.

## 10. Surveys and campaign scripts `surveys/`

These scripts produced the data in `data/`. Each begins with a docstring that says what it
does, how to call it and where its results are. `_common.py` holds shared helpers (graph
constructors, a wrapper that runs the engine and records its peak memory). Run them from
inside the `surveys/` folder.

| script | what it does |
|---|---|
| `make_jobs.py` | Builds the job list of the torsion survey (named graphs, all genus 3 and 4 classes, sun graphs, graphs with leaves), sorted by complex size. |
| `search.py` | The exact integral survey: full homology over Z of every job up to a cell limit. Restartable. Produced `data/survey/survey_exact.jsonl` (556 runs: 492 distinct pairs (Γ, n) on 150 graphs up to isomorphism). |
| `modp_survey.py` | Torsion detection through field dimensions and the universal coefficient theorem. |
| `local_survey.py` | Exact p-primary torsion through the local Smith form, plus rational Betti numbers. |
| `triangle_family.py` | Exact homology of the triangle multigraphs T(a,b,c). This is where the torsion was first found. |
| `triangle_shortcut.py` | Torsion of H_2(Conf_n T(a,b,c)) through the vertex explosion shortcut (`analysis/prod2` and `analysis/g2tors`), after checking its hypotheses. |
| `verify_case.py` | Re-checks a torsion case in different ways (exact, over F_p, after a random relabelling, on homeomorphic graphs with different complexes). |
| `colored_quotients.py` | Young subgroup quotients of K_4, K_5, K_{3,3}, K_6 (positive controls). |
| `colored_T.py` | Young subgroup quotients of Conf_6 of the triangle graphs. |
| `unordered_scan.py` | UConf_n of standard graphs with `python/unordered.py`. |
| `headline.py` | Single large exact computations (K_5, K_{3,3}, ...). |
| `nonplanar_n4.py` | The first non-planar tests at n = 4. |
| `k5_partial_doublings.py` | K_5 with some edges doubled at n = 5 (34 classes up to symmetry). |
| `nonplanar_campaign.py` | Lists every run of the non-planar campaign (`--list`), re-runs any of them with the current engine (`--run NAME`, `--run-small`) and compares with the recorded answer, and checks a local run against ranks over a large prime (`--rational NAME`), which rules out summands of order ≥ p^K that a run mod p^K cannot see (README Section 15). |

## 11. Data `data/`

| path | contents |
|---|---|
| `survey/survey_exact.jsonl` | 556 exact integral computations (492 distinct up to isomorphism of the graph) (one JSON object per line with the graph, n, cells, Betti numbers, torsion, Euler characteristic and Gal's value) |
| `survey/multigraphs.json` | all multigraphs with degrees ≥ 3 of genus 2 to 4 |
| `survey/jobs.json` | the job list of the survey |
| `survey/family.jsonl`, `family2.jsonl` | the triangle family T(a,b,c) |
| `survey/triangle_family_shortcut.jsonl` | the shortcut computations (`triangle_shortcut.py`) |
| `survey/triangle_local_snf.jsonl`, `local.jsonl`, `modp.jsonl` | T(2,2,3) by local Smith form and by field ranks |
| `survey/colored.jsonl`, `colored_T.jsonl` | Young quotients |
| `survey/unordered_scan.jsonl` | unordered configuration spaces |
| `nonplanar/campaign.json` | every run of the non-planar campaign with options and results (77 records). The `dims` field gives the dimensions of the homology of the quotient complex over the residue field and `torsion` gives the local stage counts by degree. |
| `nonplanar/raw/` | the engine's JSON output and log of every campaign run, and `run_status.txt` with return codes and times |
| `nonplanar/earlier/` | logs of the earlier experiments (partial doublings of K_5, fat graphs, group quotients of the tripled C_4, benchmarks, elimination diagnostics, vertex explosion runs) |
| `survey/k4_family/` | K_4 with some edges doubled at n = 5: exact homology of all doubling patterns, Young and odd-order symmetry quotients, independent checks of the Z/8 and Z/16 (its own README) |
| `survey/k4_family/support/` | explicit generators of the torsion of K_4 with five doubled edges, the subgraph search for the class of order 16, the cycle on the triangle with a leaf, the torus decomposition test, and the graph files to rerun all of it (its own README) |
| `accounting/` | a single accounting of every recorded computation: `accounting.csv` (one row per computation, 804 rows), `accounting_graphs.csv` (one row per graph up to isomorphism, 197 graphs), `accounting_summary.txt`, `accounting.xlsx` (the same as a workbook, with a summary sheet and the October 2026 re-verification), and the scripts `make_accounting.py` and `make_xlsx.py` that build them from the files above |

## 12. Example inputs `examples/graphs/`

Engine input files (`V E n` on the first line) for every graph used in the campaign. The
names follow the campaign: `t222_6.in` is T(2,2,2) with six particles, `t133a_6.in` is
T(1,3,3) with an apex vertex, `dk4_5.in` the doubled K_4 at n = 5, `tc4_5.in` the tripled
C_4, `k5_6.in` K_5, `k33_7.in` K_{3,3}, `dk5_5.in` the doubled K_5, `k5T133_6.in` K_5 with a
T(1,3,3) inside (some edges multiplied), and so on. `campaign.json` lists the exact edge
list of every run.

## 13. History `history/`

`history/HISTORY.md` describes the 26 engine source files (all kept in
`history/engine_versions/`) and the original driver scripts (`history/scripts/`) with their
logs. They are not needed to use the package.

## 14. Reproducing the main results

| result | command | time and memory |
|---|---|---|
| H_*(Conf_6 T(2,2,2)) | `python3 -c "import sys; sys.path.insert(0,'python'); from gconf import homology; print(homology([(0,1),(0,1),(0,2),(0,2),(1,2),(1,2)], 6))"` | 2 min, 1 GB |
| H_*(Conf_6 T(1,2,3)) | same with edges `[(0,1),(0,2),(0,2),(1,2),(1,2),(1,2)]` | 2 min, 1 GB |
| H_*(Conf_6 T(1,3,3)) | same with edges `[(0,1)]+[(0,2)]*3+[(1,2)]*3` | 8 min, 4 GB |
| T(2,2,3), n = 6 | `engine/oconf -q 2 3 -D 70000 < examples/graphs/t223_6.in`, and the same with `-q 3 2` | 30 min each, 3 GB |
| doubled K_4, n = 5 | `engine/oconf -q 2 3 < examples/graphs/dk4_5.in` | 5 min, 1.5 GB |
| tripled C_4, n = 5 | `engine/oconf -q 3 2 < examples/graphs/tc4_5.in` | 2 min, 1.3 GB |
| the certificates | `python3 certificates/verify_cert.py certificates/cert_T133_n6_p3.txt.gz` | 1 min |
| order 16, and its support | `python3 certificates/verify_cert.py certificates/cert_triangle_leaf_into_K4dbl5_order16.txt.gz` | 1 min, 3 GB |
| generators and subgraph tests | Section 8.4 | 5 to 10 min per step, 3 to 6 GB |
| the shortcut | `cd surveys && python3 triangle_shortcut.py '[[1,2,4,6]]' out.jsonl` | 2 min |
| any campaign run | `cd surveys && python3 nonplanar_campaign.py --run t133a_n6_l` | see `campaign.json` |

## 15. Limits, pitfalls and practical advice

* **Memory is the binding constraint.** The engine needs a few hundred bytes per cell of
  the complex (T(1,3,3) at n = 6, about 13 million cells, peaks at 3 to 4 GB). On an 8 GB machine, complexes up to about 15 million cells (or quotient
  complexes up to about 25 million orbit cells in local mode) are feasible.
* **Exact mode can stall.** Exact elimination over Z can suffer from entry growth in the
  residual block (for T(2,2,3) it stops after about 40 minutes with exit code 3, a 64-bit overflow). Use `-q p K` for the p-primary
  part together with a run over a large prime for the rational Betti numbers.
* **Local mode sees only orders below p^K.** A run with `-q p K` counts exactly the cyclic
  summands of order p, p², ..., p^(K−1). Elementary divisors divisible by p^K are 0 mod p^K and
  are invisible; they are not reported anywhere in the output. To rule them out, compare the
  total number of pivots of each map (the sum of its `stages`) with its rank over a large prime
  (`-p 1000003`): the two agree exactly when no elementary divisor is divisible by p^K. In
  particular, a run with K = 2 (Z/4, Z/9) that finds nothing proves only that there are no
  summands of order exactly p.
* **Quotients see only one isotypic part** (Section 5). A torsion-free quotient does not prove
  that the full space is torsion-free.
* **Never use a group whose order is divisible by p.** The engine refuses (exit code 10). An
  early run without this guard used a group of order 24 at p = 3 and produced fake 3-torsion.
* **Choose invariant base half-edges.** With automorphisms the engine does this itself, and
  it reports its choice.
* **Use `-F` with every group.** It gives the same complex and is much faster.
* **Watch out for `pkill -f`.** Killing engine runs with `pkill -f oconf` can also kill the
  shell that started them, if its command line contains the same word. Kill by process id.

## 16. References

The works the code relies on directly are

* B. Hainaut, B. Knudsen and A. Wawrykow, *Representation asymptotics in the homology of pure
  graph braid groups*, arXiv:2510.00201 (2025). Definition 3.12 and Theorem 3.14 give the
  reduced ordered Świątkowski complex; Conjecture 1.7 is the folklore conjecture of
  torsion-freeness.
* J. Świątkowski, *Estimates for homological dimension of configuration spaces of graphs*,
  Colloquium Mathematicum 89 (2001), 69–79.
* A. Abrams, *Configuration spaces and braid groups of graphs*, PhD thesis, University of
  California, Berkeley (2000).
* Ś. R. Gal, *Euler characteristic of the configuration space of a complex*, Colloquium
  Mathematicum 89 (2001).
* K. H. Ko and H. W. Park, *Characteristics of graph braid groups*, Discrete and Computational
  Geometry 48 (2012).
* B. H. An, G. C. Drummond-Cole and B. Knudsen, *Subdivisional spaces and graph braid groups*,
  Documenta Mathematica 24 (2019).
* The FLINT team, *FLINT: Fast Library for Number Theory*, https://flintlib.org, and its
  Python bindings python-flint.
