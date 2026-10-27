# graph-braid-torsion

Software for computing the integral homology of unordered configuration spaces of graphs,
with particular attention to torsion, together with every script and data file from a study
of which finite cyclic groups occur in these homology groups, and in which graphs.

The short version of what was found is the following. For a finite simple graph Γ write
B_k Γ = UConf_k(Γ) for the unordered configuration space of k points on Γ. In degree one the
torsion of H_1(B_k Γ; Z) is classical (Ko and Park): for k ≥ 2 it is a single Z/2 when Γ is
not planar, and zero otherwise. Above degree one the computations here show that

* **torsion other than Z/2 occurs.** The smallest example found is

      H_2(B_3 K_{4,5}; Z) = Z^155 ⊕ (Z/2)^4 ⊕ Z/3

  and the first non-squarefree example is

      tors H_3(B_4 K_{5,7}; Z) = (Z/2)^14 ⊕ (Z/3)^6 ⊕ Z/4.

  Both come from the same source. For the complete bipartite graph K_{n+1,b} with b ≥ n, the
  torsion of H_{n−1}(B_n K_{n+1,b}) is the torsion of the cokernel of the 0/1 inclusion matrix
  of (n−1)-subsets against n-subsets of a b-set, so Wilson's diagonal form for that matrix gives
  it in closed form. The scripts check this identification on every case they can build.
* **Z/4 occurs in H_2, without any large bipartite minor.** K_{4,5} with a triangle added on
  its 5-vertex side (9 vertices, 23 edges) has

      tors H_2(B_3 Γ; Z) = (Z/2)^11 ⊕ Z/3 ⊕ (Z/4)^2,

  confirmed by two independent computations. Wilson's mechanism alone gives nothing above Z/3
  in H_2. A survey of all 313 supergraphs of K_{4,5} on its nine vertices finds Z/4 in exactly
  12 of them: those where the edges added on the 5-vertex side form a triangle, a triangle with
  a pendant edge, K_4 minus an edge, or K_4 (so one vertex of that side stays untouched), and
  the edges added on the 4-vertex side are none, one, or two disjoint edges. A star K_{1,4} on the 5-vertex side gives (Z/3)^3 instead. No summand of
  order 5 or more, and no Z/8, has been found anywhere.
* **planar graphs carry torsion in H_2, and so far only Z/2.** The octahedron, the dodecahedron
  and the icosahedron all have exactly Z/2 in H_2(B_4). A search over 5,714 planar graphs,
  including every planar graph of minimum degree 3 on 8 and on 9 vertices, found Z/2 and
  nothing else, and no planar graph in the data has more than a single Z/2 in H_2.
* **all odd torsion found lives on a K_{4,5} minor.** Of the 207 graphs with a summand other
  than Z/2, 206 contain K_{4,5} as a subgraph and the remaining one contains it as a minor.
* **complete graphs carry only Z/2.** For 8 ≤ n ≤ 12, H_2(B_3 K_n) has torsion (Z/2)^m with
  m = C(n−1, 2) = rank H_1.
* **the smallest carriers are known.** Among all 853 connected graphs on seven vertices, K_7 is
  the only one with torsion in H_2(B_3). At k = 4 the octahedron is the carrier of H_2 torsion
  with the fewest vertices, and no graph on at most seven vertices with fewer than 12 edges has
  any torsion in H_2(B_4). The octahedron is not the only minimal carrier, though: three
  seven-vertex graphs with 12 edges, and a seven-vertex graph with a pendant edge are all
  minimal under edge deletion, and none of them contains the octahedron.

All told the stored data hold 10,653 computations on 8,140 graphs up to isomorphism (10,062
distinct pairs (Γ, k), with k between 2 and 12). Section 12 explains how that accounting is
made and how to regenerate it.

This README is the reference manual: it explains how to install the package, how to use the
Python library, and what every file in the package does. **If you are new to the package, start
with [TUTORIAL.md](TUTORIAL.md)**, a sequence of short lessons for a mathematician who is
comfortable at a terminal but is not a programmer; every command in it works as written. The
README is long because it is meant to be usable by someone who has never seen the code.
Sections 1 to 3 are enough to start computing. Section 4 is the
reference for the library, and Section 5 explains how to choose between its three pipelines.
Sections 6 to 10 go through the experiment scripts one file at a time, Sections 11 and 12
describe the stored data, Sections 13 and 14 are for reproducing the results and for
avoiding the mistakes that are easy to make, and Section 15 describes the example scripts, the
Makefile and the logs of our own validation runs.

## Contents

1. Installation
2. A first computation
3. The mathematical model in brief
4. The library `gbt/`: every module and function
5. Choosing a pipeline, and how the pipelines work
6. The experiment scripts: conventions shared by all of them
7. Experiments on K_{n+1,b}: the inclusion-matrix mechanism
8. Sweeps over families of graphs
9. The structure of the torsion classes
10. Verification and regression
11. Data `data/`
12. The inventory of all computations
13. Reproducing the main results
14. Limits, pitfalls and practical advice
15. Example scripts `examples/`, the Makefile, and the logs `logs/`
16. References

---

## 1. Installation

The package needs

* Python 3.8 or later. It has been run with Python 3.11. The library `gbt/` uses only the
  standard library: there is nothing to compile, and all arithmetic is in Python integers, so
  there is no overflow;
* the Python package `networkx` (version 2.6 or later), but only for the scripts that
  enumerate graphs or test planarity, isomorphism or subgraph containment. Install it with
  `pip install networkx` (on some systems `pip install --user networkx` or
  `pip install --break-system-packages networkx`). The scripts that need it are listed in
  `requirements.txt`: `sweep_small_graphs`, `minimal_carriers`, `petersen_family`,
  `linkless_controls`, `planar_odd_search`, `genuine_higher_torsion`, `higher_degree_torsion`,
  `prime_census`, `theta_class`, `inventory`, `supergraph_walks` and `torsion_hunt`, together
  with the two library functions `gbt.graphs.from_networkx` and `gbt.graphs.petersen_family`;
* only for `torsion_hunt.py`: `numpy` (`pip install numpy`) and the graph generator `nauty`
  (its programs `nauty-geng` and `nauty-planarg`; `sudo apt-get install nauty` on Debian or
  Ubuntu, `brew install nauty` on macOS, where the programs are called `geng` and `planarg` and
  the two names in the script must be changed accordingly).

There is no installation step. Unpack the folder and run the scripts from anywhere:

    cd graph-braid-torsion
    pip install -r requirements.txt       # only networkx
    python3 experiments/validate.py       # 60 regression checks; prints ALL CHECKS PASSED
    make check                            # the whole fast validation suite, about 5 minutes

`make` is optional: the `Makefile` only collects commands (Section 15), and each of its targets
is a few `python3 ...` lines that can be typed by hand.

Every script in `experiments/` puts the repository root on `sys.path` itself, so
`python3 experiments/NAME.py` works from any working directory. To use the library from your
own code, either run Python from the repository root or add that folder to `PYTHONPATH`.

The scripts in `experiments/` write their results into `data/`, which is located relative to
the package (`gbt/paths.py`), not relative to the current folder. They never write anywhere
else. The example scripts in `examples/` write nothing, except `sweep_template.py`, which
writes to `examples/output/`.

## 2. A first computation

From the repository root:

    >>> from gbt import complete, complete_bip, octahedron, fast, fmt, primary
    >>> H, raw, reduced = fast.compute(complete(5), 2)     # K_5, two particles
    >>> H
    {0: (1, []), 1: (6, [2]), 2: (0, [])}

The answer `H` is a dictionary from degree d to a pair (free rank, list of invariant factors
larger than 1). So H_0 = Z, H_1 = Z^6 ⊕ Z/2 and H_2 = 0, which is Abrams' computation that
B_2 K_5 is the closed non-orientable surface of genus 7. The other two return values are the
number of cells in each degree before and after the discrete Morse reduction:

    >>> raw, reduced
    ({0: 55, 1: 150, 2: 90}, {0: 1, 1: 7, 2: 1})

`fmt` prints a group, and `primary` turns invariant factors into a primary decomposition,
which is usually what you want, because the question is which primes occur:

    >>> H, _, _ = fast.compute(complete_bip(4, 5), 3)      # K_{4,5}, three particles
    >>> fmt(*H[2])
    'Z^155 + (Z/2)^3 + Z/6'
    >>> primary(H[2][1])
    '(Z/2)^4 + Z/3'

    >>> H, _, _ = fast.compute(octahedron(), 4)            # planar, 26,142 cells, under a second
    >>> {d: fmt(*v) for d, v in H.items()}
    {0: 'Z', 1: 'Z^8', 2: 'Z^95 + Z/2', 3: 'Z^10', 4: '0'}

To get the torsion of a single degree of a large complex without computing everything else,
use the slice pipeline (Section 5):

    >>> from gbt import lean
    >>> primary(lean.torsion_in_degree(complete_bip(5, 7), 4, 3))   # about 15 s, 0.5 GB
    '(Z/2)^14 + (Z/3)^6 + Z/4'

Before building a complex you can ask how big it will be. These are closed formulas, so they
cost nothing:

    >>> from gbt import reduced_cell_count, cell_count
    >>> reduced_cell_count(complete(7), 3), cell_count(complete(7), 3)
    (25256, 46704)

## 3. The mathematical model in brief

This section gives the definitions needed to use the software.

**Configuration spaces.** For a finite graph Γ, Conf_k(Γ) is the space of k-tuples of pairwise
distinct points of Γ and B_k Γ = UConf_k(Γ) = Conf_k(Γ)/S_k is the unordered version. It is
aspherical, and its fundamental group is the graph braid group (Abrams 2000; Ghrist 2001).
Only simple graphs (no loops, no multiple edges) are handled by this package.

**The Świątkowski complex.** Świątkowski (2001) built a small chain complex computing
H_*(B_k Γ), and An, Drummond-Cole and Knudsen (2019) recast it as

    S(Γ) = Z[E] ⊗ ⊗_{v ∈ V} S_v,

bigraded by (homological degree, weight). In words:

* Z[E] is the polynomial ring on the edges. The exponent of an edge e counts the particles
  *parked* in the interior of e. Parked particles have degree 0 and weight 1 each.
* Each vertex module S_v is free abelian on three kinds of generator:
  `∅_v` (nothing at v; degree 0, weight 0), `v` (a particle *sits* at v; degree 0, weight 1),
  and `h_{v,e}` for each edge e at v (a particle at v *moves* out along e; degree 1, weight 1).
* The differential is ∂h_{v,e} = e·∅_v − v, extended over the tensor product as a derivation
  with Koszul signs. The moving generators are always read in increasing order of vertex, so
  the i-th moving vertex (counting from 1) contributes the sign (−1)^{i−1}.
* The summand of weight k computes H_*(B_k Γ) for every graph, with no subdivision hypothesis.

This is what makes dense graphs with few vertices tractable. Abrams' discretised model needs
the graph subdivided until every cycle has length at least k+1, which is where it explodes.

**Reduction.** At a vertex v one may pick a basepoint edge e₀(v) (the code takes the edge of
least index), change basis to the differences h̃_{v,e} = h_{v,e} − h_{v,e₀(v)}, and cancel the
pair (h_{v,e₀(v)}, v), which is joined by a unit entry of ∂. Then ∂h̃_{v,e} = e − e₀(v): two
parked terms and no sitting term. Doing this at every vertex gives the *fully reduced*
complex, in which every occupied vertex is moving, so

    homological degree = number of occupied vertices.

The number of cells of degree j in the fully reduced complex of weight k is

    e_j(d_v − 1 : v ∈ V) · C(M + k − j − 1, k − j),

where d_v is the degree of v, e_j is the j-th elementary symmetric polynomial and M is the
number of edges. In the unreduced complex each occupied vertex has d_v + 1 states instead of
d_v − 1. For K_7 at k = 3 that is 25,256 reduced cells against 46,704 unreduced ones.
Because the reduction happens inside one tensor factor, it can be done at any subset of the
vertices. That is the `reduce_at` option of Section 4.

**What "torsion" means computationally.** H_d = ker ∂_d / im ∂_{d+1}. The free rank of H_d is
(number of d-cells) − rank ∂_d − rank ∂_{d+1}, and the torsion of H_d is read off the invariant
factors of ∂_{d+1}, the differential coming *into* degree d. In fact more is true:

    0 → H_d → coker ∂_{d+1} → im ∂_d → 0

is exact, and it splits because im ∂_d is a subgroup of a free group. So
tors coker ∂_{d+1} = tors H_d, and the torsion of one degree can be computed from a single
matrix without knowing ∂_d at all. `gbt.lean.torsion_in_degree` is built on this.

**Two facts used throughout.** Wilson (1990) gives a diagonal form for the inclusion matrix
W_{t,k}(m) of t-subsets against k-subsets of an m-set: for t ≤ min(k, m − k) it is equivalent
to a diagonal matrix with entries C(k−i, t−i), each repeated C(m,i) − C(m,i−1) times,
i = 0, …, t. Gottlieb (1966) shows that for t ≤ k, W_{t,k}(m) has full rank min(C(m,t), C(m,k)). For
K_{n+1,b} with b ≥ n, the top differential ∂_n of the unreduced complex is, in a suitable
basis, the block matrix [[I, 0], [0, W_{n−1,n}(b)], [0, 0]]. So

    H_n(B_n K_{n+1,b}) = ker W_{n−1,n}(b),      of rank max(0, C(b,n) − C(b,n−1)),
    tors H_{n−1}(B_n K_{n+1,b}) = tors coker W_{n−1,n}(b)
                                = ⊕_{i=0}^{M−2} (Z/(M−i))^{C(b,i) − C(b,i−1)},  M = min(n, b−n+1).

In particular H_n vanishes exactly when b ≤ 2n−1. Section 7 lists the scripts that check each
piece of this.

## 4. The library `gbt/`: every module and function

All modules have a docstring at the top explaining the mathematics and the conventions, and
every public function has its own docstring. `from gbt import *` brings in everything listed
in `gbt/__init__.py`. That covers all the functions below except those in `gbt.linalg`,
`gbt.report`, `gbt.paths` and `gbt.flipchains`, which are imported from their modules.

### `gbt.graphs`

A graph is a `Graph` object. Its vertices are 0, …, N−1, and the *index* of an edge in
`G.edges` is that edge's name everywhere else in the package. Both conventions matter: the
Świątkowski signs depend on the vertex order, and partial reduction depends on the choice of a
basepoint edge at each vertex.

| name | what it does |
|---|---|
| `Graph(N, edges)` | a finite simple graph. Attributes `N`, `M`, `edges` (a list of pairs (u, v) with u < v). Methods `deg(v)`, `nbrs(v)` (in the order of the incident edge indices), `other(v, ei)` (the other endpoint of edge `ei`), `edge_index(u, v)` (raises `KeyError` if absent), `is_connected()`, `delete_edge(ei)` (a copy with that edge removed, same vertex set) |
| `complete(n)` | K_n |
| `complete_bip(a, b)` | K_{a,b} with the a-side first: A = {0, …, a−1}, B = {a, …, a+b−1}, edges in lexicographic order of (a, b). The sign bookkeeping in `gbt.flipchains` assumes this order |
| `from_edges(edges)` | a graph from any edge list, with the vertices relabelled to 0, …, N−1 |
| `from_networkx(Gx)` | conversion from a networkx graph with any hashable labels (needs networkx) |
| `disjoint_union(G, H)` | G ⊔ H, with H's vertices shifted up by G.N |
| `subdivide(G, t)` | replaces every edge by a path of t edges (t = 1 changes nothing) |
| `tetrahedron()`, `octahedron()`, `dodecahedron()`, `icosahedron()` | the tetrahedron, octahedron, dodecahedron and icosahedron. `OCTAHEDRON_EDGES`, `DODECAHEDRON_EDGES`, `ICOSAHEDRON_EDGES` are the fixed edge lists, and the octahedron's antipodal pairs are {0,1}, {2,3}, {4,5} |
| `platonic_solids()` | these graphs as (name, Graph) pairs, in order of edge count |
| `prism_graph(n)`, `prism()` | the n-prism (two n-cycles joined by a matching), and the triangular prism |
| `antiprism(n)` | the n-antiprism (two n-cycles joined by a zigzag band; 4-regular) |
| `theta(paths=3, length=2)` | two hub vertices joined by `paths` internally disjoint paths of `length` edges |
| `wheel(n)` | the wheel W_n: an n-cycle plus a hub joined to every rim vertex |
| `petersen()` | the Petersen graph, as the Kneser graph on the 2-subsets of a 5-set |
| `petersen_family()` | the seven graphs of the Petersen family as (name, Graph) pairs, generated as the closure of K_6 under ΔY and YΔ moves, with isomorphism tests to remove duplicates (needs networkx) |

### `gbt.swiatkowski`

| name | what it does |
|---|---|
| `swiatkowski(G, k, reduce_at=())` | the Świątkowski complex of weight k as a `ChainComplex`, reduced at the vertices in `reduce_at`. `()` is the unreduced complex and `range(G.N)` the fully reduced one. Reducing at a vertex of degree 1 deletes all of its generators, so subsets containing a leaf are skipped |
| `edge_monomials(M, k)` | {r: all monomials of weight r in M edge variables}, as sorted tuples of edge indices with repetition |

A basis element is encoded as a pair `(mono, st)`. `mono` is a sorted tuple of edge indices
(the parked particles). `st` is a tuple of length N with `st[v] = 0` (v empty), `1` (a particle
sits at v) or `2 + e` (a particle at v moves along edge e). For example, for K_{3,3} at k = 2:

    >>> C = swiatkowski(complete_bip(3, 3), 2)
    >>> C.sizes()
    {0: 114, 1: 252, 2: 135}
    >>> c = ((0,), (0, 0, 0, 5, 0, 0))   # one particle parked on edge 0, one moving at 3 along edge 3
    >>> C.B[1][c]
    {((0, 3), (0, 0, 0, 0, 0, 0)): 1, ((0,), (0, 0, 0, 1, 0, 0)): -1}

which is ∂h_{3,e₃} = e₃·∅ − 3, multiplied by the parked particle on e₀.

### `gbt.chain`

The reference pipeline. It keeps the generators as readable labels, which makes it easy to
inspect, at the cost of a great deal of Python object overhead.

| name | what it does |
|---|---|
| `ChainComplex(gens, B)` | a chain complex of free abelian groups. `gens[d]` is a set of hashable labels; `B[d]` maps each generator of degree d to its boundary, a dict {generator of degree d−1: coefficient} |
| `.sizes()` | {degree: number of generators} |
| `.check_d2()` | True if ∂∘∂ = 0. A cheap and effective test of the sign conventions |
| `.reduce_fully()` | algebraic discrete Morse reduction (Forman): cancels every pair joined by a ±1 entry, top degree downwards, choosing pivots by least Markowitz fill (\|column\|−1)(\|row\|−1). Modifies the complex in place |
| `.dense(d)` | the matrix of ∂_d as (rows, number of rows, number of columns), with generators sorted by `repr` |
| `.homology(degrees=None)` | {degree: (free rank, sorted non-unit invariant factors)}. Call `reduce_fully()` first for anything but small complexes |
| `fmt(free, tors)` | a group as a string, e.g. `'Z^155 + (Z/2)^3 + Z/6'` |
| `primary(tors)` | invariant factors to primary decomposition, e.g. `[2,2,2,6]` → `'(Z/2)^4 + Z/3'` |
| `odd_part(x)` | the largest odd divisor of x |

### `gbt.smith`

| name | what it does |
|---|---|
| `smith(M, nr, nc)` | exact Smith normal form of an nr × nc integer matrix given as a list of rows. Returns (rank, invariant factors). Pivots of least absolute value, Euclidean reduction of row and column, and the extra step that makes the result a divisibility chain rather than merely diagonal. Meant for the small residual matrices left by the Morse reduction |

### `gbt.fast`

The same mathematics as `swiatkowski` + `chain`, for the fully reduced complex only, with
every generator replaced by an integer index into a per-degree list. About 100 times faster:
K_6 at k = 4 (71,220 cells) takes about two seconds here against several minutes through the
reference pipeline.

| name | what it does |
|---|---|
| `compute(G, k, verbose=False)` | build, Morse-reduce, Smith-normalise. Returns (H, raw sizes, reduced sizes) as in Section 2 |
| `build(G, k)` | the fully reduced complex: (sizes, B), where `B[d]` is a list of columns, each a dict {row index: coefficient} |
| `morse(sizes, B)` | unit-pivot elimination, top degree downwards. Returns the reduced (sizes, B) |
| `homology(sizes, B)` | {degree: (free rank, non-unit invariant factors)} of a (reduced) complex |

Every vertex must have degree at least 1. A graph with an isolated vertex has to go through
`swiatkowski(G, k, reduce_at=())`.

### `gbt.lean`

The same answers as `gbt.fast` in roughly a third of the memory. Cells are packed into single
Python integers (digits in base `BASE` = 2^11), and the matrices are built one degree at a
time from the top down. Because degree equals the number of occupied vertices, each degree
slice can be enumerated without the others, and only two index dictionaries are ever alive at
once.

| name | what it does |
|---|---|
| `compute(G, k, verbose=False)` | as `fast.compute` |
| `torsion_in_degree(G, k, d, verbose=False)` | the invariant factors of tors H_d(B_k G), from the two slices of degrees d and d+1 only, using tors H_d = tors coker ∂_{d+1}. The unit pivots are peeled with `linalg.peel_sweeps` and the rest goes to `linalg.sparse_smith`. This reaches much further than any full computation |
| `reduce_top_down(G, k, verbose=False, stop_degree=1)` | build and Morse-reduce one degree at a time. Returns (raw sizes, sizes, B) |
| `homology(sizes, B)` | as `fast.homology` |
| `slice_cells(G, k, d, mon)`, `boundary(G, k, d, mon, idx_low)` | the packed keys of the degree-d cells, and the columns of ∂_d in slice order |
| `encode(mov, mon, M)`, `decode(key, d, kd, M)` | (moving slots, parked edges) ↔ packed key |
| `BASE` | 2^11. Must exceed N·M for the graph in hand. The largest graph stored, K_13, has 13 · 78 = 1014 |

### `gbt.linalg`

Sparse integer and modular linear algebra. A matrix is a list of columns, each a dict
{row index: coefficient} with no zero entries. This is the format that all three builders
produce.

| name | what it does |
|---|---|
| `indexed_columns(C, d)` | the columns of ∂_d of a `ChainComplex` with integer indices. Returns (columns, {row generator: index}, column generators in order); generators are sorted by `repr`, so the indexing is the same on every run |
| `sparse_smith(cols, label='', verbose=True)` | exact (rank, non-unit invariant factors): unit-pivot peel, then `smith` on the remainder |
| `coker_torsion(cols, label='', verbose=True)` | the non-unit invariant factors of the cokernel, i.e. the torsion of Z^rows / image |
| `peel_sweeps(cols, nrows=None, verbose=False)` | the unit-pivot peel driven by repeated cheap sweeps rather than a priority queue. Much leaner on memory than a Markowitz heap, which holds one entry per nonzero and grows with every fill-in. Returns (residual columns, number of pivots) |
| `kernel_basis(cols, nrows)` | a Z-basis of the kernel lattice {x : Σ x_j col_j = 0} |
| `gf2_kernel(cols)` | a basis of the kernel mod 2, as integers whose bits are column indices |
| `rank_mod2(cols)` | rank over F_2 by bitset elimination |
| `rank_mod_p(cols, p)`, `rank_mod_p_sparse(cols, p, verbose=False)` | rank over F_p (odd p), the second by sparse Markowitz elimination |
| `p_primary_count(cols, p=2, big_prime=1000003, verbose=False)` | how many invariant factors are divisible by p: rank over a large prime minus rank over F_p |
| `betti_mod_p(C, p)` | {degree: dim H_d(C; F_p)} for a `ChainComplex` |
| `local_smith(cols, p, K, verbose=False)` | the local Smith form at p, computed modulo p^K: a list whose entry s is the number of invariant factors of p-adic valuation exactly s. Pivots only on units and divides by p between stages, so it never takes a gcd and never suffers entry growth; compare `sum(stages)` with the rank over a large prime to make sure nothing of order ≥ p^K is hidden. Written independently of `sparse_smith`, and the tool of choice when `sparse_smith` stalls (Section 14) |

### `gbt.discretised`

| name | what it does |
|---|---|
| `abrams(G, k)` | Abrams' discretised configuration space UD^k(G) as a `ChainComplex`: cells are k-sets of pairwise disjoint closed cells of G, and the degree is the number of edges used. It shares no code with the Świątkowski builders except the `Graph` class, which is the point: it is an independent check. It is homotopy equivalent to B_k G only when G is sufficiently subdivided (every path between vertices of degree ≠ 2 of length at least k−1, every cycle of length at least k+1; Abrams, Prue–Scrimshaw). **The function does not check this; the caller must subdivide first.** |

### `gbt.sizes`

| name | what it does |
|---|---|
| `reduced_cell_count(G, k)` | number of cells of the fully reduced complex of weight k (formula in Section 3) |
| `reduced_cell_count_by_degree(G, k)` | the same, as {degree: count} |
| `cell_count(G, k)` | number of cells of the unreduced complex |

### `gbt.flipchains`

Explicit chains in the unreduced complex of K_{n+1,b}, with A = {0, …, n} and
B = {n+1, …, n+b}. The flip chain of an edge e = {a, w} (a ∈ A, w ∈ B) is
z(e) = h_{w,e} − h_{a,e}, with ∂z(e) = a − w. Chains on disjoint vertex sets multiply.

| name | what it does |
|---|---|
| `z(n, b, a, w)` | the flip chain z({a, w}) as a factor |
| `X_B0(n, b, B0)` | X_{B₀} = Σ_g sgn(g) z(e₁)⋯z(e_n), over bijections g : B₀ ⊔ {∗} → A, for an n-subset B₀ of B |
| `Y_T(n, b, T)` | Y_T for an (n−1)-subset T of B. One has ∂X_{B₀} = (−1)^{n−1} Σ_{T ⊂ B₀} Y_T |
| `X(n, b)` | Σ_{B₀} X_{B₀}. Its boundary is ±(b−n+1)·y with y primitive, so (b−n+1)[y] = 0 in H_{n−1}. At b = 2n−1, where H_n = 0, the class [y] has order exactly n |
| `closed_formula_X(n, b)` | the coefficients of X from a closed formula, without expanding the products: the coefficient of a matching configuration (S, t) is (−1)^{\|S_A\|} sgn σ(S, t), where σ lists A as [S_A increasing \| t(v) for v ∈ S_B in increasing order of v \| the unused vertex] |
| `perm_sign`, `edge_index`, `word_to_cell`, `expand`, `add_into` | helpers: sign of a sequence, index of the edge {a, w}, an ordered product of generators turned into a (cell, sign), multiplying out factors, accumulating chains |

Chains are dicts {cell: coefficient} with cells in the encoding of `gbt.swiatkowski`.

### `gbt.report` and `gbt.paths`

| name | what it does |
|---|---|
| `report.homology_line(H, skip_degree_zero=True)` | `"H_1=Z^8  H_2=Z^95 + Z/2  H_3=Z^10  H_4=0"` |
| `report.torsion_line(H)`, `report.torsion_summary(H)` | the torsion only, as a line or as {degree: string} |
| `report.has_torsion(H, min_degree=0)`, `report.has_odd_torsion(H)`, `report.exponents(H)` | tests and summaries |
| `paths.DATA`, `paths.data_path(name)` | the data folder (`data/` next to `gbt/`), and a file inside it |
| `paths.load_json(name, default=None)`, `paths.save_json(name, obj)` | read and write a data file. `save_json` writes with sorted keys, so files diff cleanly |
| `paths.record(H)` | a homology dict made JSON-safe |

## 5. Choosing a pipeline, and how the pipelines work

| reduced cells | use | notes |
|---|---|---|
| up to about 5·10³ | `swiatkowski` + `ChainComplex` | the reference implementation; generators are readable labels, and it is the only one that supports an arbitrary `reduce_at` |
| up to about 6·10⁵ | `gbt.fast` | about 100× faster; fully reduced complex only |
| up to about 2·10⁶ | `gbt.lean.compute` | about a third of the memory of `fast` |
| torsion of one degree | `gbt.lean.torsion_in_degree` | builds only two slices; reached K_{5,8} at k = 4 (4.6 million cells in all) |

`gbt.reduced_cell_count(G, k)` tells you which row you are in before you build anything.
These limits are for a machine with 7 GB of memory.

**How a full computation runs.** The cells of each degree are enumerated and indexed, and the
boundary matrices are built column by column. Then the discrete Morse reduction cancels every
pair joined by a ±1 entry, which is a unimodular change of basis and so changes neither the
homology nor the non-unit invariant factors. Cancelling a pair (c in degree d, r in degree
d−1) takes three updates: delete column c from ∂_d, clear row r from the other columns of
∂_d, and delete column r from ∂_{d−1}. Forgetting the third gives a complex whose homology is
simply wrong. What is left is usually a few hundred rows across, and its Smith form is
computed exactly by `gbt.smith`. For K_7 at k = 3, 25,256 cells reduce to 370.

**How a single-degree computation runs.** `lean.torsion_in_degree(G, k, d)` builds the slices
of degrees d and d+1 only, forms ∂_{d+1}, peels unit pivots with `linalg.peel_sweeps`, and
takes the Smith form of what remains. It never looks at ∂_d, which only changes the free rank
(Section 3). For the icosahedron at k = 4 (1.3 million cells) this takes about three minutes,
while a heap-based peel runs out of 7 GB.

**Partial reduction.** `swiatkowski(G, k, reduce_at=R)` reduces at the vertices in R only.

| `reduce_at` | what you get |
|---|---|
| `()` | the unreduced complex: biggest, but every state is visible. For K_{n+1,b} the top differential visibly contains the inclusion matrix W_{n−1,n}(b) beside an identity block |
| `range(G.N)` | the fully reduced complex of An–Drummond-Cole–Knudsen |
| one side of a bipartite graph | for K_{n+1,b}: reducing at the B side keeps the inclusion-matrix rows; reducing at the A side destroys them |

All choices give the same homology, and `model_sanity.py` checks that.

## 6. The experiment scripts: conventions shared by all of them

Each file in `experiments/` is one self-contained question. It starts with a docstring that
says what is being tested and why, what output to expect, how to call it and roughly what it
costs. Every option has a default, so `python3 experiments/NAME.py` always does something
sensible.

* **No script imports another**, with one exception: `inventory.py` loads the graph lists
  defined in five other scripts rather than duplicating them. Everything else that is shared
  lives in `gbt/`.
* **Output.** Progress and results go to standard output as plain text tables. Groups are
  printed as primary decompositions (`(Z/2)^4 + Z/3`) unless the invariant factors themselves
  matter.
* **Data.** A script that computes homology stores its results in `data/NAME.json` (Section
  11), keyed so that a rerun recognises what is already there. The sweeps write after every
  completed job and skip stored jobs on restart, so an interrupted run loses at most one case.
  A rerun with a smaller `--cap` adds to the stored results and never discards larger cases.
  `make_tables.py`, `prime_census.py --skip-compute` and `inventory.py` only read the stored
  data.
* **Exit codes.** A script that tests a prediction exits with 1 when the prediction fails and
  0 otherwise; the failure is also printed. This holds for `validate`, `validate_slow`,
  `model_sanity`, `verify_headlines`, `crosscheck_mod_p`, `wilson_block`, `wilson_cokernel`,
  `top_degree_vanishing`, `torsion_class_X`, `prime_census`, `planar_odd_search`,
  `higher_degree_torsion`, `theta_class`, `octahedron_generator`, `decomposability_K8`,
  `large_cases`, and `supergraph_walks` (1 when a counterexample to "no 2-torsion is
  minor-closed" turns up). Searches that merely report what they find (for instance
  `genuine_higher_torsion` and `torsion_hunt`) exit with 0 either way.
* **Long runs** should be detached, and run unbuffered so that the log is current:

      nohup python3 -u experiments/sweep_complete.py --cap 1600000 > sweep.log 2>&1 &

The times below were measured on one core of a cloud machine with 7 GB of memory, and are
rough.

## 7. Experiments on K_{n+1,b}: the inclusion-matrix mechanism

| script | usage | what it computes | time |
|---|---|---|---|
| `top_degree_vanishing.py` | `[--cap 400000] [--amax 7] [--bmax 10] [--nmax 5]` | H_n(B_n K_{a,b}) for every triple (a, b, n) under the cap, checked against the prediction: zero if a ≤ n, or if a = n+1 and b ≤ 2n−1; rank exactly max(0, C(b,n) − C(b,n−1)) when a = n+1; and nonzero in every other case. That last statement is a conjecture, and this script is the evidence for it. Then prints the a = n+1 family with C(b,n−1), C(b,n) side by side, and the torsion one degree down | minutes |
| `torsion_class_X.py` | `[n,N ...]` (default: nine cases) | builds the chain X of `gbt.flipchains` for K_{n+1,N} straight from the closed formula, applies ∂ term by term without building the complex, and checks (P1) the parked part of ∂X vanishes and (P2) the rest is (N−n+1)·y with all coefficients of y equal to ±1. At N = 2n−1 this gives Z/n in H_{n−1}. It also runs the competing ordering in the definition of σ (by A-value instead of by B-source) and shows it fails (18 parked terms survive on K_{3,3}; on K_{4,5} the gcd collapses to 1) | seconds to minutes; K_{5,8} is the slow case |
| `wilson_block.py` | `[n,b ...]` (default (2,3), (2,4), (2,5), (3,4), (3,5), (3,6)) | checks, case by case on the unreduced complex: (1) ∂X_{B₀} = (−1)^{n−1} Σ_T Y_T for every B₀; (2) X_{B₀} has the closed-form coefficients; (3) the parked and B-sitting rows have all invariant factors 1, and kernel of rank C(b,n), so spanned by the X_{B₀}; (4) the Y_T have disjoint supports and entries ±1; (5) the non-unit invariant factors of ∂_n are those of W_{n−1,n}(b); (6) the completion of the row basis can be taken to be standard basis vectors in these cases. Prints `ALL CHECKS PASSED` | about a minute |
| `wilson_cokernel.py` | `[--cap 150000] [--matrix-cap 20000]` | (1) Wilson's diagonal form against an exact Smith form of W_{t,k}(m), compared on elementary divisors because Wilson's list is not a divisibility chain (for W_{3,5}(9) it reads 10, 6^8, 3^27, 1^48); (2) Gottlieb's rank theorem, and with it b ≤ 2n−1 ⇔ ker W_{n−1,n}(b) = 0; (3) coker W_{n−1,n}(b) against the torsion of H_{n−1}(B_n K_{n+1,b}) computed from the complex, for every (n, b) under the cap | a couple of minutes |
| `prime_census.py` | `[--skip-compute]` | (1) computes tors H_3(B_4 K_{5,7}) and tors H_3(B_4 K_{5,8}) directly with `lean.torsion_in_degree` (2.6 and 4.6 million cells in the full complexes) and compares them with coker W_{3,4}(b). These are the Z/4 cases; (2) takes a census of every torsion summand stored in `data/`, by primary type and graph family: of 955 stored summands, 29 are not Z/2, and the script checks by subgraph search that every one occurs in a graph containing K_{4,5} | under a minute |
| `large_cases.py` | `--list`, or `NAME ...` [`--force`] | the individually expensive full computations, one at a time with `gbt.lean`: `K55` (K_{5,5}, k = 4, 617,235 cells), `K56` (K_{5,6}, k = 4, 1,355,675), `K8-4`, `K12-3`, and two that do not fit in 7 GB, `K57` (2,615,835 cells) and `K13-3` (1,639,573). `--list` shows sizes, whether each fits, and whether it is stored | hours each |
| `density_interpolation.py` | `[--k 3]` | adds edges inside the large side of K_{4,5}, one at a time, up to K_9 on the same nine vertices, and watches the Z/3 disappear. The Z/3 is gone once the large side is filled in, although K_9 contains K_{4,5} | seconds |

`top_degree_vanishing.py` and `sweep_bipartite.py` overlap. The first asks only about the top
degree, and so reaches further.

## 8. Sweeps over families of graphs

| script | usage | what it computes | time |
|---|---|---|---|
| `sweep_complete.py` | `[--cap 400000] [--nmax 13] [--kmax 12]` | full H_*(B_k K_n) for every (n, k) under the cap, cheapest first. Stored: K_3 to K_12, k from 2 up to 12 for small n. Every torsion summand is Z/2; for 8 ≤ n ≤ 12 at k = 3 the multiplicity in H_2 is exactly rank H_1 = C(n−1, 2) | minutes at the default cap; hours at `--cap 1600000` |
| `sweep_bipartite.py` | `[--cap 400000] [--amax 7] [--bmax 8] [--kmax 7]` | full H_*(B_k K_{a,b}). Lists every case with odd torsion, and the "background": in most cases the 2-torsion of H_{k−1} has multiplicity exactly rank H_1, and the interesting graphs are those that deviate. The stored data used `--cap 700000` | minutes; about an hour at 700,000 |
| `sweep_small_graphs.py` | `[--vertices 6] [--k 3,4] [--mindeg 3]` | every connected graph on the given number of vertices with minimum degree at least `--mindeg`, from the networkx graph atlas (complete up to seven vertices), at each k. Writes `data/small_graphs_n{N}_k{K}_d{D}.json`. `--mindeg 2` loses only subdivisions. `--mindeg 3` also loses graphs with leaves, which **do** matter (Section 9), so the stored exhaustive sweeps use `--mindeg 1` | seconds (≤ 6 vertices) to minutes (7 vertices) |
| `planar_torsion.py` | `[--k 3,4]` | planar graphs with torsion in H_2 (the octahedron at k = 4), next to the three planar "atoms" of An–Knudsen (the cycle, the star and the theta graph), whose H_2 is torsion-free | seconds |
| `platonic_solids.py` | `[--kmax 3] [--include-k4] [--extras]` | the tetrahedron, octahedron, dodecahedron and icosahedron. With `--include-k4`, the k = 4 cases, using the single-degree route for the large ones: the octahedron, dodecahedron and icosahedron have Z/2 in H_2(B_4), the tetrahedron none. With `--extras`, the prism and antiprism families at k = 4 | seconds at k ≤ 3; minutes (dodecahedron) to much longer with several GB (icosahedron) at k = 4 |
| `planar_odd_search.py` | `[--kmax 4] [--cap 600000]` | a search for torsion other than Z/2 in a planar graph, in the places where the inclusion-matrix mechanism comes closest to being available: generalised thetas K_{2,m}, wheels, bipyramids, grids, maximal planar graphs, the tetrahedron, octahedron, dodecahedron and icosahedron, prisms and antiprisms (28 graphs). Every invariant factor found is 2. Exits with 1 if a counterexample turns up | a few minutes |
| `higher_degree_torsion.py` | `[--max-degree 4]` | Z/2 in H_3 and H_4 of planar graphs, from Künneth products of the octahedron's class with disjoint triangles, including a connected example (octahedron and triangle joined by an edge, Z/2 in H_3(B_5)). The octahedron alone stays torsion-free in H_3 at k = 5, 6 | about a minute |
| `genuine_higher_torsion.py` | `[--k 5] [--degree 3] [--cap 400000]`, or `--deep` | torsion in H_3 of planar graphs on at most 8 vertices, where no Künneth product can produce it: all 55 planar graphs of minimum degree 3 on 6 and 7 vertices, and a sample of 10 on 8 vertices. 59 fit the cap at k = 5 and 52 at k = 6. `--deep` takes the octahedron and the 4-antiprism further. Nothing found. The search does not cover graphs with leaves | minutes at k = 5; about an hour at k = 6 |
| `petersen_family.py` | `[--k 3]` | the seven Petersen-family graphs: H_2 torsion-free at k = 3, exactly Z/2 at k = 4. So H_2 torsion does not detect intrinsic linking | minutes at k = 3; hours at k = 4 |
| `supergraph_walks.py` | `[--start NAME ...] [--walks 4] [--steps 8] [--cap 600000] [--nmax 10]` | random walks that enlarge a known 2-torsion carrier (the octahedron, the 7-vertex carriers, the pendant-edge example, the 4-antiprism, the 5-prism, the 9-vertex triangulation, K_{4,5}, K_7) one random edge or new vertex at a time, keeping the start as a subgraph, and recompute tors H_2 after each step. Tests whether "no 2-torsion" is minor-closed. 470 graphs up to 10 vertices: the 2-torsion never died. This is also how the Z/4 in H_2 was found | about two hours as run |
| `torsion_hunt.py` | `[--family F ...] [--sample 600] [--cap 400000] [--timeout 600] [--workers 2]` | the large search for (1) odd torsion in planar graphs and (2) summands of order ≥ 4 in H_2. Families: all planar graphs of minimum degree 3 on 8 vertices (k = 4, 5) and on 9 vertices (k = 4), 600 on 10 vertices; K_{4,5} plus any graphs on its two sides (k = 3, 4); K_{4,6}, K_{5,5}, K_{5,6}, K_{4,7}, K_{5,7} plus graphs on one side; complete multipartite graphs; random graphs on 9 and 10 vertices; multigraphs (via subdivision). For each graph: build the degree-2 and degree-3 slices, peel unit pivots, then for every prime p ≤ 97 compare ranks over F_p and F_1000003 and take the local Smith form where they differ. Writes `data/torsion_hunt.jsonl`, resumable; hits are marked `!!` | about 12 hours on 2 cores for everything stored |
| `linkless_controls.py` | `[--k 3]` | linklessly embeddable graphs: K_5, K_{3,3}, K_{3,3}+e, apex graphs, maximal planar graphs, the octahedron, and K_6 minus a perfect matching (which is the octahedron again, a consistency check). The apex graph apex(planar 8) has (Z/2)^4 in H_2 already at k = 3 | seconds at k = 3; a few minutes at k = 4 |

## 9. The structure of the torsion classes

These scripts look at individual torsion classes rather than groups. Most use the
"adjoin and test" idiom of Section 14.

| script | usage | what it shows | time |
|---|---|---|---|
| `two_torsion_anatomy.py` | none | (A) the star (exchange) class of a tripod, pushed into a graph at k = 2, has order exactly 2 in K_{3,3} and K_5 and infinite order in planar graphs; (B) in a theta graph with branches of length 3 the product of the two hub star classes generates H_2(B_4 Θ) = Z; (C) that theta class, pushed into the octahedron, has order 2, and has infinite order in the torsion-free triangular prism | a few minutes |
| `theta_class.py` | `[--graph octahedron] [--k 4] [--max-thetas N]` | pushes the generator of H_2 of each embedded theta subgraph into the graph and tests whether it kills the Z/2. For the octahedron it does: the planar torsion is the image of a **free** class, created by relations rather than carried by a new generator. `--graph` accepts any zero-argument constructor in `gbt` | about 20 s for the octahedron; longer for larger graphs |
| `octahedron_generator.py` | `[--strategy random\|prefix\|bisect] [--seed 4] [--tries 12]` | an explicit cycle z = ∂x/2 generating the octahedron's Z/2 in H_2(B_4), found by halving (Section 14), and the fact that its image in K_6 generates the Z/2 of H_2(B_4 K_6). `prefix` is included to document that the obvious search fails here | a few minutes |
| `decomposability_K8.py` | `[--n 8]` | builds every Künneth product class in H_2(B_3 K_n), i.e. K_5 with two particles times a disjoint triangle with one, with the Koszul sign of the merge, and adjoins all of them to im ∂_3. For K_8 the 9,296 products kill the whole (Z/2)^21, so every torsion class there is decomposable. For n ≤ 7 there are no products at all | about twenty minutes and several GB for K_8 |
| `kunneth_K6.py` | `[--n 6] [--k 4]` | lists every splitting of K_6 into two vertex-disjoint complete graphs and every k₁ + k₂ = 4, with the homology of the pieces, and shows that none of the 75 combinations can contribute torsion to H_2. So the Z/2 of H_2(B_4 K_6) is indecomposable | seconds |
| `edge_minimality.py` | `[--case K7-3\|K6-4\|OCT-4\|descent\|stab\|all]` | (1) every single-edge deletion kills the torsion of K_7 at k = 3 (21 deletions) and of the octahedron at k = 4 (12), but none of the 15 deletions from K_6 at k = 4 does; (2) a randomised greedy descent from K_6 that keeps the torsion lands on the octahedron; (3) K_6 has torsion-free H_2 at k = 2, 3 and Z/2 from k = 4 on | seconds to a minute per case |
| `minimal_carriers.py` | none | tests edge-minimality and octahedron containment for the octahedron, the three 7-vertex 12-edge carriers (read from `data/small_graphs_n7_k4_d1.json`), and the graph K_5 + a vertex joined to two vertices of K_5 + a pendant edge (edges (0,1), (0,2), (0,3), (0,4), (1,2), (1,3), (1,4), (2,3), (2,4), (3,4), (3,5), (4,5), (5,6)). All are edge-minimal and only the octahedron contains the octahedron. The last example has Z/2 in H_2(B_4) and loses it when the pendant edge is deleted, so a leaf can create torsion | about a minute |

Edge-minimality is not indecomposability. A graph can be an edge-minimal carrier of a torsion
class and that class can still be the image of a free class on a subgraph, which is exactly
what `theta_class.py` finds for the octahedron.

## 10. Verification and regression

Each script checks part of the software against something independent of it and ends with a
verdict line.

| script | what it checks | time |
|---|---|---|
| `model_sanity.py` | five checks of the machinery, independent of any torsion result: (1) ∂∘∂ = 0; (2) the closed-form cell counts against the built complexes; (3) `reduce_at` does not change the homology; (4) Abrams' model agrees with Świątkowski's on small cases, and reproduces Abrams' surfaces B_2 K_5 = N_7 and B_2 K_{3,3} = N_5 (H_1 = Z^6 ⊕ Z/2 and Z^4 ⊕ Z/2); (5) tors coker ∂_{d+1} = tors H_d in every degree, with free ranks differing by exactly rank ∂_d. Prints `ALL CHECKS PASSED`. Run it first if anything looks surprising | under a minute here; the docstring allows three |
| `validate.py` | 60 regression checks of every headline number that is cheap to recompute, in eight groups: known values from the literature; the top-degree vanishing; the torsion of H_{n−1}(B_n K_{n+1,b}) against coker W; the rank formula for H_n; complete graphs; the planar carriers; the An–Knudsen atoms; internal consistency (`fast` against `lean` against the reference pipeline against Abrams). Run it after changing anything in `gbt/`. Prints `ALL CHECKS PASSED` | under a minute here; up to five on a slower machine |
| `validate_slow.py` | `--list`, or case numbers: the four million-cell cases (K_{5,5} and K_{5,6} at k = 4, K_8 at k = 4, K_12 at k = 3), one at a time | minutes to hours each |
| `crosscheck_mod_p.py` | `[--cap 200000]`: Betti numbers over F_2 and over F_1000003 computed from the **raw** complex, with no Morse reduction and no Smith form, against the prediction of the integral answer through universal coefficients, dim H_d(F_p) = f_d + #{p \| q_{d,i}} + #{p \| q_{d−1,i}}. First on named cases, then on every stored sweep result under the cap | a few seconds to two minutes |
| `verify_headlines.py` | `[--skip-large]`: recomputes the headline torsion groups by two further routes, the single-degree pipeline (23 cases) and Abrams' model on the subdivided graph (6 cases: K_5, K_{3,3}, the octahedron at k = 3 and 4, K_{4,5}, K_7), and compares with the expected values. Prints `ALL HEADLINE VALUES CONFIRMED by both routes.` | about ten minutes; one with `--skip-large` |
| `make_tables.py` | `[--latex]`: regenerates the summary tables (complete graphs by (n, k); the multiplicity against C(n−1, 2); every bipartite case with odd torsion; the a = n+1 family against coker W; the background count; the top-degree trichotomy; and the other stored files) from `data/` alone. With `--latex` the tables come out as LaTeX tabular bodies. With an empty `data/` it prints "(no data -- run the corresponding sweep first)" for each table | instant |

## 11. Data `data/`

All files are JSON, written with sorted keys. In every file a homology dictionary `H` maps a
degree (as a string) to `[free rank, [invariant factors > 1]]`, and `torsion` is a list of
invariant factors > 1. Keys use `|` as a separator.

| file | key | fields | records |
|---|---|---|---|
| `sweep_complete.json` | `K{n}\|{k}` | `H`; `raw`, `red`: total cells before and after the Morse reduction | 46 (K_3 to K_12) |
| `sweep_bipartite.json` | `K{a},{b}\|{k}` | as above | 92 (27 graphs) |
| `large_cases.json` | `K55`, `K56`, `K8-4`, `K12-3` | `H`, `case`, `k`, `raw`, `red` (total cells before and after the Morse reduction; `null` for `K55` and `K56`, whose runs did not record them), `summary` (the torsion as one line) | 4 |
| `top_degree_vanishing.json` | `K{a},{b}\|{n}` | `H_n` (H_n only), `cells`, `predicted` (null where only the conjecture applies) | 81 |
| `prime_census.json` | `K{a},{b}\|{k}\|{d}` | `a`, `b`, `k`, `degree`, `torsion` (of H_d only) | 5 |
| `density_interpolation.json` | `{name}\|{k}` | `H`, `edges`, `note` | 6 |
| `planar_torsion.json` | `{name}\|{k}` | `H`, `vertices`, `edges`, `note` | 18 |
| `planar_odd_search.json` | `{name}\|{k}` | `H`, `N`, `M`, `note` | 100 (28 graphs) |
| `platonic_solids.json` | `{name}\|{k}` | `H` (full computation) or `H2_torsion` (single degree), `N`, `M` | 19 |
| `petersen_family.json` | `{name}\|{k}` | `H` | 14 |
| `linkless_controls.json` | `{name}\|{k}` | `H` | 9 |
| `higher_degree_torsion.json` | `{name}\|{k}\|{d}` | `torsion` (of H_d), `N`, `M`, `k`, `degree`, `cells`, `planar` | 5 |
| `genuine_higher_torsion.json` | `{N}\|{edge list}\|{k}\|{d}`, or `deep\|{name}\|{k}\|{d}` | `torsion` (of H_d), `N`, `M`, `k`, `degree`, `edgelist` | 116 |
| `supergraph_walks.json` | `{start}\|{walk}\|{step}` | `k`, `N`, `M`, `edges`, `slices` (cells in the two slices), `torsion` (`{"2": invariant factors of tors H_2}`, plus higher degrees when H_2 lost its 2-torsion), `seconds` | 470 |
| `torsion_hunt.jsonl` | one JSON object per line | `family`, `name`, `k`, `N`, `M`, `planar`, `edges`, `rank` (of the residual), `primes_checked` (`"<= 97"`), `torsion` (`{p: [number of invariant factors of valuation 0, 1, 2, …]}`), `summary` (tors H_2 as a primary decomposition), `seconds` | 7,291 |
| `two_torsion_anatomy.json` | `A\|{graph}`, `C\|{graph}`, and `B` | `order` of the pushed class (`"2"` or `"infinite"`); for `B`, `is_cycle` and `generates` | 8 |
| `small_graphs_n{N}_k{K}_d{D}.json` | (a list) | one entry per (graph, k): `edgelist`, `degrees`, `edges`, `k`, `planar`, `torsion` (a string, `none` if torsion-free), `torsion_above_H1` | see below |
| `inventory.csv` | (a table) | Section 12 | 10,653 |

The small-graph sweeps stored are

| file | graphs | k | carriers of torsion in some H_i, i ≥ 2 |
|---|---|---|---|
| `small_graphs_n4_k2-3-4-5-6_d1.json` | all 6 connected graphs on 4 vertices | 2–6 | none |
| `small_graphs_n5_k2-3-4-5_d1.json` | all 21 on 5 vertices | 2–5 | none |
| `small_graphs_n5_k2-3-4-5_d2.json` | the 11 of minimum degree ≥ 2 | 2–5 | none |
| `small_graphs_n6_k2-3-4_d1.json` | all 112 on 6 vertices | 2–4 | 5 graphs, all at k = 4: the octahedron (the only one with 12 edges), its three supergraphs on six vertices, and one 13-edge graph with a vertex of degree 3 |
| `small_graphs_n6_k4_d3.json` | the 19 of minimum degree ≥ 3 | 4 | 5 |
| `small_graphs_n7_k3_d1.json` | all 853 on 7 vertices | 3 | 1 (K_7) |
| `small_graphs_n7_k4_d1.json` | all 853 on 7 vertices | 4 | 130 |
| `small_graphs_n7_k4_d3.json` | the 150 of minimum degree ≥ 3 | 4 | 97 |

So the exhaustive sweeps cover 6 + 21 + 112 + 853 = 992 graphs.

## 12. The inventory of all computations

`experiments/inventory.py` reads every file in `data/`, rebuilds the graph behind every
record, identifies graphs up to isomorphism (Weisfeiler–Lehman hash, then a full isomorphism
test), and writes `data/inventory.csv` with one row per stored computation. The columns are

| column | meaning |
|---|---|
| `source` | the script or sweep file that produced the record |
| `graph` | the label used by that script |
| `class` | an integer naming the isomorphism class, shared across sources (−1 for the one record not attached to a single graph) |
| `N`, `M` | vertices and edges |
| `k` | number of particles |
| `computed` | `all degrees`, `tors H_d` (the torsion of one degree), `H_n (top) rank only`, or the order of one explicit class |
| `result` | the torsion found, by degree, as a primary decomposition, or `torsion-free` |

It also prints a summary per source. The current totals:

| | |
|---|---|
| records | 10,653 |
| distinct graphs (up to isomorphism) | 8,140 |
| distinct pairs (graph, k) | 10,062, with k from 2 to 12 |
| graphs with torsion in some H_i, i ≥ 2 | 3,991 |
| graphs with a torsion summand other than Z/2 | 207: 206 contain K_{4,5} as a subgraph and one as a minor; the only orders seen are 2, 3 and 4 |
| planar graphs | 5,714, of which 2,032 have torsion in some H_i with i ≥ 2, and none has a summand other than Z/2 |

Of the 8,140 graphs, 992 come from the exhaustive sweeps of Section 11 and most of the rest from
`torsion_hunt.py`. The script takes a
few seconds; rerun it after any new computation.

## 13. Reproducing the main results

Run from the repository root. `data/` already holds the results, so `make_tables.py` and
`inventory.py` reproduce every summary without recomputing anything. The commands below
recompute.

| result | command | time and memory |
|---|---|---|
| everything cheap, in one go | `python3 experiments/validate.py` | under five minutes |
| the machinery is sound | `python3 experiments/model_sanity.py` | under a minute |
| H_2(B_3 K_{4,5}) = Z^155 ⊕ (Z/2)^4 ⊕ Z/3 | `python3 -c "from gbt import *; print(fast.compute(complete_bip(4,5),3)[0])"` | a second |
| Z/4 in H_3(B_4 K_{5,7}) and H_3(B_4 K_{5,8}) | `python3 experiments/prime_census.py` | under a minute |
| torsion of H_{n−1}(B_n K_{n+1,b}) is coker W | `python3 experiments/wilson_cokernel.py`, `python3 experiments/wilson_block.py` | minutes |
| the explicit class X | `python3 experiments/torsion_class_X.py` | minutes |
| top-degree vanishing and its converse | `python3 experiments/top_degree_vanishing.py` | minutes |
| complete graphs | `python3 experiments/sweep_complete.py --cap 400000`, then `--cap 1600000`, then `python3 experiments/large_cases.py K8-4 K12-3` | minutes; hours; hours and a few GB |
| complete bipartite graphs | `python3 experiments/sweep_bipartite.py --cap 700000`, then `python3 experiments/large_cases.py K55 K56` | about an hour; hours and a few GB |
| planar torsion | `python3 experiments/planar_torsion.py --k 3,4` and `python3 experiments/platonic_solids.py --include-k4 --extras` | seconds; up to hours and several GB (icosahedron) |
| no odd planar torsion found | `python3 experiments/planar_odd_search.py --kmax 4` | a few minutes |
| smallest carriers | `python3 experiments/sweep_small_graphs.py --vertices N --k K --mindeg 1` for (N, K) = (4, "2,3,4,5,6"), (5, "2,3,4,5"), (6, "2,3,4"), (7, 3), (7, 4); then `python3 experiments/minimal_carriers.py` | seconds to minutes each |
| the structure of the classes | `two_torsion_anatomy.py`, `theta_class.py --graph octahedron`, `octahedron_generator.py`, `decomposability_K8.py`, `kunneth_K6.py`, `edge_minimality.py` | minutes each; `decomposability_K8` about twenty minutes |
| higher degrees | `higher_degree_torsion.py`, `genuine_higher_torsion.py --k 5`, `genuine_higher_torsion.py --deep` | minutes to an hour |
| intrinsic linking | `petersen_family.py --k 3,4`, `linkless_controls.py --k 3,4` | hours (Petersen graph at k = 4) |
| independent checks | `crosscheck_mod_p.py`, `verify_headlines.py` | minutes; about ten minutes |
| Z/4 in H_2, and the planar search | `python3 experiments/torsion_hunt.py` (needs nauty and numpy) | about 12 hours on 2 cores |
| supergraph walks | `python3 experiments/supergraph_walks.py --walks 6 --steps 10` | about 2 hours |
| summaries | `make_tables.py`, `inventory.py` | seconds |

## 14. Limits, pitfalls and practical advice

* **The peel is the ceiling for the full homology.** `ChainComplex.reduce_fully`, `fast.morse`
  and `lean.reduce_top_down` are Python dictionary surgery. Above roughly two million cells the
  peel dominates the run time, and above four million it exhausts 7 GB. K_{5,7} at k = 4 (2.6
  million cells) and K_13 at k = 3 (1.6 million) did not finish for this reason. For the torsion
  of one degree use `lean.torsion_in_degree`, which reaches much further. The Z/4 of K_{5,7}
  was computed that way, and independently from the explicit chain of `torsion_class_X.py`.
  A compiled peel is the obvious next step for anyone who wants to go further.
* **Torsion of H_d comes from ∂_{d+1}.** The torsion of H_d is read off the invariant factors
  of ∂_{d+1}, not of ∂_d. Getting this wrong shifts every reported group by one degree.
* **A Morse cancellation needs three updates.** Cancelling (c, r) requires deleting column r
  from ∂_{d−1} as well as column c from ∂_d. Forgetting this gives a complex whose homology is
  wrong. The cross-checks in Section 10 exist because this happened once.
* **Wilson's form is diagonal, not a divisibility chain.** Compare it with a Smith form on
  elementary divisors (the primary decomposition), never on the raw lists.
* **Adjoin and test.** To decide whether a cycle z is a nonzero torsion class in H_d, compare
  `coker_torsion` of the columns of ∂_{d+1} with and without z adjoined. If the torsion drops,
  [z] was nonzero; if it does not move, z was already a boundary. With a family of columns
  this tests whether the family *generates* the torsion. `decomposability_K8.py`,
  `theta_class.py`, `two_torsion_anatomy.py` and `octahedron_generator.py` use it.
* **Halving, and why half a basis is not enough.** To produce a torsion class by hand, find x
  with ∂x ≡ 0 mod 2 but ∂x ≠ 0 over Z; then z = ∂x/2 is a cycle with 2z a boundary. The map
  x ↦ [∂x/2] from ker(∂ mod 2) onto H_d[2] is linear and onto, so half of the kernel works,
  but not half of a *basis*: in the octahedron the first fourteen basis vectors and the first
  twelve prefix sums all landed in the wrong coset. Use random kernel elements.
* **Leaves matter.** A vertex of degree 2 can be smoothed away without changing B_k Γ, but a
  vertex of degree 1 cannot. Restricting a sweep to minimum degree 3 loses graphs whose torsion
  needs a pendant edge (`minimal_carriers.py`). The stored exhaustive sweeps use `--mindeg 1`.
* **`gbt.fast` and `gbt.lean` need every vertex to have degree at least 1.** Graphs with an
  isolated vertex must go through `swiatkowski(G, k, reduce_at=())`.
* **`gbt.lean.BASE`** must exceed N·M. It is 2^11, enough for every graph here; raise it for
  larger graphs.
* **Abrams' model does not check its hypothesis.** Subdivide first: every path between
  vertices of degree ≠ 2 needs length at least k−1 and every cycle length at least k+1.
  `verify_headlines.py` shows the subdivisions that suffice for the cases it runs.
* **Edge-minimal does not mean indecomposable** (Section 9).
* **The exact Smith form can stall.** `gbt.smith` (used by `sparse_smith`, and so by
  `lean.torsion_in_degree` and the full pipelines) can run for hours on a residual of only a
  dozen columns when its entries grow, as on K_{4,5} plus a triangle at k = 3. The local Smith
  form `gbt.linalg.local_smith` never takes a gcd and finishes such cases in a fraction of a
  second, one prime at a time; `torsion_hunt.py` uses it for every prime p ≤ 97.
* **Stored data and `--cap`.** The sweeps skip stored jobs, so to recompute a case delete its
  key from the JSON file, or for `large_cases.py` pass `--force`.

## 15. Example scripts `examples/`, the Makefile, and the logs `logs/`

### `examples/`

Four short scripts, written to be read as well as run, and used by the lessons of
[TUTORIAL.md](TUTORIAL.md). Each starts with a docstring explaining what it is for, how to call
it, and what its output means.

| script | usage | what it does | time |
|---|---|---|---|
| `compute.py` | `GRAPH K [--degree D] [--invariant-factors]` | the homology of B_K of a named graph (`complete:7`, `bipartite:4,5`, `octahedron`, `prism:5`, `wheel:6`, `theta:3,2`, `petersen`, …) or of an edge list such as `"0-1 1-2 2-0"`. Prints the cell counts first, chooses `gbt.fast` or `gbt.lean` by size, and with `--degree D` computes only tors H_D from two slices. The quickest way to get an answer with no programming | seconds for most graphs |
| `look_at_the_complex.py` | `[GRAPH] [K] [--max-cells N]` | prints every generator of the Świątkowski complex in words ("parked: e2={0,3} \| moving: 0 along e1={0,2}") with its boundary, checks ∂∘∂ = 0, computes the homology, then does the same for the complex reduced at every vertex of degree ≥ 2 and checks that the homology agrees. Default: the tripod with two particles, whose H_1 = Z is the exchange loop | a second |
| `explicit_torsion_class.py` | none | the three steps used for explicit classes, on the Z/2 of H_1(B_2 K_5): halving (z = ∂x/2 for x in ker ∂ mod 2), adjoin and test (the torsion of the cokernel drops when z is added as a column), and push forward along K_5 ⊂ K_6 | under a second |
| `sweep_template.py` | `[--nmax 8] [--kmax 4] [--cap 300000] [--out PATH]` | a template for a new experiment: a family of graphs (wheels, as shipped) sorted by size, a cap, results saved after every job and skipped on rerun, and a summary. Edit `jobs()` to study another family. Writes to `examples/output/` (ignored by git), never to `data/` | two seconds as shipped |

`look_at_the_complex.py` exports `describe_cell(G, cell, reduced_at=())`, which turns a cell
`(mono, st)` into the text above; `compute.py` exports `parse_graph(text)`. Both can be imported
by your own scripts if `examples/` is on `sys.path`.

### The Makefile

The package is pure Python and nothing is compiled; the Makefile only collects commands.

| target | runs | time |
|---|---|---|
| `make help` | prints the list below | |
| `make examples` | the four example scripts (and removes their scratch output) | seconds |
| `make check` | `model_sanity`, `validate`, `wilson_block`, `wilson_cokernel`, `crosscheck_mod_p` (default cap) and `verify_headlines --skip-large` | about 5 minutes |
| `make check-full` | `check`, then `verify_headlines` in full and `validate_slow 0` (K_{5,5} at k = 4) | 20 minutes or more |
| `make tables` | `make_tables.py` | instant |
| `make inventory` | `inventory.py` (rebuilds `data/inventory.csv`) | seconds |
| `make hunt` | `torsion_hunt.py` with its default (first-phase) families; resumes from `data/torsion_hunt.jsonl`, so with the shipped data it only prints the list of hits | about two minutes with the shipped data; hours from scratch |
| `make logs` | the scripts of `check` (with `verify_headlines` in full), `inventory` and `prime_census --skip-compute`, each writing its output to `logs/NAME.txt` | about 10 minutes |
| `make clean` | removes `__pycache__` folders and `examples/output/` | |

Every target accepts `PYTHON=...` to choose the interpreter, e.g. `make check PYTHON=python3.11`.

### `logs/`

The output of `make logs` from the final run of the package: `model_sanity.txt`,
`validate.txt`, `wilson_block.txt`, `wilson_cokernel.txt`, `crosscheck_mod_p.txt`,
`verify_headlines.txt`, `inventory.txt` and `prime_census.txt`. Each ends with its verdict line.
Rerun `make logs` and compare with `diff` to check a new installation or a change to `gbt/`
(timings in the logs will differ; nothing else should).

## 16. References

The works the code relies on directly are

* J. Świątkowski, *Estimates for homological dimension of configuration spaces of graphs*,
  Colloquium Mathematicum 89 (2001), 69–79. The chain complex.
* B. H. An, G. C. Drummond-Cole and B. Knudsen, *Subdivisional spaces and graph braid groups*,
  Documenta Mathematica 24 (2019), 1513–1583. The complex as Z[E] ⊗ ⊗_v S_v, the reduced
  version, and Proposition 5.21: edge stabilization H_i(B_k Γ) → H_i(B_{k+1} Γ) is injective,
  so torsion present at k persists at every larger k.
* A. Abrams, *Configuration spaces and braid groups of graphs*, PhD thesis, University of
  California, Berkeley (2000). The discretised model and the surfaces B_2 K_5 and B_2 K_{3,3}.
* P. Prue and T. Scrimshaw, *Abrams's stable equivalence for graph braid groups*, Topology and
  its Applications 178 (2014), 136–145. The subdivision hypothesis for Abrams' model.
* R. Ghrist, *Configuration spaces and braid groups on graphs in robotics*, in Knots, Braids,
  and Mapping Class Groups (New York, 1998), AMS/IP Studies in Advanced Mathematics 24 (2001),
  29–40.
* K. H. Ko and H. W. Park, *Characteristics of graph braid groups*, Discrete and Computational
  Geometry 48 (2012), 915–963. Torsion in H_1.
* B. H. An and B. Knudsen, *On the second homology of planar graph braid groups*, Journal of
  Topology 15 (2022), 666–691. The three planar atoms.
* R. M. Wilson, *A diagonal form for the incidence matrices of t-subsets vs. k-subsets*,
  European Journal of Combinatorics 11 (1990), 609–615.
* D. H. Gottlieb, *A certain class of incidence matrices*, Proceedings of the American
  Mathematical Society 17 (1966), 1233–1237.
* R. Forman, *Morse theory for cell complexes*, Advances in Mathematics 134 (1998), 90–145.
  The discrete Morse reduction.
* N. Robertson, P. Seymour and R. Thomas, *Sachs' linkless embedding conjecture*, Journal of
  Combinatorial Theory, Series B 64 (1995), 185–227. The Petersen family.
* A. Hagberg, D. Schult and P. Swart, *Exploring network structure, dynamics, and function
  using NetworkX*, Proceedings of the 7th Python in Science Conference (2008), 11–15. The
  graph atlas, planarity and isomorphism tests.
