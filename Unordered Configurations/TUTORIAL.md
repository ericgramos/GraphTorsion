# Tutorial: computing and understanding torsion in graph braid groups

This tutorial is a guided path through the package. It is written for a reader who knows the
mathematics (configuration spaces of graphs, cellular homology, Smith normal forms) and is
comfortable typing commands into a terminal and editing a short Python script, but who is not
a programmer.

Each lesson says what question it answers, gives commands that work exactly as written from the
top folder of the package, shows what the output looks like, and explains how to read it. The
README is the reference manual; this file is the walk-through. Times and memory are for one
core of an ordinary computer; nothing in Lessons 1 to 10 needs more than 1 GB.

Contents

1. Setting up
2. Graphs: how to describe the graph you have in mind
3. The homology of B_k(Γ), from the command line and from Python
4. Looking inside the chain complex
5. Large cases: the torsion of a single degree
6. Checking an answer by an independent route
7. Writing down a torsion class explicitly
8. The bipartite mechanism: torsion from an inclusion matrix
9. The stored results, and how to search them
10. Writing your own experiment
11. Checking the software yourself
12. Practical advice and common mistakes

---

## 1. Setting up

You need Python 3.8 or later. To find out which version you have, type

    python3 --version

The library itself needs nothing else. A few scripts (graph enumeration, planarity and
isomorphism tests) also need the package `networkx`:

    pip install networkx            # or: pip install --user networkx

Unpack the package and **work from its top folder** (the one containing `README.md`,
`TUTORIAL.md`, `gbt/`, `experiments/`, `examples/` and `data/`). Every command in this tutorial
is typed there:

    cd graph-braid-torsion
    make examples                   # runs the four example scripts; a few seconds

If `make` is not installed (some Windows machines), every `make` target is a short list of
`python3 ...` commands that you can read in the file `Makefile` and type yourself.

There are three ways to use the package, and the lessons use all of them.

* **Run a script.** `python3 examples/compute.py complete:5 2` runs a file and prints its
  answer. Every script explains itself at the top of the file (open it in any text editor), and
  `python3 SCRIPT --help` lists its options.
* **Type a one-line Python command.** `python3 -c "from gbt import *; print(...)"` runs one line
  of Python. The quotes matter: everything between them is the program.
* **Work interactively.** Type `python3` on its own. The prompt changes to `>>>`, and every line
  you type is executed at once. This is the best way to explore. Leave with `exit()`. In this
  tutorial, lines starting with `>>>` are typed at that prompt and the lines after them are what
  Python prints.

Python finds the package `gbt` because you started it from the top folder. If you start Python
somewhere else, `import gbt` fails with `ModuleNotFoundError`; either go back to the top folder or
set `export PYTHONPATH=/path/to/graph-braid-torsion` first. The scripts in `experiments/` and
`examples/` do not have this problem: they find the package themselves.

## 2. Graphs: how to describe the graph you have in mind

**The question.** How do I tell the program which graph I mean?

A graph is a `Graph` object. Its vertices are 0, 1, …, N−1, and its edges are a list of pairs.
The **position of an edge in that list is its name** everywhere else in the package (edge 0,
edge 1, …). The usual families are built in:

    >>> from gbt import *
    >>> G = complete_bip(2, 3)
    >>> G.N, G.M
    (5, 6)
    >>> G.edges
    [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4)]
    >>> G.deg(0), G.nbrs(0)
    (3, [2, 3, 4])

So in K_{2,3} the first side is {0, 1}, the second side is {2, 3, 4}, and edge 4 is {1, 3}.
The constructors are

| constructor | graph |
|---|---|
| `complete(n)` | K_n |
| `complete_bip(a, b)` | K_{a,b}, with A = {0, …, a−1} first, then B |
| `octahedron()`, `tetrahedron()`, `dodecahedron()`, `icosahedron()` | Platonic solids |
| `prism_graph(n)`, `antiprism(n)` | the n-prism and n-antiprism |
| `wheel(n)`, `theta(paths, length)`, `petersen()` | wheels, theta graphs, the Petersen graph |
| `petersen_family()` | the seven Petersen-family graphs, as (name, graph) pairs (needs networkx) |

For anything else, give the edge list. Vertices may be any integers; they are renumbered
0, …, N−1 in increasing order:

    >>> G = from_edges([(10, 20), (20, 30), (30, 10), (30, 40)])     # a triangle with a tail
    >>> G.N, G.edges
    (4, [(0, 1), (1, 2), (0, 2), (2, 3)])

Two more tools: `subdivide(G, t)` replaces every edge by a path of t edges (the configuration
space does not change), and `disjoint_union(G, H)` puts two graphs side by side. If you have a
`networkx` graph `Gx` already, `from_networkx(Gx)` converts it.

**What the package does not handle.** Loops and multiple edges: graphs here are simple.

**Why the order matters.** The chain complex depends on the numbering: the signs come from the
order of the vertices, and the reduction at a vertex uses its lowest-numbered edge. The homology
does not depend on any of this, but explicit chains (Lessons 4 and 7) do. If you compare chains
between two graphs, build them with the same numbering.

## 3. The homology of B_k(Γ), from the command line and from Python

**The question.** What is H_*(B_k Γ; Z) for this graph and this k?

From the command line, with no programming at all:

    python3 examples/compute.py bipartite:4,5 3

prints

    graph: 9 vertices, 20 edges;  k = 3
    cells of the fully reduced complex by degree: {0: 1540, 1: 6510, 2: 8520, 3: 3406}  (19,976 in all)
    pipeline: gbt.fast (full homology)
    H_0 = Z
    H_1 = Z^12 + Z/2
    H_2 = Z^155 + (Z/2)^4 + Z/3
    H_3 = 0
    (0.2 s)

The graph can also be an edge list, in quotes: `python3 examples/compute.py "0-1 1-2 2-0 2-3" 3`.
Run `python3 examples/compute.py --help` for the list of names (`complete:7`, `bipartite:4,5`,
`octahedron`, `prism:5`, `wheel:6`, `theta:3,2`, …).

**From Python.** The same computation:

    >>> from gbt import *
    >>> H, raw, reduced = fast.compute(complete_bip(4, 5), 3)
    >>> H
    {0: (1, []), 1: (12, [2]), 2: (155, [2, 2, 2, 6]), 3: (0, [])}

`H` is a dictionary from the degree d to a pair (free rank, list of invariant factors bigger
than 1). So H_2 = Z^155 ⊕ Z/2 ⊕ Z/2 ⊕ Z/2 ⊕ Z/6. The invariant factors form a divisibility
chain, which hides the primes: Z/6 is Z/2 ⊕ Z/3. Two functions print groups:

    >>> fmt(*H[2])                  # invariant factors
    'Z^155 + (Z/2)^3 + Z/6'
    >>> primary(H[2][1])            # primary decomposition: usually what you want
    '(Z/2)^4 + Z/3'

The other two return values, `raw` and `reduced`, are the number of cells in each degree before
and after the discrete Morse reduction (README Section 5): the computation shrinks a complex of
19,976 cells to a few hundred before taking a Smith normal form.

**How long will it take?** Ask before you build anything; the count is a closed formula:

    >>> reduced_cell_count(complete(8), 4)
    1065001

As a rule of thumb, `gbt.fast` handles up to about 600,000 cells in seconds to minutes, and
`gbt.lean.compute` (same answers, a third of the memory) up to about 2 million. Beyond that, use
Lesson 5. Some reference points:

| graph, k | cells | time |
|---|---|---|
| K_{4,5}, 3 | 19,976 | 0.2 s |
| octahedron, 4 | 26,142 | 0.3 s |
| K_6, 4 | 71,220 | 2 s |
| K_8, 3 | 63,868 | 2 s |
| K_7, 4 | 307,636 | 30 s |

## 4. Looking inside the chain complex

**The question.** What are the cells, and what is the differential, concretely?

    python3 examples/look_at_the_complex.py

builds the Świątkowski complex of the tripod K_{1,3} with two particles and prints every
generator in words with its boundary. A typical entry is

    [deg 1]  sitting: 2 | moving: 3 along e2={0,3}
             d = -1 * (sitting: 2, 3)
                 +1 * (parked: e2={0,3} | sitting: 2)

Read it as: one particle sits at vertex 2, and a particle at vertex 3 moves out along edge 2.
Its boundary is "the move is finished" (the particle is parked inside edge 2) minus "the move
has not started" (the particle sits at 3). That is the formula ∂h_{v,e} = e·∅_v − v of the
README, Section 3. The script then checks ∂∘∂ = 0, prints H_1 = Z (the exchange loop of two
particles at the centre of a "Y"), reduces the complex at the central vertex, prints the
reduced cells, and checks that the homology has not changed. Try it on other small graphs:

    python3 examples/look_at_the_complex.py "0-1 1-2 2-0" 2       # a triangle: B_2 is a circle

**Inside Python.** The reference pipeline keeps the cells as readable labels:

    >>> C = swiatkowski(complete_bip(3, 3), 2)        # reduce_at=() : the unreduced complex
    >>> C.sizes()
    {0: 114, 1: 252, 2: 135}
    >>> c = ((0,), (0, 0, 0, 5, 0, 0))
    >>> C.B[1][c]
    {((0, 3), (0, 0, 0, 0, 0, 0)): 1, ((0,), (0, 0, 0, 1, 0, 0)): -1}

A cell is a pair `(mono, st)`. `mono` lists the edges carrying a parked particle (with
repetition). `st` has one entry per vertex: 0 = empty, 1 = a particle sits there, 2 + e = a
particle there moves along edge e. So `c` is "one particle parked on edge 0, and the particle at
vertex 3 moves along edge 3", and `C.B[1][c]` is its boundary. `C.check_d2()` tests ∂∘∂ = 0,
and `C.reduce_fully(); C.homology()` computes the homology.

`swiatkowski(G, k, reduce_at=R)` reduces at the vertices in R only (README Section 5). The
unreduced complex `reduce_at=()` is the one where every state is visible; for K_{n+1,b} it is
the one where the inclusion matrix of Lesson 8 can be seen. Use this pipeline for complexes up to
a few thousand cells; it is slow, but nothing is hidden.

## 5. Large cases: the torsion of a single degree

**The question.** The complex is too big for a full computation, but I only care about the
torsion of one H_d.

The torsion of H_d can be read off the single matrix ∂_{d+1}, because

    tors H_d = tors coker ∂_{d+1}.

(The sequence 0 → H_d → coker ∂_{d+1} → im ∂_d → 0 is exact, and it splits because im ∂_d is a
subgroup of a free group. So ∂_d only affects the free rank.) In the fully reduced complex
every degree can be enumerated on its own, so only the cells of degrees d and d+1 are ever built:

    python3 examples/compute.py bipartite:5,7 4 --degree 3

prints

    cells of the fully reduced complex by degree: {0: 73815, 1: 450660, 2: 967680, 3: 859600, 4: 264080}  (2,615,835 in all)
    pipeline: gbt.lean.torsion_in_degree (torsion of H_3 only; 1,123,680 cells in the two slices)
    tors H_3 = (Z/2)^14 + (Z/3)^6 + Z/4
    (9.3 s)

The full complex (2.6 million cells) does not fit in 7 GB; the two slices take 10 seconds and
half a gigabyte. In Python: `lean.torsion_in_degree(complete_bip(5, 7), 4, 3)` returns the list
of invariant factors.

This route reached complexes of 4.6 million cells (K_{5,8} at k = 4). Its limit is the size of
the two slices, roughly 2 million cells together on a 7 GB machine.

## 6. Checking an answer by an independent route

**The question.** How do I know the computer is right?

There are three independent checks, and you can run each on your own graph.

**(a) Universal coefficients.** Betti numbers over a field need no Smith form and no reduction.
If H_d = Z^{f_d} ⊕ ⊕ Z/q_i, then dim H_d(F_p) = f_d + #{i : p | q_i in degree d} + #{i : p | q_i
in degree d−1}. For K_5 with two particles, H_1 = Z^6 ⊕ Z/2 and H_2 = 0, so over F_2 one expects
dimensions 7 in degree 1 and 1 in degree 2, and over F_3 no change:

    >>> from gbt.linalg import betti_mod_p
    >>> C = swiatkowski(complete(5), 2)
    >>> betti_mod_p(C, 2), betti_mod_p(C, 3)
    ({0: 1, 1: 7, 2: 1}, {0: 1, 1: 6, 2: 0})

`python3 experiments/crosscheck_mod_p.py` does this for every stored result small enough.

**(b) A different model.** Abrams' discretised configuration space is a different cell complex
for the same space, valid once every edge is subdivided enough (cycles of length at least k+1,
paths between vertices of degree ≠ 2 of length at least k−1). It shares no code with the
Świątkowski builders:

    >>> C = abrams(subdivide(octahedron(), 2), 3)       # every edge cut into 2
    >>> C.reduce_fully(); C.homology()
    {0: (1, []), 1: (8, []), 2: (47, []), 3: (0, [])}

which agrees with `fast.compute(octahedron(), 3)`. The function does not check the subdivision
condition for you, and the model grows very quickly; use it on small cases.

**(c) A different pipeline.** `fast.compute`, `lean.compute` and the reference pipeline
(`swiatkowski` + `reduce_fully` + `homology`) must agree; `experiments/validate.py` checks this on
a spread of graphs.

## 7. Writing down a torsion class explicitly

**The question.** H_1(B_2 K_5) has a Z/2. Which cycle is it, and where does it go in a bigger
graph?

    python3 examples/explicit_torsion_class.py

does three things, all in under a second.

1. **Halving.** It takes an element x of the kernel of ∂_2 modulo 2 whose integral boundary is
   not zero. Every coefficient of ∂x is then even, z = ∂x/2 is an integral cycle, and 2z = ∂x is
   a boundary, so [z] has order 1 or 2.
2. **Adjoin and test.** It adds z as an extra column to ∂_2. Adding a column divides out the
   subgroup generated by z, so the torsion of the cokernel drops from Z/2 to 0 exactly when [z]
   generates the Z/2:

       torsion of H_1 (from coker d_2): Z/2
         try 1: d x has 100 nonzero entries, all even -> z = d x / 2 is a 1-cycle: True
         torsion after adjoining z: 0   => [z] generates the Z/2

3. **Push forward.** K_5 sits inside K_6. The inclusion renames edges, so z gives a cycle of
   K_6, and the same test shows it generates the Z/2 of H_1(B_2 K_6).

It also prints the first few terms of z in words, as in Lesson 4: z is a sum of 100 cells in
which one particle sits still and the other moves.

The same three steps, on bigger complexes, are what `experiments/octahedron_generator.py`,
`experiments/theta_class.py` and `experiments/two_torsion_anatomy.py` do. Read the example first;
those scripts are the same idea with more bookkeeping.

**A warning that cost real time.** The map x ↦ [∂x/2] from ker(∂ mod 2) onto the 2-torsion is
linear and onto, so half of all kernel elements work. That does not mean half of a basis works:
for the octahedron the first fourteen basis vectors all landed in the wrong coset. Use random
combinations, as the example does.

## 8. The bipartite mechanism: torsion from an inclusion matrix

**The question.** Where does the Z/3 of K_{4,5} come from?

Let W_{t,k}(m) be the 0/1 matrix of t-subsets against k-subsets of an m-set (entry 1 when one
set contains the other). For K_{n+1,b} with b ≥ n, the top differential ∂_n of the unreduced
complex is, in a suitable basis, an identity block next to W_{n−1,n}(b). Therefore

    tors H_{n−1}(B_n K_{n+1,b}) = tors coker W_{n−1,n}(b).

You can check this yourself in a few lines, with no chain complex at all for the right-hand side:

    >>> import itertools
    >>> def W(t, k, m):
    ...     rows = list(itertools.combinations(range(m), t))
    ...     cols = list(itertools.combinations(range(m), k))
    ...     return [[1 if set(r) <= set(c) else 0 for c in cols] for r in rows], len(rows), len(cols)
    ...
    >>> rank, factors = smith(*W(2, 3, 5))
    >>> primary([d for d in factors if d > 1])
    '(Z/2)^4 + Z/3'
    >>> primary(fast.compute(complete_bip(4, 5), 3)[0][2][1])
    '(Z/2)^4 + Z/3'

(At the `>>>` prompt, a line ending in `:` continues on lines starting with `...`; press Enter
on an empty `...` line to finish.) Wilson's diagonal form then gives the answer in closed form:
summands Z/j for every j ≤ min(n, b − n + 1). The scripts that test this are

| script | what it checks | time |
|---|---|---|
| `experiments/wilson_cokernel.py` | Wilson's diagonal form, Gottlieb's rank theorem, and coker W against the homology, for every (n, b) small enough | 10 s |
| `experiments/wilson_block.py` | the identity block, step by step, on the unreduced complex | 5 s |
| `experiments/torsion_class_X.py` | the explicit chain X with ∂X = (b−n+1)·y (from `gbt.flipchains`), which gives Z/n when b = 2n−1 | seconds |
| `experiments/prime_census.py` | the Z/4 in H_3(B_4 K_{5,7}) and H_3(B_4 K_{5,8}) by Lesson 5, against coker W | under a minute |

## 9. The stored results, and how to search them

**The question.** Has this already been computed?

Every experiment stores its answers in `data/` as JSON, so they can be read back without
recomputing. In each file a homology group is `[free rank, [invariant factors]]`:

    >>> import json
    >>> d = json.load(open("data/sweep_complete.json"))
    >>> fmt(*d["K8|3"]["H"]["2"])
    'Z^734 + (Z/2)^21'

The keys are `GRAPH|k` (README Section 11 lists every file and its format).

**All computations in one table.** `data/inventory.csv` has one row per stored computation
(2,892 rows), with columns `source`, `graph`, `class` (an isomorphism-class number shared across
files), `N`, `M`, `k`, `computed`, `result`. Open it in a spreadsheet, or search it from Python:

    >>> import csv
    >>> rows = list(csv.DictReader(open("data/inventory.csv")))
    >>> odd = [r for r in rows if "Z/3" in r["result"] or "Z/4" in r["result"]]
    >>> len(odd)
    16
    >>> odd[0]["graph"], odd[0]["k"], odd[0]["result"]
    ('K4,5', '3', 'H1:Z/2; H2:(Z/2)^4 + Z/3')

After any new computation, `make inventory` (or `python3 experiments/inventory.py`) rebuilds the
table and prints totals. `make tables` prints the summary tables (complete graphs, complete
bipartite graphs, the top-degree vanishing, …) from the stored data.

## 10. Writing your own experiment

**The question.** I want to run a family of graphs that nobody has run.

Copy `examples/sweep_template.py` and edit one function. As shipped it computes the wheels W_3
to W_8 with 2 to 4 particles:

    python3 examples/sweep_template.py

    W_3    k=2          93 cells  H_1=Z^4  H_2=0                         torsion: none   [0.0s]
    ...
    W_8    k=4      77,076 cells  H_1=Z^9  H_2=Z^124  H_3=0  H_4=0         torsion: none   [0.9s]
    summary: 18 jobs stored; torsion above H_1 in: none

To study another family, change `jobs()`. For example, all antiprisms with three to five
particles:

    def jobs(nmax, kmax):
        out = []
        for n in range(3, nmax + 1):
            for k in range(3, kmax + 1):
                out.append(("antiprism_%d" % n, antiprism(n), k))
        return out

(and add `antiprism` to the line `from gbt import ...` at the top). The rest of the template does
the bookkeeping that every sweep needs: it sorts the jobs by size, skips those larger than
`--cap`, saves after every job, and skips stored jobs when you rerun it, so an interrupted run
resumes where it stopped. Results go to `examples/output/`, not to `data/`, so experimenting
never disturbs the project's stored results; delete that folder to start again.

For a long run, start it so that it survives closing the terminal and writes a log:

    nohup python3 -u examples/sweep_template.py --nmax 12 --cap 2000000 > wheels.log 2>&1 &
    tail -f wheels.log              # watch it; Ctrl-C stops watching, not the run

## 11. Checking the software yourself

| command | what it checks | time |
|---|---|---|
| `make check` | the model (∂∘∂ = 0, cell counts, reductions, Abrams' model), 60 regression values, the bipartite mechanism step by step, mod-p Betti numbers against the integral answers, and the headline groups by two further routes | about 5 minutes |
| `make check-full` | the above, plus the dodecahedron, the icosahedron and K_{5,5} by the second route, and the first million-cell regression case | 20 minutes or more |
| `python3 experiments/validate_slow.py --list` | the four million-cell cases, one at a time | minutes to hours each |

Each script ends with a verdict line (`ALL CHECKS PASSED`, `ALL HEADLINE VALUES CONFIRMED by
both routes.`, …) and exits with a nonzero code if anything disagrees. The outputs of our own
runs are in `logs/`; `make logs` regenerates them, so you can compare line by line with
`diff`.

Run `python3 experiments/model_sanity.py` first whenever something looks surprising: it tests
the machinery without assuming any of the results.

## 12. Practical advice and common mistakes

* **Count before you compute.** `reduced_cell_count(G, k)` is instant. Memory is the limit: a few
  hundred bytes per cell for `gbt.fast`. If a run is killed with no message, or the shell
  reports `Killed`, the computer ran out of memory; use `gbt.lean`, or Lesson 5.
* **The torsion of H_d comes from ∂_{d+1},** the differential coming in, not the one going out.
  Mixing them up shifts every answer by one degree.
* **Invariant factors are not primes.** `[2, 2, 2, 6]` is (Z/2)^4 ⊕ Z/3. Use `primary(...)`.
* **Leaves matter.** A vertex of degree 2 can be removed without changing B_k(Γ), but a vertex of
  degree 1 cannot. Do not restrict a search to minimum degree 3 if you want every graph.
* **A disconnected graph has disconnected configuration spaces.** H_0 counts the ways of
  distributing the k particles among the components, so it can be Z^m with m > 1.
* **Torsion persists in k.** Edge stabilization is injective, so once H_i(B_k Γ) has an element
  of order q, so does H_i(B_{k+1} Γ). Searching small k first loses nothing.
* **Subgraph inclusions are not injective.** K_9 contains K_{4,5}, but the Z/3 of K_{4,5} dies
  in K_9 (`experiments/density_interpolation.py`).
* **Abrams' model needs subdivision**, which the code does not check (Lesson 6).
* **Stored results are reused.** Most experiments skip anything already in `data/`. To force a
  recomputation, delete the entry (or the file), or use the script's `--force` option where it
  has one.
* **Do not kill programs with `pkill -f name`** from a shell whose command line contains the same
  word: it kills the shell too. Use the process number: `ps aux | grep python3`, then
  `kill 12345`.
