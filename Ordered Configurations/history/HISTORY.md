# History of the engine and the working scripts

This folder keeps every intermediate version of the engine and the one-off scripts that
drove the computations, exactly as they were used. Nothing here is needed to use the
package. The current engine is `engine/oconf.cpp`, which is a commented and tidied copy of
version 26 with the same logic (the regression test `validation/regression_vs_reference.py`
checks that the two give identical output on 33 test cases covering every mode).

The old versions are kept for two reasons. First, the results recorded in `data/` were
produced by specific versions, and a reader who wants to repeat a computation exactly as it
was done can build the version that produced it. Second, the sequence shows which ideas
mattered for performance, which is useful to anyone who wants to push these computations
further.

Every version builds with

    g++ -O3 -march=native -std=c++17 -o oconfNN oconfNN.cpp            # versions 1 to 17
    g++ -O3 -march=native -fopenmp -std=c++17 -o oconfNN oconfNN.cpp   # versions 18 and later

Version 26 also needs `lgr.h` (and through it `grdense.h`) from this folder in the same
directory. The two headers were merged into `engine/gr42.h` in the release.

## The engine versions

The names follow the files in `engine_versions/`. Line counts are approximate.

| file | lines | what changed |
|---|---|---|
| `oconf.cpp` | 310 | First version. Builds the reduced ordered Świątkowski complex (Hainaut, Knudsen and Wawrykow, Definition 3.12) with cells stored as strings in a hash map, and eliminates every boundary map over Z with unit pivots. The residual block is written to a file and its Smith form is computed in Python with FLINT. |
| `oconf2.cpp` | 310 | Cells encoded as 128-bit integers (12 bits per particle: 8 bits for the location, 4 bits for the position of the particle on its edge). About five times less memory than strings. |
| `oconf3.cpp` | 320 | Column dropping. Rows of d_k used as unit pivots are removed from the columns of d_{k+1} before d_{k+1} is built. This is valid because after a unit-pivot reduction, ker d_k projects isomorphically onto the remaining coordinates. It roughly halves the size of the second matrix. |
| `oconf4.cpp` | 375 | Option `-z` (check that d∘d = 0 on every cell) and option `-k` (coloured particles, giving the quotient by a Young subgroup of the symmetric group). |
| `oconf6.cpp` | 375 | Matrix entries stored as 64-bit integers with overflow detection, instead of arbitrary precision. (Version 5 was a short-lived experiment whose source was not kept.) |
| `oconf7.cpp` | 405 | Option `-q p K`: local Smith form over Z/p^K. Entries stay bounded, so the computation cannot stall on entry growth, and the stage counts give the exact p-primary torsion up to exponent K-1. Entries stored in 16 bytes. |
| `oconf8.cpp` | 420 | Option `-g` with a single cyclic permutation of the labels: the coinvariant complex for a cyclic group of label permutations. |
| `oconf9.cpp` | 390 | Lean 8-byte matrix entries (32-bit column, 32-bit value). |
| `oconf10.cpp` | 470 | A separate elimination routine over F_2 (class `F2E`) with bit operations. |
| `oconf11.cpp` | 520 | A separate routine over F_p for small primes (class `SPE`, entries in one byte). |
| `oconf12.cpp` | 535 | Option `-T`: eliminate the transpose (pivot on rows instead of columns), which is sometimes much faster. |
| `oconf13.cpp` | 545 | `-g` accepts several generators; the engine closes them up to the whole group. |
| `oconf14.cpp` | 625 | Local engine (class `LPE`) with packed entries `col<<5 | val`, which halves memory in the local mode. |
| `oconf15.cpp` | 625 | Only canonical cells (orbit representatives) are enumerated, instead of enumerating everything and discarding the rest. |
| `oconf16.cpp` | 630 | Elimination diagnostics printed every 200,000 pivots (live rows, fill per row). These logs showed that fill explodes once fewer than about 40,000 live columns remain. |
| `oconf17.cpp` | 720 | The dense finishing phase. When the remaining live part becomes small and dense (option `-D` sets the column threshold, `-C` the fill), the engine switches to a dense streaming row reduction. |
| `oconf18.cpp` | 725 | The dense phase parallelised with OpenMP. |
| `oconf19.cpp` | 810 | Graph automorphisms (option `-a`), including automorphisms that reverse edges (written `~e`), and twisted coefficients (option `-X`) for the sign characters of labels and edges. The group action is applied directly to the 128-bit key. |
| `oconf20.cpp` | 825 | Rows allocated with exact capacity; priority-queue rebuilds; compaction of the column lists; group elements shuffled so that canonicalisation aborts early. `oconf20s.cpp` is a variant used for one timing experiment. |
| `oconf21.cpp` | 895 | Bit-sliced dense phase over Z/4 (two bit planes), the common case p = 2, K = 2. |
| `oconf22.cpp` | 960 | Option `-F`: fast orbit enumeration (canonical shapes under the automorphisms, then canonical labellings under the label group). Canonicalisation factored over the product A×L. The matrix build runs in parallel. The guard that refuses a group whose order is divisible by p (exit code 10) was added here after an earlier mistake produced fake 3-torsion from a group of order 24. |
| `oconf23.cpp` | 1115 | Bit-sliced dense phase over Z/9 (four bit planes, with a ripple-carry adder that reduces mod 9). |
| `oconf24.cpp` | 1095 | Option `-M`: start the dense phase when the sparse part would exceed a memory budget. Option `-Y`: an arbitrary sign for each generator of the group. |
| `oconf25.cpp` | 1115 | Online compression of leftover rows in the dense phase: each new leftover row is reduced against the span of the previous ones mod p. Before this, the Z/9 phase ran out of memory on the T(1,3,3)+apex run at n = 7. |
| `oconf26.cpp` | 1130 | Option `-W`: coefficients in the Galois ring GR(4,2) = Z/4[ω] with ω³ = 1, for characters of order 3 (uses `lgr.h` and `grdense.h`). |

## The scripts in `scripts/`

These are the driver scripts as they were run. They contain hard-coded paths
(`/home/claude/gconf/...`) and refer to engine binaries by version (`./oconf4`,
`./oconf22`, and so on). The documented and path-independent versions are in `surveys/`,
`validation/` and `python/`. The logs these scripts produced are in `logs/` and in
`data/nonplanar/earlier/`.

| script | purpose |
|---|---|
| `headline.py`, `verify.py` | single large integral computations and the independent re-checks of the triangle examples (now `surveys/headline.py` and `surveys/verify_case.py`) |
| `g2.py` | first version of the vertex-explosion shortcut for triangle graphs (now inside `surveys/triangle_shortcut.py`) |
| `fat.py` | the "fat" (multiplied-edge) cycle graphs: doubled and tripled C_4, doubled K_{2,3}, with Z/5 quotients |
| `k5sub.py`, `k5run.py`, `k5run2.py`, `k5run3.py`, `k5sum.py`, `k5sum2.py` | K_5 with some edges doubled at n = 5 (now `surveys/k5_partial_doublings.py`) |
| `loc.py`, `loc17.py` | local Smith form runs on the triangle family |
| `bench.py`, `timing.py`, `vexp_bench.py` | timing comparisons between engine versions and between the direct computation and the vertex-explosion route |
| `check_auto.py`, `check_auto_b.py`, `check_twist.py` | the first checks of automorphism quotients and twisted characters (merged into `validation/check_quotients.py`) |
| `colored.py` | coloured (Young subgroup) quotients (now `surveys/colored_quotients.py`) |
| `cutimg.cpp`, `cutimg2.cpp` | earlier versions of `analysis/cutimg3.cpp` |
| `abrams_dt.py` | timing of the Abrams discretised model against the Świątkowski complex |
| `np_queue.py` to `np_queue4.py`, `q17.sh` to `q40.sh`, `run2.sh`, `run3.sh`, `apex_run.sh`, `grp.sh`, `k33c.sh` | the queues of the non-planar campaign. Each queue ran a list of engine invocations one after another, so that only one large computation used the memory at a time. Every run they launched is recorded in `data/nonplanar/campaign.json` |
