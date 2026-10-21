// =====================================================================================================
//  oconf.cpp  --  homology of ordered configuration spaces of graphs
// =====================================================================================================
//
//  WHAT THIS PROGRAM COMPUTES
//  --------------------------
//  Let G be a finite graph (multi-edges and loops allowed) and n >= 1.  The ordered configuration
//  space Conf_n(G) is the space of n-tuples of pairwise distinct points of G.  This program builds a
//  finite chain complex whose homology is H_*(Conf_n(G); Z), namely the REDUCED ORDERED SWIATKOWSKI
//  COMPLEX (Swiatkowski 2001 for the unordered version; the ordered version is written out in
//  Hainaut-Knudsen-Wawrykow, arXiv:2510.00201, Construction 3.11 / Definition 3.12, Theorem 3.14),
//  and then computes ranks and torsion of its boundary maps by sparse Gaussian elimination.
//
//  Optionally it works with the quotient ("coinvariant") complex C_G for a finite group G that acts
//  on the complex by permuting particle labels and/or by graph automorphisms, possibly twisted by a
//  character of G.  If |G| is invertible in the coefficient ring, H_*(C_G) is the G-coinvariants
//  (equivalently the chi-isotypic part) of H_*(Conf_n(G)); see Section 5 of the README.
//
//  THE CELLS
//  ---------
//  Number the edges 0..E-1.  Edge e = (u, v) is oriented u -> v.  It has two half-edges:
//      half-edge 2e   = the end of e at its TAIL u,
//      half-edge 2e+1 = the end of e at its HEAD v.
//  Each vertex v with at least one half-edge has a chosen BASE half-edge h0(v).  By default h0(v) is the
//  lowest-numbered half-edge at v.  When a group of graph automorphisms is used, h0 is chosen so that the
//  group permutes base half-edges among themselves (see choose_invariant_bases in main()).
//
//  A cell is a placement of the n labelled particles 0..n-1 such that each particle is either
//    (a) AT A VERTEX v, carrying the reduced generator (h - h0(v)) for some half-edge h != h0(v) at v;
//        at most one particle per vertex; or
//    (b) IN THE INTERIOR OF AN EDGE e; the particles on e are linearly ordered from tail to head.
//  The DEGREE of a cell is the number of particles sitting at vertices.  The degree is at most
//  min(n, number of vertices of valence >= 2).
//
//  The boundary of a cell c is
//        d(c) = sum_{i} (-1)^i [ c with particle p_i pushed into edge(h_i) at the h_i end
//                              - c with particle p_i pushed into edge(h0) at the h0 end ],
//  where i = 0, 1, 2, ... indexes the OCCUPIED vertices of c in increasing vertex number, p_i is the
//  particle at the i-th occupied vertex, and h_i is its half-edge.  "Pushed in at the tail end" means the
//  particle becomes the FIRST particle of that edge's list; "at the head end" means it becomes the LAST.
//
//  HOW HOMOLOGY IS EXTRACTED
//  -------------------------
//  For each k >= 1 the matrix of d_k is built with ROWS = k-cells and COLUMNS = (k-1)-cells (so it is the
//  transpose of the usual matrix; rank and Smith form do not care).  The elimination engines are:
//
//    SparseElim  exact over Z with unit (+-1) pivots only; whatever cannot be pivoted on (the "residual
//                block") is written to a text file, and the Python driver computes its Smith normal form
//                exactly with FLINT.  Also used for large primes P (-p P with P > 13).
//    F2E         over F_2 (-p 2), entries are just column indices.
//    SPE         over F_p for 3 <= p <= 13 (-p p), packed entries.
//    LPE         LOCAL Smith form over Z/p^K (-q p K): pivots only on units (entries prime to p).  When no
//                unit is left, all remaining entries are divisible by p; they are divided by p and the
//                process continues ("stage" 1, 2, ...).  The number of pivots found in stage s equals the
//                number of elementary divisors of d_k that are exactly p^s (for s < K).  So the
//                p-primary torsion of H_{k-1} is read off as  (Z/p)^{stage1} + (Z/p^2)^{stage2} + ...
//                Elementary divisors divisible by p^K are invisible ("leftover"), see README.
//    LGR         local Smith form over the Galois ring GR(4,2) = Z/4[w]/(w^2+w+1) (-W ...), used for the
//                omega-isotypic part of a Z/3-symmetric complex at the prime 2 (file gr42.h).
//
//  When the sparse elimination runs out of cheap pivots (fill-in explodes) but only a moderate number
//  m of columns is still alive, LPE and LGR switch to a DENSE phase: the remaining rows are streamed
//  through a reduced-row-echelon basis kept densely (one byte per entry for general Z/p^K, or bit-sliced
//  for Z/4, Z/9 and GR(4,2)).  See the comments at LPE::dense_finish.
//
//  An important memory saving ("column dropping"): if d_k was fully reduced using unit pivots, the
//  k-cells that served as pivot rows of d_k can be deleted from the columns of d_{k+1}.  Reason: after a
//  unit-pivot reduction, the kernel of d_k projects isomorphically onto the coordinates of the non-pivot
//  k-cells, so the image of d_{k+1} (which lies in the kernel) is faithfully represented by those
//  coordinates.  The flag canDrop records whether this is allowed for the next degree.
//
//  INPUT (standard input)
//  ----------------------
//      V E n
//      u_0 v_0
//      ...
//      u_{E-1} v_{E-1}
//  Vertices are 0..V-1; edge e goes from u_e to v_e.  Limits: n <= 10, 3E < 255 (so a half-edge or edge
//  location fits in 8 bits of the cell key).
//
//  COMMAND-LINE OPTIONS  (full explanations with examples are in README.md, Section 4)
//  --------------------
//    -p P        coefficients: P = 0 (default) exact over Z; P = 2: F_2; 3 <= P <= 13: F_P; else F_P via
//                SparseElim (P must be prime).
//    -q p K      local Smith form over Z/p^K (p prime, p^K < 32).  Reports "stages".
//    -o PREFIX   prefix for residual files PREFIX_k.txt (exact mode); default "res".
//    -d m        only compute degrees 0..m (maps d_1..d_m), i.e. H_0..H_{m-1} completely.
//    -c          only count cells (no matrices).
//    -z          check d o d = 0 on every cell of degree >= 2 and exit.
//    -T          (with -p 2) rank-only treatment of the top map via its transpose.
//    -k COLORS   quotient by a Young subgroup: COLORS is a word like "aabc" assigning a colour to each
//                particle; particles of equal colour are indistinguishable.
//    -g GENS     label-permutation group: generators separated by ';', each a comma-separated permutation
//                of 0..n-1 (shorter lists are extended by the identity).  Example: -g "1,2,3,4,0".
//    -a GENS     graph automorphisms: generators separated by ';'.  Each is "v_0,...,v_{V-1}" (image of
//                each vertex) optionally followed by "|e_0,...,e_{E-1}" (image of each edge; prefix "~"
//                marks an edge mapped with reversed orientation).  Without the edge part, parallel edges
//                are matched in input order, reversing orientation when needed.
//    -X TW       character twist: TW contains 'l' for the sign of the label permutation and/or 'e' for
//                the sign of the edge permutation.
//    -Y s1,s2,.. a sign (+1/-1) for each automorphism generator, extended to a character (checked).
//    -W t1,t2,.. switch to GR(4,2) mode: exponent t_i in Z/3 for each automorphism generator, i.e. the
//                generator acts by w^{t_i}; the group must have odd order.
//    -F          fast orbit enumeration (strongly recommended whenever -g or -a is used).
//    -D m        LPE dense phase allowed when at most m columns are alive (0 = never; -1 = immediately).
//    -M a b      dense phase: switch immediately once live columns <= a, otherwise only when the sparse
//                part holds more than b nonzero entries.  Defaults a = 30000, b = 120000000.
//    -C c        dense phase only considered once the cheapest pivot column has >= c entries (default 48).
//    -G m        LGR dense-phase threshold (default 70000).
//    -m          generic mode: read a sparse integer matrix "rows cols" + "r c v" triples from stdin and run
//                SparseElim on it (used by the Python driver for auxiliary matrices).
//
//  OUTPUT (standard output, one JSON object)
//  ------------------------------------------
//    {"V":..,"E":..,"n":..,"P":..,"dim":D,"cells":[|C_0|,..,|C_D|],
//     "maps":[{"k":k,"rows":|C_k|,"cols":..,"pivots":rank found,"res_rows":..,"res_cols":..,
//              "res_nnz":..,"maxabs":..,"res_file":"..", ["stages":[..],"leftover_nnz":..]}, ...]}
//  Progress messages go to standard error.
//
//  EXIT CODES
//  ----------
//    1 bad input / graph too large      2 a boundary face was not found among the (k-1)-cells (bug guard)
//    3 integer overflow in exact mode   4 local mode: a non-unit entry was not divisible by p (bug guard)
//    5 dense phase: unsupported modulus 6 dense phase: entry in a dead column (bug guard)
//    7 dense phase: internal invariant violated (bug guard)
//    8 invalid automorphism, or no invariant base half-edge exists
//    9 out of memory in a dense phase   10 the group order is divisible by p (coinvariants invalid)
//   11 group does not factor as labels x automorphisms (bug guard)
//   12 -Y or -W values do not define a character      13 GR mode needs a group of odd order
//
//  BUILD:   g++ -O3 -march=native -fopenmp -std=c++17 -o oconf oconf.cpp      (see Makefile)
// =====================================================================================================

#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

// ----------------------------------------------------------------------------------------------------
//  Global description of the graph and of the run
// ----------------------------------------------------------------------------------------------------
int V, E, N;                    // number of vertices, number of edges, number of particles n
vector<int> etail, ehead;       // etail[e], ehead[e]: endpoints of edge e (oriented tail -> head)
vector<int> hvert;              // hvert[h]: the vertex at which half-edge h sits (h = 2e tail, 2e+1 head)
vector<vector<int>> vhalf;      // vhalf[v]: list of half-edges at v, in increasing order
vector<int> base;               // base[v]: the base half-edge h0(v), or -1 if v has no half-edges

bool checksq = false;           // -z
bool transposeTop = false;      // -T
bool colored = false;           // -k
vector<int> pcolor;             // pcolor[p] = colour (0,1,2,..) of particle p, for -k
bool grouped = false;           // -g or -a given
vector<vector<int>> gelts;      // label-permutation generators (-g)
vector<string> agens;           // automorphism generators exactly as given (-a)
long long LQ = 0;               // -q p K : the prime p (0 if not in local mode)
int LK = 0;                     //          and the exponent K
int twistLab = 0, twistEdge = 0;// -X l / -X e
vector<int> agSigns;            // -Y signs, one per automorphism generator
bool grMode = false;            // -W given: GR(4,2) arithmetic
vector<int> agW;                // -W exponents (mod 3), one per automorphism generator
bool fastEnum = false;          // -F

// ----------------------------------------------------------------------------------------------------
//  Group elements
//
//  Every element of the group G = L x A is stored as a GE ("group element"):
//    lab[p]   image of particle label p           (the L component)
//    vm[v]    image of vertex v                   (the A component, a graph automorphism)
//    em[e]    image of edge e, fl[e] = 1 if the edge is mapped with reversed orientation
//    hm[h]    image of half-edge h (precomputed from em and fl)
//    chi      the value (+1 or -1) of the chosen sign character on this element
//    w        the exponent in Z/3 of the chosen omega character (GR mode only)
// ----------------------------------------------------------------------------------------------------
struct GE {
    vector<int> lab, vm, em;
    vector<char> fl;
    int chi = 1;
    int w = 0;
    int hm[256];
};
vector<GE> gall;                // all NON-identity elements of G (shuffled; see main)
vector<GE> Aelts, Lelts;        // A (automorphisms, identity first) and L (label perms, identity first)
vector<const GE *> gLab, gAut;  // pointers into gall: pure label elements, pure automorphism elements
thread_local int lastW = 0;     // omega exponent of the last canonicalisation (GR mode), per thread

// Sign of a permutation given as an image list p[0..n-1]: product over cycles of (-1)^(length-1).
int permsign(const vector<int> &p) {
    int n = p.size(), s = 1;
    vector<char> seen(n, 0);
    for (int i = 0; i < n; i++)
        if (!seen[i]) {
            int len = 0, j = i;
            while (!seen[j]) {
                seen[j] = 1;
                j = p[j];
                len++;
            }
            if (len % 2 == 0) // a cycle of even length is an odd permutation
                s = -s;
        }
    return s;
}

// ----------------------------------------------------------------------------------------------------
//  Cells and their 128-bit keys
//
//  A Cell is the "readable" form:
//     vocc[v] = particle sitting at vertex v (or -1),  vh[v] = the half-edge h of its generator h - h0(v),
//     elist[e] = particles on edge e listed from tail to head.
//  A Key is a compact, totally ordered encoding used for storage, sorting and lookup.  Particle p owns
//  the 12 bits  [12p, 12p+12)  of the key:
//     low 8 bits  = location: the half-edge h (0 <= h < 2E) if p sits at a vertex,
//                             or 2E + e if p lies on edge e;
//     high 4 bits = rank of p on its edge (0 = closest to the tail); always 0 for vertex particles.
//  This needs 12n <= 128 bits (n <= 10) and 2E + E < 256.
// ----------------------------------------------------------------------------------------------------
struct Cell {
    vector<int> vocc;
    vector<int> vh;
    vector<vector<int>> elist;
};
typedef unsigned __int128 Key;

// Pack a cell into a 128-bit key: 12 bits per particle, the low 8 its location (half-edge if at a vertex, 2E + e if on edge e) and the high 4 its rank in the word of that edge.
Key encode(const Cell &c) {
    Key s = 0;
    for (int v = 0; v < V; v++)
        if (c.vocc[v] >= 0) {
            int p = c.vocc[v];
            s |= (Key)(unsigned)(c.vh[v]) << (12 * p);
        }
    for (int e = 0; e < E; e++)
        for (int r = 0; r < (int)c.elist[e].size(); r++) {
            int p = c.elist[e][r];
            s |= (Key)(unsigned)((2 * E + e) | (r << 8)) << (12 * p);
        }
    return s;
}

// Inverse of encode: rebuild the occupants of the vertices, their half-edges and the edge words from a key.
Cell decode(Key s) {
    Cell c;
    c.vocc.assign(V, -1);
    c.vh.assign(V, -1);
    c.elist.assign(E, {});
    vector<vector<pair<int, int>>> tmp(E); // (rank, particle) per edge, sorted afterwards
    for (int p = 0; p < N; p++) {
        unsigned w = (unsigned)((s >> (12 * p)) & 0xFFF);
        int loc = w & 0xFF, r = w >> 8;
        if (loc < 2 * E) {
            int v = hvert[loc];
            c.vocc[v] = p;
            c.vh[v] = loc;
        } else
            tmp[loc - 2 * E].push_back({r, p});
    }
    for (int e = 0; e < E; e++) {
        sort(tmp[e].begin(), tmp[e].end());
        for (auto &x : tmp[e])
            c.elist[e].push_back(x.second);
    }
    return c;
}

// ----------------------------------------------------------------------------------------------------
//  The action of a group element on a key
//
//  g acts on a cell by: relabelling particle p as lab[p]; moving a vertex particle with half-edge h to
//  half-edge hm[h] (at vertex vm[v]); moving the particles of edge e to edge em[e], reversing their order
//  if fl[e] = 1 (an orientation-reversing image).
//
//  The induced map on chains is a SIGNED permutation of cells.  The sign has two factors:
//    * the Koszul sign: the boundary formula numbers occupied vertices in increasing order, so a cell is
//      really an ordered (exterior) product of its vertex generators.  If g sends the occupied vertices
//      v_1 < ... < v_k to vm[v_1], ..., vm[v_k], re-sorting those images costs the sign of the sorting
//      permutation, i.e. (-1)^(number of inversions);
//    * the character value chi(g), if a twist was requested.
//  A label permutation alone never changes signs (the differential does not look at labels).
//  Returns the image key and stores the sign in 'sign'.
// ----------------------------------------------------------------------------------------------------
static inline Key applyKey(const GE &g, Key s, int &sign) {
    int len[128];  // number of particles on each edge (needed to reverse ranks on flipped edges)
    unsigned w[16];
    for (int e = 0; e < E; e++)
        len[e] = 0;
    for (int p = 0; p < N; p++) {
        w[p] = (unsigned)((s >> (12 * p)) & 0xFFF);
        unsigned loc = w[p] & 0xFF;
        if ((int)loc >= 2 * E)
            len[loc - 2 * E]++;
    }
    Key out = 0;
    int occ[16], k = 0;
    for (int p = 0; p < N; p++) {
        unsigned loc = w[p] & 0xFF, r = w[p] >> 8, nl, nr = 0;
        if ((int)loc < 2 * E) {          // vertex particle: move its half-edge
            nl = g.hm[loc];
            occ[k++] = hvert[loc];
        } else {                         // edge particle: move its edge, maybe reverse its rank
            int e = loc - 2 * E;
            nl = 2 * E + g.em[e];
            nr = g.fl[e] ? len[e] - 1 - r : r;
        }
        out |= (Key)(nl | (nr << 8)) << (12 * g.lab[p]); // particle p becomes particle lab[p]
    }
    sort(occ, occ + k);                  // occupied vertices in increasing order
    int inv = 0;
    for (int i = 0; i < k; i++)
        for (int j = i + 1; j < k; j++)
            if (g.vm[occ[i]] > g.vm[occ[j]])
                inv++;
    sign = ((inv & 1) ? -1 : 1) * g.chi;
    return out;
}

// ----------------------------------------------------------------------------------------------------
//  Canonical orbit representatives
//
//  The basis of the coinvariant complex C_G is one cell per G-orbit; we choose the orbit element with the
//  SMALLEST KEY as representative.  If g maps c to the representative b = g.c with total sign s (Koszul
//  sign times character), then in C_G we have [c] = s * [b]   (for GR mode also a factor w^{w(g)}).
//
//  A cell c that is fixed by some g acting with a nontrivial scalar (sign -1, or w^t with t != 0) is ZERO
//  in C_G after inverting |G|: [c] = lambda [c] with 1 - lambda a unit.  (For lambda = -1 this needs 2
//  invertible, i.e. p odd; for odd-order groups lambda = -1 cannot occur on a fixed cell at all.)
//
//  canonFactored walks over the factorisation G = A x L: it applies each automorphism a once, and then
//  each label permutation l (which only moves 12-bit fields around), which is much cheaper than applying
//  all |A|*|L| elements from scratch.  Return value: s (+1/-1), or 0 if c is zero in C_G.  The omega
//  exponent of the winning element is left in lastW.
// ----------------------------------------------------------------------------------------------------
int canonFactored(Key kc, Key &best) {
    best = kc;
    int bs = 1, bw = 0;
    for (auto &a : Aelts) {
        int sa;
        Key ka = applyKey(a, kc, sa);
        unsigned f[16];
        for (int p = 0; p < N; p++)
            f[p] = (unsigned)((ka >> (12 * p)) & 0xFFF);
        for (auto &l : Lelts) {
            Key k = 0;
            for (int p = 0; p < N; p++)
                k |= (Key)f[p] << (12 * l.lab[p]);
            int sg = sa * l.chi;
            if (k == kc && (sg < 0 || a.w != 0)) {   // c is fixed by a.l with a nontrivial scalar
                if (sg < 0 && a.w == 0 && grMode) {
                    fprintf(stderr, "gr: sign stabiliser (group must be odd)\n");
                    exit(13);
                }
                return 0;
            }
            if (k < best) {
                best = k;
                bs = sg;
                bw = a.w;
            }
        }
    }
    lastW = bw;
    return bs;
}

// Replace c by its canonical representative and return the coefficient s with [c] = s [rep]
// (0 if c vanishes).  Without a group it does nothing (returns 1).  With -k (Young subgroup) it relabels
// the particles of each colour in increasing order of appearance; Young-subgroup quotients have no signs.
int canonicalize(Cell &c) {
    lastW = 0;
    if (grouped) {
        Key kc = encode(c), best;
        int bs = canonFactored(kc, best);
        if (!bs)
            return 0;
        c = decode(best);
        return bs;
    }
    if (!colored)
        return 1;
    map<int, vector<int>> pool; // labels of each colour, increasing
    for (int p = 0; p < N; p++)
        pool[pcolor[p]].push_back(p);
    map<int, int> used;
    auto relabel = [&](int p) {
        int col = pcolor[p];
        return pool[col][used[col]++];
    };
    for (int v = 0; v < V; v++)
        if (c.vocc[v] >= 0)
            c.vocc[v] = relabel(c.vocc[v]);
    for (int e = 0; e < E; e++)
        for (auto &p : c.elist[e])
            p = relabel(p);
    return 1;
}

// Is the key kc its own orbit representative, and nonzero?  (Used by the slow enumerator.)
bool isCanonicalKey(Key kc) {
    for (auto &g : gall) {
        int sg;
        Key k = applyKey(g, kc, sg);
        if (k < kc)
            return false;
        if (k == kc && sg < 0)
            return false;
    }
    return true;
}

// ----------------------------------------------------------------------------------------------------
//  Enumeration of the k-cells
// ----------------------------------------------------------------------------------------------------

// FAST ORBIT ENUMERATION (-F).  A cell = a "shape" (which vertices are occupied and by which half-edge,
// and how many particles lie on each edge) together with a labelling (which particle goes into which
// slot).  Steps:
//   1. enumerate shapes; keep a shape only if it is lexicographically minimal (as a sorted list of
//      12-bit fields) among its images under the automorphism group A;
//   2. for each kept shape run through all n! labellings, keep those that are minimal under the label
//      group L alone (cheap early-abort test);
//   3. canonicalise each survivor under the full group and record it (if nonzero);
//   4. sort and remove duplicates.
// Every G-orbit is reached: some element of it has an A-minimal shape, and within that shape (labels do
// not change shapes) some labelling is L-minimal.  The output is exactly the set of key-minimal nonzero
// orbit representatives, i.e. the same list as the slow enumerator produces.
void enumerate_fast(int k, vector<Key> &out) {
    vector<int> occh(V, -1);    // occh[v] = half-edge used at v, or -1
    vector<int> lens(E, 0);     // number of particles on each edge
    vector<unsigned> fields;    // the n slot codes of the current shape
    fields.reserve(N);
    long long shapes = 0, canonShapes = 0;
    function<void(int, int)> recV, recE;
    auto process = [&]() {
        fields.clear();
        for (int v = 0; v < V; v++)
            if (occh[v] >= 0)
                fields.push_back((unsigned)occh[v]);
        for (int e = 0; e < E; e++)
            for (int r = 0; r < lens[e]; r++)
                fields.push_back((unsigned)((2 * E + e) | (r << 8)));
        shapes++;
        // step 1: A-minimality of the shape
        vector<unsigned> sig(fields);
        sort(sig.begin(), sig.end());
        Key tk = 0;
        for (int i = 0; i < N; i++)
            tk |= (Key)fields[i] << (12 * i);
        for (auto *g : gAut) {
            int sg;
            Key kk = applyKey(*g, tk, sg);
            vector<unsigned> s2(N);
            for (int i = 0; i < N; i++)
                s2[i] = (unsigned)((kk >> (12 * i)) & 0xFFF);
            sort(s2.begin(), s2.end());
            if (s2 < sig)
                return;
        }
        canonShapes++;
        // step 2 and 3: labellings
        vector<int> perm(N);
        iota(perm.begin(), perm.end(), 0);
        do {
            Key kc = 0;
            for (int i = 0; i < N; i++)
                kc |= (Key)fields[i] << (12 * perm[i]);   // slot i receives particle perm[i]
            bool ok = true;
            for (auto *g : gLab) {
                int sg;
                if (applyKey(*g, kc, sg) < kc) {
                    ok = false;
                    break;
                }
            }
            if (!ok)
                continue;
            Key best;
            if (canonFactored(kc, best))
                out.push_back(best);
        } while (next_permutation(perm.begin(), perm.end()));
    };
    recE = [&](int e, int rem) {        // distribute rem particles over edges e..E-1
        if (e == E - 1) {
            lens[e] = rem;
            process();
            lens[e] = 0;
            return;
        }
        for (int l = 0; l <= rem; l++) {
            lens[e] = l;
            recE(e + 1, rem - l);
        }
        lens[e] = 0;
    };
    recV = [&](int v, int left) {       // choose 'left' more occupied vertices among v..V-1
        if (v == V) {
            if (left == 0) {
                if (E > 0)
                    recE(0, N - k);
                else if (N == k)
                    process();
            }
            return;
        }
        if (V - v > left)
            recV(v + 1, left);          // leave v empty
        if (left > 0)
            for (int h : vhalf[v])
                if (h != base[v]) {
                    occh[v] = h;
                    recV(v + 1, left - 1);
                    occh[v] = -1;
                }
    };
    recV(0, k);
    sort(out.begin(), out.end());
    out.erase(unique(out.begin(), out.end()), out.end());
    fprintf(stderr, "  [enum] deg %d: shapes %lld, A-canonical %lld, orbit cells %zu\n", k, shapes,
            canonShapes, out.size());
}

// Enumerate all k-cells (or their orbit representatives).  The slow path places particles 0,1,...,n-1
// one at a time, each either at a free vertex (with any non-base half-edge) or at any position on any
// edge, and keeps cells with exactly k vertex particles.
void enumerate(int k, vector<Key> &out) {
    if (grouped && fastEnum) {
        enumerate_fast(k, out);
        return;
    }
    Cell c;
    c.vocc.assign(V, -1);
    c.vh.assign(V, -1);
    c.elist.assign(E, {});
    int moving = 0;                     // number of particles currently at vertices
    function<void(int)> rec = [&](int p) {
        if (p == N) {
            if (moving == k) {
                if (grouped) {
                    Key kc = encode(c);
                    if (isCanonicalKey(kc))
                        out.push_back(kc);
                } else if (colored) {
                    Cell d = c;
                    canonicalize(d);
                    out.push_back(encode(d));   // duplicates removed by the caller
                } else
                    out.push_back(encode(c));
            }
            return;
        }
        int remaining = N - p;
        if (moving < k) {               // particle p at a vertex
            for (int v = 0; v < V; v++)
                if (c.vocc[v] < 0)
                    for (int h : vhalf[v])
                        if (h != base[v]) {
                            c.vocc[v] = p;
                            c.vh[v] = h;
                            moving++;
                            rec(p + 1);
                            moving--;
                            c.vocc[v] = -1;
                            c.vh[v] = -1;
                        }
        }
        if (k - moving <= remaining - 1) { // particle p on an edge (only if k is still reachable)
            for (int e = 0; e < E; e++) {
                int L = c.elist[e].size();
                for (int pos = 0; pos <= L; pos++) {
                    c.elist[e].insert(c.elist[e].begin() + pos, p);
                    rec(p + 1);
                    c.elist[e].erase(c.elist[e].begin() + pos);
                }
            }
        }
    };
    rec(0);
}

// =====================================================================================================
//  SPARSE ELIMINATION ENGINES
//
//  All engines share the same strategy, a simple form of Markowitz pivoting (Markowitz 1957):
//    * a priority queue holds the live columns keyed by their current number of nonzero entries;
//    * the column with the fewest entries is taken; among its rows that carry an admissible pivot
//      (a unit) the SHORTEST row is chosen; the pivot row is subtracted from every other row having an
//      entry in that column; pivot row and column are then deleted;
//    * colrows[c] lists the rows having an entry in column c.  It is updated lazily: entries for rows that
//      lost their entry in c are not removed immediately, they are filtered when the list becomes more
//      than twice as long as the true count colcnt[c].
//  Each engine stores a row as a sorted vector of (column, value) entries.  The number of pivots is the
//  rank found.  rowdead[r] = 1 marks the rows used as pivot rows (needed for column dropping).
// =====================================================================================================

// ---------------- exact elimination over Z (unit pivots only), or over F_P for a large prime P -------
struct SparseElim {
    int nr, nc;
    vector<vector<pair<int, int>>> rows;    // entries (column, value); values fit in 32 bits
    vector<vector<int>> colrows;
    vector<int> colcnt;
    vector<char> rowdead, coldead;
    ll P;                                   // 0 = integers (pivots must be +-1); else a prime modulus
    long long pivots = 0;
    ll maxabs = 1;                          // largest absolute value seen (exact mode diagnostics)

    SparseElim(int nr_, int nc_, ll P_)
        : nr(nr_), nc(nc_), rows(nr_), colrows(nc_), colcnt(nc_, 0), rowdead(nr_, 0), coldead(nc_, 0), P(P_) {}

    // Merge duplicate entries of each row, reduce mod P, drop zeros, build the column lists.
    void finalize() {
        for (int r = 0; r < nr; r++) {
            auto &R = rows[r];
            sort(R.begin(), R.end());
            vector<pair<int, ll>> M;
            for (auto &x : R) {
                if (!M.empty() && M.back().first == x.first)
                    M.back().second += x.second;
                else
                    M.push_back(x);
            }
            vector<pair<int, int>> M2;
            for (auto &x : M) {
                ll v = x.second;
                if (P) {
                    v %= P;
                    if (v < 0)
                        v += P;
                }
                if (llabs(v) > INT_MAX) {
                    fprintf(stderr, "overflow\n");
                    exit(3);
                }
                if (v != 0)
                    M2.push_back({x.first, (int)v});
            }
            R.swap(M2);
            for (auto &x : R) {
                colrows[x.first].push_back(r);
                colcnt[x.first]++;
            }
        }
    }
    // Is this value usable as a pivot?  Over Z only +-1; over F_P any nonzero value.
    bool isunit(ll v) const {
        if (P)
            return v != 0;
        return v == 1 || v == -1;
    }
    // Entry (r, c) of the sparse matrix, found by binary search in row r; 0 if absent.
    ll getval(int r, int c) const {
        auto &R = rows[r];
        auto it = lower_bound(R.begin(), R.end(), make_pair(c, INT_MIN));
        if (it != R.end() && it->first == c)
            return it->second;
        return 0;
    }
    // a^e mod m by repeated squaring (128-bit intermediate products).
    static ll powmod(ll a, ll e, ll m) {
        ll r = 1;
        a %= m;
        while (e) {
            if (e & 1)
                r = (__int128)r * a % m;
            a = (__int128)a * a % m;
            e >>= 1;
        }
        return r;
    }
    // Pivot on entry (r, c): for every other row r2 with an entry a in column c,
    // row r2 <- row r2 - (a / pivot) * row r.  Overflow is detected and aborts (exit 3).
    void pivot(int r, int c) {
        ll u = getval(r, c);
        ll uinv = P ? powmod(u, P - 2, P) : u;       // over Z, u = +-1 is its own inverse
        vector<pair<int, int>> &R = rows[r];
        vector<int> targets;
        for (int r2 : colrows[c])
            if (r2 != r && !rowdead[r2])
                targets.push_back(r2);
        sort(targets.begin(), targets.end());
        targets.erase(unique(targets.begin(), targets.end()), targets.end());
        for (int r2 : targets) {
            ll a = getval(r2, c);
            if (a == 0)
                continue;
            ll f = P ? (ll)((__int128)a * uinv % P) : a * uinv;
            auto &S = rows[r2];
            vector<pair<int, int>> out;      // merge of S and -f*R, both sorted by column
            out.reserve(S.size() + R.size());
            size_t i = 0, j = 0;
            while (i < S.size() || j < R.size()) {
                if (j == R.size() || (i < S.size() && S[i].first < R[j].first)) {
                    out.push_back(S[i]);
                    i++;
                } else if (i == S.size() || R[j].first < S[i].first) {   // new entry (fill-in)
                    ll v;
                    if (P) {
                        v = (P - (__int128)f * R[j].second % P) % P;
                    } else {
                        if (__builtin_mul_overflow(f, R[j].second, &v)) {
                            fprintf(stderr, "overflow\n");
                            exit(3);
                        }
                        v = -v;
                    }
                    if (llabs(v) > INT_MAX) {
                        fprintf(stderr, "overflow\n");
                        exit(3);
                    }
                    if (v != 0) {
                        out.push_back({R[j].first, (int)v});
                        colrows[R[j].first].push_back(r2);
                        colcnt[R[j].first]++;
                    }
                    j++;
                } else {                                                  // both rows have this column
                    ll v;
                    if (P) {
                        v = ((S[i].second - (__int128)f * R[j].second % P) % P + P) % P;
                    } else {
                        ll t;
                        if (__builtin_mul_overflow(f, R[j].second, &t) ||
                            __builtin_sub_overflow(S[i].second, t, &v)) {
                            fprintf(stderr, "overflow\n");
                            exit(3);
                        }
                    }
                    if (llabs(v) > INT_MAX) {
                        fprintf(stderr, "overflow\n");
                        exit(3);
                    }
                    if (v != 0) {
                        out.push_back({S[i].first, (int)v});
                        if (!P)
                            maxabs = max(maxabs, llabs(v));
                    } else
                        colcnt[S[i].first]--;                             // cancellation
                    i++;
                    j++;
                }
            }
            S.swap(out);
        }
        for (auto &x : R)
            colcnt[x.first]--;
        rowdead[r] = 1;
        coldead[c] = 1;
        colcnt[c] = 0;
        vector<pair<int, int>>().swap(rows[r]);
        vector<int>().swap(colrows[c]);
        pivots++;
    }
    void compact_col(int c) {          // drop stale entries of colrows[c]
        auto &L = colrows[c];
        vector<int> M;
        for (int r : L)
            if (!rowdead[r] && getval(r, c) != 0)
                M.push_back(r);
        sort(M.begin(), M.end());
        M.erase(unique(M.begin(), M.end()), M.end());
        L.swap(M);
    }
    // One pass over the columns in order of increasing count.  Returns true if some pivot was made.
    bool sweep() {
        typedef pair<int, int> PI;
        priority_queue<PI, vector<PI>, greater<PI>> pq;
        for (int c = 0; c < nc; c++)
            if (!coldead[c] && colcnt[c] > 0)
                pq.push({colcnt[c], c});
        bool any = false;
        while (!pq.empty()) {
            auto [cnt, c] = pq.top();
            pq.pop();
            if (coldead[c] || colcnt[c] != cnt || cnt == 0)
                continue;                                     // stale queue entry
            if (colrows[c].size() > 2 * (size_t)colcnt[c] + 8)
                compact_col(c);
            int best = -1;
            size_t bestlen = SIZE_MAX;
            for (int r : colrows[c]) {
                if (rowdead[r])
                    continue;
                ll v = getval(r, c);
                if (v != 0 && isunit(v) && rows[r].size() < bestlen) {
                    bestlen = rows[r].size();
                    best = r;
                }
            }
            if (best < 0)
                continue;                                     // no unit in this column (yet)
            vector<int> touched;
            for (auto &x : rows[best])
                touched.push_back(x.first);
            pivot(best, c);
            any = true;
            for (int c2 : touched)                            // their counts changed: re-queue them
                if (!coldead[c2] && colcnt[c2] > 0)
                    pq.push({colcnt[c2], c2});
        }
        return any;
    }
    // Elimination driver: repeatedly pivot in the column with the fewest usable entries (Markowitz order) on its shortest row, until no usable pivot is left.
    void run() {
        while (sweep()) {
        }
    }
    // Write the remaining (non-unit) block as "rows cols" followed by "r c v" lines, renumbered densely.
    void dump_residual(const string &fn, int &rr, int &rc, long long &rnnz) {
        vector<int> rmap(nr, -1), cmap(nc, -1);
        rr = 0;
        rc = 0;
        rnnz = 0;
        for (int r = 0; r < nr; r++)
            if (!rowdead[r] && !rows[r].empty())
                rmap[r] = rr++;
        for (int c = 0; c < nc; c++)
            if (!coldead[c] && colcnt[c] > 0)
                cmap[c] = rc++;
        FILE *f = fopen(fn.c_str(), "w");
        fprintf(f, "%d %d\n", rr, rc);
        for (int r = 0; r < nr; r++)
            if (rmap[r] >= 0)
                for (auto &x : rows[r]) {
                    fprintf(f, "%d %d %lld\n", rmap[r], cmap[x.first], (long long)x.second);
                    rnnz++;
                }
        fclose(f);
    }
};

// ---------------- elimination over F_2: a row is the sorted list of its nonzero columns --------------
struct F2E {
    int nr, nc;
    vector<vector<int>> rows;
    vector<vector<int>> colrows;
    vector<int> colcnt;
    vector<char> rowdead, coldead;
    long long pivots = 0;
    ll maxabs = 1;
    F2E(int r, int c) : nr(r), nc(c), rows(r), colrows(c), colcnt(c, 0), rowdead(r, 0), coldead(c, 0) {}
    void finalize() {                  // entries occurring an even number of times cancel mod 2
        for (int r = 0; r < nr; r++) {
            auto &R = rows[r];
            sort(R.begin(), R.end());
            vector<int> M;
            for (size_t i = 0; i < R.size();) {
                size_t j = i;
                while (j < R.size() && R[j] == R[i])
                    j++;
                if ((j - i) & 1)
                    M.push_back(R[i]);
                i = j;
            }
            R.swap(M);
            for (int c : R) {
                colrows[c].push_back(r);
                colcnt[c]++;
            }
        }
    }
    // Over F_2 a row is the set of its nonzero columns: is (r, c) nonzero?
    bool has(int r, int c) const {
        auto &R = rows[r];
        return binary_search(R.begin(), R.end(), c);
    }
    void pivot(int r, int c) {         // row r2 <- row r2 XOR row r, for all r2 containing c
        vector<int> &R = rows[r];
        vector<int> targets;
        for (int r2 : colrows[c])
            if (r2 != r && !rowdead[r2])
                targets.push_back(r2);
        sort(targets.begin(), targets.end());
        targets.erase(unique(targets.begin(), targets.end()), targets.end());
        for (int r2 : targets) {
            if (!has(r2, c))
                continue;
            auto &S = rows[r2];
            vector<int> out;
            out.reserve(S.size() + R.size());
            size_t i = 0, j = 0;
            while (i < S.size() || j < R.size()) {
                if (j == R.size() || (i < S.size() && S[i] < R[j])) {
                    out.push_back(S[i]);
                    i++;
                } else if (i == S.size() || R[j] < S[i]) {
                    out.push_back(R[j]);
                    colrows[R[j]].push_back(r2);
                    colcnt[R[j]]++;
                    j++;
                } else {
                    colcnt[S[i]]--;
                    i++;
                    j++;
                }
            }
            S.swap(out);
        }
        for (int x : R)
            colcnt[x]--;
        rowdead[r] = 1;
        coldead[c] = 1;
        colcnt[c] = 0;
        pivots++;
        vector<int>().swap(rows[r]);
        vector<int>().swap(colrows[c]);
    }
    // Elimination driver: repeatedly pivot in the column with the fewest usable entries (Markowitz order) on its shortest row, until no usable pivot is left.
    void run() {
        bool any = true;
        while (any) {
            any = false;
            typedef pair<int, int> PI;
            priority_queue<PI, vector<PI>, greater<PI>> pq;
            for (int c = 0; c < nc; c++)
                if (!coldead[c] && colcnt[c] > 0)
                    pq.push({colcnt[c], c});
            while (!pq.empty()) {
                auto [cnt, c] = pq.top();
                pq.pop();
                if (coldead[c] || colcnt[c] != cnt || cnt == 0)
                    continue;
                if (colrows[c].size() > 2 * (size_t)colcnt[c] + 8) {
                    vector<int> M;
                    for (int r : colrows[c])
                        if (!rowdead[r] && has(r, c))
                            M.push_back(r);
                    sort(M.begin(), M.end());
                    M.erase(unique(M.begin(), M.end()), M.end());
                    colrows[c].swap(M);
                }
                int best = -1;
                size_t bl = SIZE_MAX;
                for (int r : colrows[c])
                    if (!rowdead[r] && rows[r].size() < bl && has(r, c)) {
                        bl = rows[r].size();
                        best = r;
                    }
                if (best < 0)
                    continue;
                vector<int> touched(rows[best].begin(), rows[best].end());
                pivot(best, c);
                any = true;
                for (int c2 : touched)
                    if (!coldead[c2] && colcnt[c2] > 0)
                        pq.push({colcnt[c2], c2});
            }
        }
    }
};

// ---------------- elimination over F_p, 3 <= p <= 13: entry packed as (column << 4) | value -----------
struct SPE {
    int nr, nc;
    uint32_t p;
    vector<vector<uint32_t>> rows;
    vector<vector<int>> colrows;
    vector<int> colcnt;
    vector<char> rowdead, coldead;
    long long pivots = 0;
    ll maxabs = 1;
    SPE(int r, int c, uint32_t p_) : nr(r), nc(c), p(p_), rows(r), colrows(c), colcnt(c, 0), rowdead(r, 0), coldead(c, 0) {}
    // Packed entry: the column index (high bits).
    static uint32_t col(uint32_t x) {
        return x >> 4;
    }
    // Packed entry: the value (low bits).
    static uint32_t val(uint32_t x) {
        return x & 15;
    }
    uint32_t inv(uint32_t a) const {   // brute-force inverse (p is tiny)
        for (uint32_t b = 1; b < p; b++)
            if (a * b % p == 1)
                return b;
        return 0;
    }
    // Sort every sparse row, merge repeated columns, drop zero entries, and build the column structures (rows per column, number of usable pivot entries per column).
    void finalize() {
        for (int r = 0; r < nr; r++) {
            auto &R = rows[r];
            sort(R.begin(), R.end());
            vector<uint32_t> M;
            for (size_t i = 0; i < R.size();) {
                uint32_t c = col(R[i]);
                uint32_t v = 0;
                while (i < R.size() && col(R[i]) == c) {
                    v = (v + val(R[i])) % p;
                    i++;
                }
                if (v)
                    M.push_back(c << 4 | v);
            }
            R.swap(M);
            for (auto x : R) {
                colrows[col(x)].push_back(r);
                colcnt[col(x)]++;
            }
        }
    }
    // Entry (r, c) of the sparse matrix, found by binary search in row r; 0 if absent.
    uint32_t get(int r, int c) const {
        auto &R = rows[r];
        auto it = lower_bound(R.begin(), R.end(), (uint32_t)c << 4);
        return (it != R.end() && col(*it) == (uint32_t)c) ? val(*it) : 0;
    }
    // Pivot on the usable entry (r, c): subtract the right multiple of row r from every other live row with an entry in column c, then mark row r and column c as dead.
    void pivot(int r, int c) {
        auto &R = rows[r];
        uint32_t ui = inv(get(r, c));
        vector<int> tg;
        for (int r2 : colrows[c])
            if (r2 != r && !rowdead[r2])
                tg.push_back(r2);
        sort(tg.begin(), tg.end());
        tg.erase(unique(tg.begin(), tg.end()), tg.end());
        for (int r2 : tg) {
            uint32_t a = get(r2, c);
            if (!a)
                continue;
            uint32_t f = a * ui % p;
            auto &S = rows[r2];
            vector<uint32_t> out;
            out.reserve(S.size() + R.size());
            size_t i = 0, j = 0;
            while (i < S.size() || j < R.size()) {
                if (j == R.size() || (i < S.size() && col(S[i]) < col(R[j]))) {
                    out.push_back(S[i]);
                    i++;
                } else if (i == S.size() || col(R[j]) < col(S[i])) {
                    uint32_t v = (p - f * val(R[j]) % p) % p;
                    if (v) {
                        out.push_back(col(R[j]) << 4 | v);
                        colrows[col(R[j])].push_back(r2);
                        colcnt[col(R[j])]++;
                    }
                    j++;
                } else {
                    uint32_t v = (val(S[i]) + p - f * val(R[j]) % p) % p;
                    if (v)
                        out.push_back(col(S[i]) << 4 | v);
                    else
                        colcnt[col(S[i])]--;
                    i++;
                    j++;
                }
            }
            S.swap(out);
        }
        for (auto x : R)
            colcnt[col(x)]--;
        rowdead[r] = 1;
        coldead[c] = 1;
        colcnt[c] = 0;
        pivots++;
        vector<uint32_t>().swap(rows[r]);
        vector<int>().swap(colrows[c]);
    }
    // Elimination driver: repeatedly pivot in the column with the fewest usable entries (Markowitz order) on its shortest row, until no usable pivot is left.
    void run() {
        bool any = true;
        while (any) {
            any = false;
            typedef pair<int, int> PI;
            priority_queue<PI, vector<PI>, greater<PI>> pq;
            for (int c = 0; c < nc; c++)
                if (!coldead[c] && colcnt[c] > 0)
                    pq.push({colcnt[c], c});
            while (!pq.empty()) {
                auto [cnt, c] = pq.top();
                pq.pop();
                if (coldead[c] || colcnt[c] != cnt || cnt == 0)
                    continue;
                if (colrows[c].size() > 2 * (size_t)colcnt[c] + 8) {
                    vector<int> M;
                    for (int r : colrows[c])
                        if (!rowdead[r] && get(r, c))
                            M.push_back(r);
                    sort(M.begin(), M.end());
                    M.erase(unique(M.begin(), M.end()), M.end());
                    colrows[c].swap(M);
                }
                int best = -1;
                size_t bl = SIZE_MAX;
                for (int r : colrows[c])
                    if (!rowdead[r] && rows[r].size() < bl && get(r, c)) {
                        bl = rows[r].size();
                        best = r;
                    }
                if (best < 0)
                    continue;
                vector<int> touched;
                for (auto x : rows[best])
                    touched.push_back(col(x));
                pivot(best, c);
                any = true;
                for (int c2 : touched)
                    if (!coldead[c2] && colcnt[c2] > 0)
                        pq.push({colcnt[c2], c2});
            }
        }
    }
};

// =====================================================================================================
//  LPE: LOCAL SMITH FORM OVER Z/p^K  (option -q p K)
//
//  Entries are integers mod p^K (< 32), packed as (column << 5) | value.  We pivot only on UNITS, i.e.
//  entries not divisible by p.  When a full sweep finds no unit, every live entry is divisible by p
//  (otherwise it would be a unit); we divide all live entries by p, replace the modulus p^K by p^(K-1),
//  and start a new STAGE.  The number of pivots of stage s is the number of elementary divisors of the
//  matrix that are exactly p^s (s = 0 counts units).  Elementary divisors divisible by p^K survive as
//  nonzero entries at the end ("leftover_nonzero") only if they are not 0 mod p^K, so a run mod p^K sees
//  every p-primary cyclic summand of order < p^K exactly; larger ones show up as missing rank compared
//  with a rational (large prime) run.  This is a standard local version of Smith normal form; compare
//  Dumas-Saunders-Villard (2001) on elimination-based Smith form computations.
//
//  Dense finishing phase.  Sparse elimination is cheap until the last few percent of pivots, where fill
//  explodes.  When at most denseMax columns remain alive and either (live columns <= denseMin) or (the
//  sparse part holds more than nnzBudget entries), and the cheapest column already has >= cntThresh
//  entries, the remaining rows are handed to dense_finish().  The dense phase streams the rows one at a
//  time through a dense REDUCED row-echelon basis B (each basis row has a unit 1 in its pivot column and
//  0 in every other pivot column):
//     * reduce the incoming row x against B: for each of x's original entries sitting in a pivot column c,
//       subtract x[c] * B[c]  (because B is reduced, the value x[c] is not changed by the other
//       subtractions, so only the original sparse entries need to be looked at);
//     * if x now has a unit entry in column j: scale x so that x[j] = 1, eliminate column j from every
//       existing basis row (keeps B reduced), append x to B;
//     * otherwise x is "left over": all its entries are divisible by p.
//  At the end of the stage the leftovers are reduced against the final basis, divided by p, and the next
//  stage starts on them.  Three implementations exist:
//     - generic: one byte per entry, any modulus < 32 (axpyD / scalD below);
//     - dense_z4: Z/4 bit-sliced (two 64-bit planes per row: low bit, high bit);
//     - dense_z9: Z/9 bit-sliced (four planes, values 0..8 stored in binary).
//  In z4 and z9 the leftovers are not stored row by row.  Instead y = x/p mod p is inserted into an
//  online F_p row-echelon basis Y (only the F_p-span of the leftovers matters, because the final stage-1
//  rows are an F_p-linear function of y: final(y) = y - sum_{pivot cols c} y[c] * (B[c] mod p)).  This
//  keeps memory bounded by the number of columns.  These two only resolve one further stage (K = 2):
//  they report the number of Z/p summands; Z/p^2 and higher are not resolved.
// =====================================================================================================
struct LPE {
    int nr, nc;
    uint32_t p0, mod;                       // p and the current modulus p^(K - stage)
    vector<vector<uint32_t>> rows;
    vector<vector<int>> colrows;
    vector<int> colcnt;
    vector<char> rowdead, coldead;
    long long pivots = 0, lastCheck = 0;
    bool goDense = false;
    static inline long long denseMax = 0;   // -D
    static inline int cntThresh = 48;       // -C
    static inline long long denseMin = 30000, nnzBudget = 120000000;   // -M
    ll maxabs = 1;
    vector<long long> stages;               // pivots found in stage 0, 1, 2, ...
    long long leftover_nonzero = 0;         // entries that survive all stages (divisors >= p^K)
    LPE(int r, int c, uint32_t p, uint32_t m)
        : nr(r), nc(c), p0(p), mod(m), rows(r), colrows(c), colcnt(c, 0), rowdead(r, 0), coldead(c, 0) {}
    static uint32_t col(uint32_t x) {
        return x >> 5;
    }
    // Packed entry: the value (low bits).
    static uint32_t val(uint32_t x) {
        return x & 31;
    }
    uint32_t inv(uint32_t a) const {        // inverse of a unit mod 'mod' by brute force
        for (uint32_t b = 1; b < mod; b++)
            if (a * b % mod == 1)
                return b;
        return 0;
    }
    // Sort every sparse row, merge repeated columns, drop zero entries, and build the column structures (rows per column, number of usable pivot entries per column).
    void finalize() {
        for (int r = 0; r < nr; r++) {
            auto &R = rows[r];
            sort(R.begin(), R.end());
            vector<uint32_t> M;
            for (size_t i = 0; i < R.size();) {
                uint32_t c = col(R[i]);
                uint32_t v = 0;
                while (i < R.size() && col(R[i]) == c) {
                    v = (v + val(R[i])) % mod;
                    i++;
                }
                if (v)
                    M.push_back(c << 5 | v);
            }
            R.swap(M);
            for (auto x : R) {
                colrows[col(x)].push_back(r);
                colcnt[col(x)]++;
            }
        }
    }
    // Entry (r, c) of the sparse matrix, found by binary search in row r; 0 if absent.
    uint32_t get(int r, int c) const {
        auto &R = rows[r];
        auto it = lower_bound(R.begin(), R.end(), (uint32_t)c << 5);
        return (it != R.end() && col(*it) == (uint32_t)c) ? val(*it) : 0;
    }
    // Is this value usable as a pivot (a unit of the coefficient ring)?
    bool unit(uint32_t v) const {
        return v % p0 != 0;
    }
    // Pivot on the usable entry (r, c): subtract the right multiple of row r from every other live row with an entry in column c, then mark row r and column c as dead.
    void pivot(int r, int c) {
        auto &R = rows[r];
        uint32_t ui = inv(get(r, c));
        vector<int> tg;
        for (int r2 : colrows[c])
            if (r2 != r && !rowdead[r2])
                tg.push_back(r2);
        sort(tg.begin(), tg.end());
        tg.erase(unique(tg.begin(), tg.end()), tg.end());
        for (int r2 : tg) {
            uint32_t a = get(r2, c);
            if (!a)
                continue;
            uint32_t f = a * ui % mod;
            auto &S = rows[r2];
            vector<uint32_t> out;
            out.reserve(S.size() + R.size());
            size_t i = 0, j = 0;
            while (i < S.size() || j < R.size()) {
                if (j == R.size() || (i < S.size() && col(S[i]) < col(R[j]))) {
                    out.push_back(S[i]);
                    i++;
                } else if (i == S.size() || col(R[j]) < col(S[i])) {
                    uint32_t v = (mod - f * val(R[j]) % mod) % mod;
                    if (v) {
                        out.push_back(col(R[j]) << 5 | v);
                        colrows[col(R[j])].push_back(r2);
                        colcnt[col(R[j])]++;
                    }
                    j++;
                } else {
                    uint32_t v = (val(S[i]) + mod - f * val(R[j]) % mod) % mod;
                    if (v)
                        out.push_back(col(S[i]) << 5 | v);
                    else
                        colcnt[col(S[i])]--;
                    i++;
                    j++;
                }
            }
            S.assign(out.begin(), out.end());   // copy (not swap) so the capacity is exact: saves memory
        }
        for (auto x : R)
            colcnt[col(x)]--;
        rowdead[r] = 1;
        coldead[c] = 1;
        colcnt[c] = 0;
        pivots++;
        vector<uint32_t>().swap(rows[r]);
        vector<int>().swap(colrows[c]);
    }
    // One Markowitz pass.  Also: (a) periodically rebuilds the priority queue to drop stale entries,
    // (b) every 200000 pops prints diagnostics and compacts the column lists if they hold too many stale
    // rows, (c) decides when to switch to the dense phase (sets goDense and returns).
    bool sweep() {
        bool any = false;
        typedef pair<int, int> PI;
        priority_queue<PI, vector<PI>, greater<PI>> pq;
        for (int c = 0; c < nc; c++)
            if (!coldead[c] && colcnt[c] > 0)
                pq.push({colcnt[c], c});
        long long sinceDiag = 0;
        while (!pq.empty()) {
            auto [cnt, c] = pq.top();
            pq.pop();
            if (pq.size() > 4 * (size_t)nc + 1000000) {     // (a)
                vector<PI> keep;
                keep.reserve(nc);
                while (!pq.empty()) {
                    auto t = pq.top();
                    pq.pop();
                    if (!coldead[t.second] && colcnt[t.second] == t.first && t.first > 0)
                        keep.push_back(t);
                }
                sort(keep.begin(), keep.end());
                keep.erase(unique(keep.begin(), keep.end()), keep.end());
                pq = priority_queue<PI, vector<PI>, greater<PI>>(keep.begin(), keep.end());
            }
            if (++sinceDiag >= 200000) {                    // (b)
                sinceDiag = 0;
                long long nnz = 0, cre = 0, lc = 0;
                for (int r2 = 0; r2 < nr; r2++)
                    if (!rowdead[r2])
                        nnz += rows[r2].size();
                for (int c2 = 0; c2 < nc; c2++) {
                    cre += colrows[c2].size();
                    if (!coldead[c2] && colcnt[c2] > 0)
                        lc++;
                }
                fprintf(stderr, "  [elim] pivots=%lld rowNNZ=%lld colrowsEntries=%lld liveCols=%lld pq=%zu cnt=%d\n",
                        pivots, nnz, cre, lc, pq.size(), cnt);
                if (cre > nnz + nnz / 2 + 1000000) {
                    for (int c2 = 0; c2 < nc; c2++)
                        if (!coldead[c2] && colrows[c2].size() > (size_t)colcnt[c2]) {
                            vector<int> M;
                            M.reserve(colcnt[c2]);
                            for (int r : colrows[c2])
                                if (!rowdead[r] && get(r, c2))
                                    M.push_back(r);
                            sort(M.begin(), M.end());
                            M.erase(unique(M.begin(), M.end()), M.end());
                            colrows[c2].assign(M.begin(), M.end());
                            colcnt[c2] = M.size();
                        }
                }
            }
            if (coldead[c] || colcnt[c] != cnt || cnt == 0)
                continue;
            if (colrows[c].size() > 2 * (size_t)colcnt[c] + 8) {
                vector<int> M;
                for (int r : colrows[c])
                    if (!rowdead[r] && get(r, c))
                        M.push_back(r);
                sort(M.begin(), M.end());
                M.erase(unique(M.begin(), M.end()), M.end());
                colrows[c].swap(M);
            }
            int best = -1;
            size_t bl = SIZE_MAX;
            for (int r : colrows[c])
                if (!rowdead[r] && rows[r].size() < bl) {
                    uint32_t v = get(r, c);
                    if (v && unit(v)) {
                        bl = rows[r].size();
                        best = r;
                    }
                }
            if (best < 0)
                continue;
            vector<int> touched;
            for (auto x : rows[best])
                touched.push_back(col(x));
            pivot(best, c);
            any = true;
            if (denseMax > 0 && cnt >= cntThresh && pivots - lastCheck >= 2000) {   // (c)
                lastCheck = pivots;
                long long lc = 0;
                for (int c2 = 0; c2 < nc; c2++)
                    if (!coldead[c2] && colcnt[c2] > 0)
                        lc++;
                long long nz = 0;
                if (lc <= denseMax && lc > denseMin)
                    for (int r2 = 0; r2 < nr; r2++)
                        if (!rowdead[r2])
                            nz += rows[r2].size();
                if (lc <= denseMax && (lc <= denseMin || nz > nnzBudget)) {
                    goDense = true;
                    fprintf(stderr, "  [dense] switching: pivots=%lld liveCols=%lld colcnt=%d\n", pivots, lc, cnt);
                    return true;
                }
            }
            for (int c2 : touched)
                if (!coldead[c2] && colcnt[c2] > 0)
                    pq.push({colcnt[c2], c2});
        }
        return any;
    }

    // ---------- generic dense arithmetic: one byte per entry, modulus known at compile time ----------
    template <uint32_t MD> static void axpy(uint8_t *__restrict x, const uint8_t *__restrict b, uint32_t f, int m) {
        uint32_t g = (MD - f) % MD;                 // x <- x - f*b = x + (MD - f)*b  (mod MD)
        for (int j = 0; j < m; j++) {
            uint16_t t = (uint16_t)x[j] + (uint16_t)g * (uint16_t)b[j];
            x[j] = (uint8_t)(t % MD);
        }
    }
    template <uint32_t MD> static void scal(uint8_t *x, uint32_t f, int m) {   // x <- f*x (mod MD)
        for (int j = 0; j < m; j++)
            x[j] = (uint8_t)(((uint16_t)x[j] * f) % MD);
    }
    // Run-time dispatch to a compile-time modulus (lets the compiler replace % by multiplications).
    void axpyD(uint8_t *x, const uint8_t *b, uint32_t f, int m) {
        switch (mod) {
        case 2: axpy<2>(x, b, f, m); break;
        case 4: axpy<4>(x, b, f, m); break;
        case 8: axpy<8>(x, b, f, m); break;
        case 16: axpy<16>(x, b, f, m); break;
        case 3: axpy<3>(x, b, f, m); break;
        case 9: axpy<9>(x, b, f, m); break;
        case 27: axpy<27>(x, b, f, m); break;
        case 5: axpy<5>(x, b, f, m); break;
        case 25: axpy<25>(x, b, f, m); break;
        case 7: axpy<7>(x, b, f, m); break;
        case 11: axpy<11>(x, b, f, m); break;
        case 13: axpy<13>(x, b, f, m); break;
        case 17: axpy<17>(x, b, f, m); break;
        case 19: axpy<19>(x, b, f, m); break;
        case 23: axpy<23>(x, b, f, m); break;
        case 29: axpy<29>(x, b, f, m); break;
        case 31: axpy<31>(x, b, f, m); break;
        default: fprintf(stderr, "dense: unsupported modulus %u\n", mod); exit(5);
        }
    }
    // Dense phase: multiply a row of bytes by the scalar f modulo the current modulus.
    void scalD(uint8_t *x, uint32_t f, int m) {
        switch (mod) {
        case 2: scal<2>(x, f, m); break;
        case 4: scal<4>(x, f, m); break;
        case 8: scal<8>(x, f, m); break;
        case 16: scal<16>(x, f, m); break;
        case 3: scal<3>(x, f, m); break;
        case 9: scal<9>(x, f, m); break;
        case 27: scal<27>(x, f, m); break;
        case 5: scal<5>(x, f, m); break;
        case 25: scal<25>(x, f, m); break;
        case 7: scal<7>(x, f, m); break;
        case 11: scal<11>(x, f, m); break;
        case 13: scal<13>(x, f, m); break;
        case 17: scal<17>(x, f, m); break;
        case 19: scal<19>(x, f, m); break;
        case 23: scal<23>(x, f, m); break;
        case 29: scal<29>(x, f, m); break;
        case 31: scal<31>(x, f, m); break;
        default: exit(5);
        }
    }

    // ---------- Z/4 bit-sliced dense phase -------------------------------------------------------------
    // An element x of Z/4 is stored as two bits x = x0 + 2 x1, in two separate 64-bit "planes", so one
    // machine word holds 64 entries.  Addition: (a0,a1) + (c0,c1) = (a0^c0, a1^c1^(a0&c0)).
    // Negation: -(c0,c1) = (c0, c1^c0).  2*(c0,c1) = (0, c0).
    static inline int getv(const uint64_t *lo, const uint64_t *hi, int c) {
        return (int)((lo[c >> 6] >> (c & 63)) & 1) | (int)(((hi[c >> 6] >> (c & 63)) & 1) << 1);
    }
    // x <- x - f*b over Z/4 (f = 1, 2, 3)
    static inline void z4axpy(uint64_t *__restrict xl, uint64_t *__restrict xh, const uint64_t *__restrict bl,
                              const uint64_t *__restrict bh, int f, int W) {
        if (f == 2) {               // -2b = 2b = (0, b0)
            for (int w = 0; w < W; w++)
                xh[w] ^= bl[w];
            return;
        }
        if (f == 1) {               // add -b = (b0, b1 ^ b0)
            for (int w = 0; w < W; w++) {
                uint64_t a0 = xl[w], a1 = xh[w], c0 = bl[w], c1 = bh[w] ^ bl[w];
                xl[w] = a0 ^ c0;
                xh[w] = a1 ^ c1 ^ (a0 & c0);
            }
            return;
        }
        for (int w = 0; w < W; w++) {   // f == 3: -3b = b, add b
            uint64_t a0 = xl[w], a1 = xh[w], c0 = bl[w], c1 = bh[w];
            xl[w] = a0 ^ c0;
            xh[w] = a1 ^ c1 ^ (a0 & c0);
        }
    }
    // Dense finishing phase over Z/4: bit-sliced rows (two planes), reduced row echelon basis built online; leftovers (all entries even) are halved to F_2.
    void dense_z4(long long sparseStagePivots) {
        vector<int> cmap(nc, -1);       // live column -> dense index 0..m-1
        int m = 0;
        for (int c = 0; c < nc; c++)
            if (!coldead[c] && colcnt[c] > 0)
                cmap[c] = m++;
        vector<vector<int>>().swap(colrows);   // not needed any more: free it
        int W = (m + 63) / 64;
        fprintf(stderr, "  [z4] m=%d, basis memory <= %.2f GB\n", m, 2.0 * m * W * 8 / 1e9);
        vector<uint64_t *> BL, BH;      // basis rows: low plane, high plane
        vector<int> rowOfPiv(m, -1);    // rowOfPiv[c] = index of the basis row whose pivot is column c
        vector<uint64_t> xl(W), xh(W);  // the row being processed
        vector<vector<uint64_t>> Yl;    // online F2 RREF of halved leftovers y = x/2 mod 2
        vector<int> yp;                 // pivot column of each Yl row
        auto insY = [&](vector<uint64_t> &y) {
            for (size_t i = 0; i < Yl.size(); i++) {
                int c = yp[i];
                if ((y[c >> 6] >> (c & 63)) & 1)
                    for (int w = 0; w < W; w++)
                        y[w] ^= Yl[i][w];
            }
            int c = -1;
            for (int w = 0; w < W && c < 0; w++)
                if (y[w])
                    c = w * 64 + __builtin_ctzll(y[w]);
            if (c < 0)
                return;                 // y already in the span
            for (size_t i = 0; i < Yl.size(); i++)
                if ((Yl[i][c >> 6] >> (c & 63)) & 1)
                    for (int w = 0; w < W; w++)
                        Yl[i][w] ^= y[w];
            Yl.push_back(y);
            yp.push_back(c);
        };
        long long done = 0, stagePiv = sparseStagePivots;
        for (int r = 0; r < nr; r++) {
            if (rowdead[r] || rows[r].empty())
                continue;
            fill(xl.begin(), xl.end(), 0);
            fill(xh.begin(), xh.end(), 0);
            for (auto e : rows[r]) {    // densify
                int c = cmap[col(e)];
                if (c < 0) {
                    fprintf(stderr, "z4: dead col\n");
                    exit(6);
                }
                int v = val(e) & 3;
                if (v & 1)
                    xl[c >> 6] |= 1ULL << (c & 63);
                if (v & 2)
                    xh[c >> 6] |= 1ULL << (c & 63);
            }
            for (auto e : rows[r]) {    // reduce against the basis (only original entries matter)
                int c = cmap[col(e)];
                int br = rowOfPiv[c];
                if (br >= 0) {
                    int f = getv(xl.data(), xh.data(), c);
                    if (f)
                        z4axpy(xl.data(), xh.data(), BL[br], BH[br], f, W);
                }
            }
            vector<uint32_t>().swap(rows[r]);  // free the sparse row as soon as it is consumed
            int j = -1;                 // first unit (odd) entry
            for (int w = 0; w < W; w++)
                if (xl[w]) {
                    j = w * 64 + __builtin_ctzll(xl[w]);
                    break;
                }
            if (j >= 0) {
                if (getv(xl.data(), xh.data(), j) == 3)   // normalise the pivot to 1: multiply by 3 = -1
                    for (int w = 0; w < W; w++)
                        xh[w] ^= xl[w];
                long long nb = BL.size();
                const uint64_t *pl = xl.data();
                const uint64_t *ph = xh.data();
#pragma omp parallel for schedule(static, 64)
                for (long long b = 0; b < nb; b++) {        // keep the basis reduced
                    int f = getv(BL[b], BH[b], j);
                    if (f)
                        z4axpy(BL[b], BH[b], pl, ph, f, W);
                }
                uint64_t *nl = (uint64_t *)malloc(W * 8);
                uint64_t *nh = (uint64_t *)malloc(W * 8);
                if (!nl || !nh) {
                    fprintf(stderr, "z4: out of memory\n");
                    exit(9);
                }
                memcpy(nl, xl.data(), W * 8);
                memcpy(nh, xh.data(), W * 8);
                rowOfPiv[j] = BL.size();
                BL.push_back(nl);
                BH.push_back(nh);
                stagePiv++;
                pivots++;
                rowdead[r] = 1;
            } else {                    // leftover: only even entries; keep x/2 mod 2 = the high plane
                bool nz = false;
                for (int w = 0; w < W; w++)
                    if (xh[w]) {
                        nz = true;
                        break;
                    }
                if (nz) {
                    vector<uint64_t> y(xh.begin(), xh.end());
                    insY(y);
                }
            }
            if (++done % 500000 == 0)
                fprintf(stderr, "  [z4] %lld rows, basis %zu, leftover span %zu\n", done, BL.size(), Yl.size());
        }
        stages.push_back(stagePiv);
        fprintf(stderr, "  [z4] stage 0 done: basis %zu, leftover span %zu\n", BL.size(), Yl.size());
        // Final reduction: x = 2y, subtracting x[c]*B[c] = 2*B[c] changes y by (B[c] mod 2) = low plane.
        vector<vector<uint64_t>> F;
        for (auto &y : Yl) {
            vector<int> hits;
            for (int w = 0; w < W; w++) {
                uint64_t t = y[w];
                while (t) {
                    int c = w * 64 + __builtin_ctzll(t);
                    t &= t - 1;
                    if (rowOfPiv[c] >= 0)
                        hits.push_back(c);
                }
            }
            for (int c : hits) {
                const uint64_t *bl = BL[rowOfPiv[c]];
                for (int w = 0; w < W; w++)
                    y[w] ^= bl[w];
            }
            bool nz = false;
            for (int w = 0; w < W; w++) {
                uint64_t t = y[w];
                while (t) {
                    int c = w * 64 + __builtin_ctzll(t);
                    t &= t - 1;
                    if (rowOfPiv[c] >= 0) {
                        fprintf(stderr, "z4: pivot col nonzero\n");
                        exit(7);
                    }
                    nz = true;
                }
            }
            if (nz)
                F.push_back(move(y));
        }
        vector<vector<uint64_t>>().swap(Yl);
        for (auto p : BL)
            free(p);
        for (auto p : BH)
            free(p);
        long long r1 = 0;               // F2 rank of the final leftovers = number of Z/2 summands
        for (int c = 0; c < m && r1 < (long long)F.size(); c++) {
            int w = c >> 6;
            uint64_t bit = 1ULL << (c & 63);
            long long piv = -1;
            for (long long i = r1; i < (long long)F.size(); i++)
                if (F[i][w] & bit) {
                    piv = i;
                    break;
                }
            if (piv < 0)
                continue;
            swap(F[piv], F[r1]);
            for (long long i = 0; i < (long long)F.size(); i++)
                if (i != r1 && (F[i][w] & bit))
                    for (int t = 0; t < W; t++)
                        F[i][t] ^= F[r1][t];
            r1++;
        }
        if (r1) {
            stages.push_back(r1);
            pivots += r1;
        }
        fprintf(stderr, "  [z4] stage 1 (Z/2 summands): %lld\n", r1);
    }

    // ---------- Z/9 bit-sliced dense phase -------------------------------------------------------------
    // An element of Z/9 is stored in binary as four bits (planes 0..3).  add9w adds two such 4-bit
    // numbers with a ripple-carry adder (result 0..16, five bits s0..s4), detects "s >= 9" as
    // s4 | (s3 & (s0|s1|s2)), and under that mask subtracts 9, implemented as adding 23 = 10111b mod 32.
    // Multiplication by a constant g is done by precomputing the multiples 0*y, 1*y, ..., 8*y.
    // F_3 vectors (for leftovers) use a different two-plane encoding: plane A marks entries equal to 1,
    // plane B marks entries equal to 2.  f3axpy computes A += f*B over F_3 in that encoding (f = 1, 2;
    // multiplying by 2 swaps the planes).
    static void f3axpy(vector<uint64_t> &A1, vector<uint64_t> &A2, const vector<uint64_t> &B1v,
                       const vector<uint64_t> &B2v, int f, int W) {
        for (int w = 0; w < W; w++) {
            uint64_t a1 = A1[w], a2 = A2[w], b1 = B1v[w], b2 = B2v[w];
            if (f == 2)
                swap(b1, b2);
            A1[w] = (a1 & ~(b1 | b2)) | (b1 & ~(a1 | a2)) | (a2 & b2);   // results equal to 1
            A2[w] = (a2 & ~(b1 | b2)) | (b2 & ~(a1 | a2)) | (a1 & b1);   // results equal to 2
        }
    }
    struct R9 {
        uint64_t *p[4];
    };
    static inline void add9w(uint64_t &x0, uint64_t &x1, uint64_t &x2, uint64_t &x3, uint64_t y0, uint64_t y1,
                             uint64_t y2, uint64_t y3) {
        uint64_t s0 = x0 ^ y0, c = x0 & y0;
        uint64_t t = x1 ^ y1;
        uint64_t s1 = t ^ c;
        c = (x1 & y1) | (c & t);
        t = x2 ^ y2;
        uint64_t s2 = t ^ c;
        c = (x2 & y2) | (c & t);
        t = x3 ^ y3;
        uint64_t s3 = t ^ c;
        uint64_t s4 = (x3 & y3) | (c & t);
        uint64_t m = s4 | (s3 & (s0 | s1 | s2));     // lanes where the sum is >= 9
        uint64_t r0 = s0 ^ m;                        // add 23 = 1,1,1,0,1 (bits 0..4) under the mask
        c = s0 & m;
        t = s1 ^ m;
        uint64_t r1 = t ^ c;
        c = (s1 & m) | (c & t);
        t = s2 ^ m;
        uint64_t r2 = t ^ c;
        c = (s2 & m) | (c & t);
        uint64_t r3 = s3 ^ c;
        x0 = r0;
        x1 = r1;
        x2 = r2;
        x3 = r3;
    }
    // Add two bit-sliced rows over Z/9 (four bit planes per entry).
    static inline void add9row(uint64_t *const *X, uint64_t *const *Y, int W) {
        uint64_t *a = X[0], *b = X[1], *c = X[2], *d = X[3];
        const uint64_t *e = Y[0], *f = Y[1], *g = Y[2], *h = Y[3];
        for (int w = 0; w < W; w++)
            add9w(a[w], b[w], c[w], d[w], e[w], f[w], g[w], h[w]);
    }
    // Read entry c of a bit-sliced Z/9 row.
    static inline int get9(uint64_t *const *X, int c) {
        int w = c >> 6, b = c & 63;
        return (int)((X[0][w] >> b) & 1) | (int)(((X[1][w] >> b) & 1) << 1) | (int)(((X[2][w] >> b) & 1) << 2) |
               (int)(((X[3][w] >> b) & 1) << 3);
    }
    // Write entry c of a bit-sliced Z/9 row.
    static inline void set9(uint64_t *const *X, int c, int v) {
        int w = c >> 6;
        uint64_t bit = 1ULL << (c & 63);
        for (int i = 0; i < 4; i++) {
            if ((v >> i) & 1)
                X[i][w] |= bit;
            else
                X[i][w] &= ~bit;
        }
    }
    // M[4g + i] = plane i of g*Y, for g = 0..8 (built by repeated addition of Y)
    static void multiples9(uint64_t *const *Y, vector<vector<uint64_t>> &M, int W) {
        M.assign(36, vector<uint64_t>(W, 0));
        auto P = [&](int g, int i) {
            return M[4 * g + i].data();
        };
        for (int i = 0; i < 4; i++)
            memcpy(P(1, i), Y[i], W * 8);
        for (int g = 2; g <= 8; g++) {
            for (int i = 0; i < 4; i++)
                memcpy(P(g, i), P(g - 1, i), W * 8);
            uint64_t *X[4] = {P(g, 0), P(g, 1), P(g, 2), P(g, 3)};
            uint64_t *Yy[4] = {P(1, 0), P(1, 1), P(1, 2), P(1, 3)};
            add9row(X, Yy, W);
        }
    }
    // Dense finishing phase over Z/9, as dense_z4 for p = 3.
    void dense_z9(long long sparseStagePivots) {
        vector<int> cmap(nc, -1);
        int m = 0;
        for (int c = 0; c < nc; c++)
            if (!coldead[c] && colcnt[c] > 0)
                cmap[c] = m++;
        vector<vector<int>>().swap(colrows);
        int W = (m + 63) / 64;
        fprintf(stderr, "  [z9] m=%d, basis memory <= %.2f GB\n", m, 4.0 * m * W * 8 / 1e9);
        vector<R9> B;
        vector<int> rowOfPiv(m, -1);
        vector<uint64_t> xs(4 * (size_t)W);
        uint64_t *X[4] = {xs.data(), xs.data() + W, xs.data() + 2 * W, xs.data() + 3 * W};
        vector<uint64_t> accs(36 * (size_t)W);   // accumulators S_g (g = 1..8) used while reducing a row
        vector<vector<uint64_t>> Mx;
        vector<vector<uint64_t>> Ya, Yb;         // online F3 RREF of halved leftovers y = x/3 mod 3
        vector<int> ypiv;
        auto insertY = [&](vector<uint64_t> &a, vector<uint64_t> &b) {
            for (size_t i = 0; i < Ya.size(); i++) {
                int c = ypiv[i];
                int w = c >> 6;
                uint64_t bit = 1ULL << (c & 63);
                int f = (a[w] & bit) ? 1 : ((b[w] & bit) ? 2 : 0);
                if (f)
                    f3axpy(a, b, Ya[i], Yb[i], 3 - f, W);
            }
            int c = -1;
            for (int w = 0; w < W && c < 0; w++) {
                uint64_t q = a[w] | b[w];
                if (q)
                    c = w * 64 + __builtin_ctzll(q);
            }
            if (c < 0)
                return;
            if ((b[c >> 6] >> (c & 63)) & 1)
                swap(a, b);                       // normalise the pivot to 1 (multiply by 2)
            for (size_t i = 0; i < Ya.size(); i++) {
                int w = c >> 6;
                uint64_t bit = 1ULL << (c & 63);
                int f = (Ya[i][w] & bit) ? 1 : ((Yb[i][w] & bit) ? 2 : 0);
                if (f)
                    f3axpy(Ya[i], Yb[i], a, b, 3 - f, W);
            }
            Ya.push_back(a);
            Yb.push_back(b);
            ypiv.push_back(c);
        };
        long long done = 0, stagePiv = sparseStagePivots;
        // Reduce the row given by its sparse entries against the basis.  The subtracted combination is
        // sum_j x[c_j] B[c_j]; basis rows are grouped by their coefficient g = -x[c_j] (mod 9), each group
        // is summed (one addition per basis row), and only then multiplied by g (g additions per group).
        auto reduce = [&](const vector<pair<int, int>> &ents) {
            fill(xs.begin(), xs.end(), 0);
            for (auto &[c, v] : ents)
                set9(X, c, v);
            bool used[9] = {false};
            for (auto &[c, v] : ents) {
                int br = rowOfPiv[c];
                if (br < 0)
                    continue;
                int f = get9(X, c);
                if (!f)
                    continue;
                int g = (9 - f) % 9;
                uint64_t *A[4] = {&accs[(4 * g) * (size_t)W], &accs[(4 * g + 1) * (size_t)W],
                                  &accs[(4 * g + 2) * (size_t)W], &accs[(4 * g + 3) * (size_t)W]};
                if (!used[g]) {
                    for (int i = 0; i < 4; i++)
                        memcpy(A[i], B[br].p[i], W * 8);
                    used[g] = true;
                } else
                    add9row(A, B[br].p, W);
            }
            for (int g = 1; g <= 8; g++)
                if (used[g]) {
                    uint64_t *A[4] = {&accs[(4 * g) * (size_t)W], &accs[(4 * g + 1) * (size_t)W],
                                      &accs[(4 * g + 2) * (size_t)W], &accs[(4 * g + 3) * (size_t)W]};
                    for (int t = 0; t < g; t++)
                        add9row(X, A, W);
                }
        };
        auto unitcol = [&]() {                   // first column holding a unit (value not divisible by 3)
            for (int w = 0; w < W; w++) {
                uint64_t nz = X[0][w] | X[1][w] | X[2][w] | X[3][w];
                while (nz) {
                    int c = w * 64 + __builtin_ctzll(nz);
                    nz &= nz - 1;
                    if (get9(X, c) % 3)
                        return c;
                }
            }
            return -1;
        };
        static const int inv9[9] = {0, 1, 5, 0, 7, 2, 0, 4, 8};   // inverses of the units mod 9
        for (int r = 0; r < nr; r++) {
            if (rowdead[r] || rows[r].empty())
                continue;
            vector<pair<int, int>> ents;
            ents.reserve(rows[r].size());
            for (auto e : rows[r]) {
                int c = cmap[col(e)];
                if (c < 0) {
                    fprintf(stderr, "z9: dead col\n");
                    exit(6);
                }
                ents.push_back({c, (int)(val(e) % 9)});
            }
            vector<uint32_t>().swap(rows[r]);
            reduce(ents);
            int j = unitcol();
            if (j >= 0) {
                int u = inv9[get9(X, j)];
                if (u != 1) {                    // normalise: x <- u*x
                    multiples9(X, Mx, W);
                    for (int i = 0; i < 4; i++)
                        memcpy(X[i], Mx[4 * u + i].data(), W * 8);
                }
                multiples9(X, Mx, W);            // multiples of the new pivot row
                long long nb = B.size();
#pragma omp parallel for schedule(static, 64)
                for (long long b = 0; b < nb; b++) {   // b <- b - b[j] * x  (keeps B reduced)
                    int f = get9(B[b].p, j);
                    if (f) {
                        int g = (9 - f) % 9;
                        uint64_t *Y[4] = {Mx[4 * g].data(), Mx[4 * g + 1].data(), Mx[4 * g + 2].data(),
                                          Mx[4 * g + 3].data()};
                        add9row(B[b].p, Y, W);
                    }
                }
                R9 nr9;
                for (int i = 0; i < 4; i++) {
                    nr9.p[i] = (uint64_t *)malloc(W * 8);
                    if (!nr9.p[i]) {
                        fprintf(stderr, "z9: out of memory\n");
                        exit(9);
                    }
                    memcpy(nr9.p[i], X[i], W * 8);
                }
                rowOfPiv[j] = B.size();
                B.push_back(nr9);
                stagePiv++;
                pivots++;
                rowdead[r] = 1;
            } else {                             // leftover: entries are 0, 3 or 6; store y = x/3 mod 3
                vector<uint64_t> ya(W, 0), yb(W, 0);
                bool any = false;
                for (int w = 0; w < W; w++) {
                    uint64_t nz = X[0][w] | X[1][w] | X[2][w] | X[3][w];
                    while (nz) {
                        int c = w * 64 + __builtin_ctzll(nz);
                        nz &= nz - 1;
                        int v = get9(X, c);
                        if (v % 3) {
                            fprintf(stderr, "z9: unit in leftover\n");
                            exit(7);
                        }
                        if (v == 3)
                            ya[w] |= 1ULL << (c & 63);
                        else
                            yb[w] |= 1ULL << (c & 63);
                        any = true;
                    }
                }
                if (any)
                    insertY(ya, yb);
            }
            if (++done % 200000 == 0)
                fprintf(stderr, "  [z9] %lld rows, basis %zu, leftover span %zu\n", done, B.size(), Ya.size());
        }
        stages.push_back(stagePiv);
        fprintf(stderr, "  [z9] stage 0 done: basis %zu, leftover span %zu\n", B.size(), Ya.size());
        // final reduction of each leftover against the stage-0 basis taken mod 3, then F3 rank
        map<int, pair<vector<uint64_t>, vector<uint64_t>>> bcache;   // basis rows mod 3, computed lazily
        auto bmod3 = [&](int br) -> pair<vector<uint64_t>, vector<uint64_t>> & {
            auto it = bcache.find(br);
            if (it != bcache.end())
                return it->second;
            vector<uint64_t> a(W, 0), b(W, 0);
            for (int w = 0; w < W; w++) {
                uint64_t nz = B[br].p[0][w] | B[br].p[1][w] | B[br].p[2][w] | B[br].p[3][w];
                while (nz) {
                    int c = w * 64 + __builtin_ctzll(nz);
                    nz &= nz - 1;
                    int v = get9(B[br].p, c) % 3;
                    if (v == 1)
                        a[w] |= 1ULL << (c & 63);
                    else if (v == 2)
                        b[w] |= 1ULL << (c & 63);
                }
            }
            return bcache[br] = {move(a), move(b)};
        };
        vector<vector<uint64_t>> F0, F1;
        for (size_t i = 0; i < Ya.size(); i++) {
            vector<uint64_t> a = Ya[i], b = Yb[i];
            vector<pair<int, int>> hits;         // (pivot column, value of y there)
            for (int w = 0; w < W; w++) {
                uint64_t q = a[w] | b[w];
                while (q) {
                    int c = w * 64 + __builtin_ctzll(q);
                    q &= q - 1;
                    if (rowOfPiv[c] >= 0)
                        hits.push_back({c, ((a[w] >> (c & 63)) & 1) ? 1 : 2});
                }
            }
            for (auto &[c, f] : hits) {
                auto &bb = bmod3(rowOfPiv[c]);
                f3axpy(a, b, bb.first, bb.second, 3 - f, W);
            }
            for (int w = 0; w < W; w++) {
                uint64_t q = a[w] | b[w];
                while (q) {
                    int c = w * 64 + __builtin_ctzll(q);
                    q &= q - 1;
                    if (rowOfPiv[c] >= 0) {
                        fprintf(stderr, "z9: pivot col nonzero after final reduction\n");
                        exit(7);
                    }
                }
            }
            bool nz = false;
            for (int w = 0; w < W && !nz; w++)
                if (a[w] | b[w])
                    nz = true;
            if (nz) {
                F0.push_back(move(a));
                F1.push_back(move(b));
            }
        }
        bcache.clear();
        for (auto &r9 : B)
            for (int i = 0; i < 4; i++)
                free(r9.p[i]);
        long long nrw = F0.size(), r1 = 0;
        auto getf = [&](long long i, int c) {
            int w = c >> 6;
            uint64_t bit = 1ULL << (c & 63);
            return (F0[i][w] & bit) ? 1 : ((F1[i][w] & bit) ? 2 : 0);
        };
        for (int c = 0; c < m && r1 < nrw; c++) {
            long long piv = -1;
            for (long long i = r1; i < nrw; i++)
                if (getf(i, c)) {
                    piv = i;
                    break;
                }
            if (piv < 0)
                continue;
            swap(F0[piv], F0[r1]);
            swap(F1[piv], F1[r1]);
            int pv = getf(r1, c);
            for (long long i = 0; i < nrw; i++)
                if (i != r1) {
                    int v = getf(i, c);                  // kill v using pivot pv: add (-v/pv) * pivot row
                    if (v)
                        f3axpy(F0[i], F1[i], F0[r1], F1[r1], (v == pv) ? 2 : 1, W);
                }
            r1++;
        }
        if (r1) {
            stages.push_back(r1);
            pivots += r1;
        }
        fprintf(stderr, "  [z9] stage 1 (Z/3 summands): %lld\n", r1);
    }

    // ---------- generic dense phase (byte per entry), all stages ----------------------------------------
    void dense_finish(long long sparseStagePivots) {
        if (p0 == 2 && mod == 4) {
            dense_z4(sparseStagePivots);
            return;
        }
        if (p0 == 3 && mod == 9) {
            dense_z9(sparseStagePivots);
            return;
        }
        vector<int> cmap(nc, -1);
        int m = 0;
        for (int c = 0; c < nc; c++)
            if (!coldead[c] && colcnt[c] > 0)
                cmap[c] = m++;
        vector<vector<int>>().swap(colrows);
        vector<pair<int, vector<pair<int, uint8_t>>>> in;   // (original row index, sparse dense-indexed row)
        for (int r = 0; r < nr; r++)
            if (!rowdead[r] && !rows[r].empty()) {
                vector<pair<int, uint8_t>> v;
                v.reserve(rows[r].size());
                for (auto x : rows[r]) {
                    int cc = cmap[col(x)];
                    if (cc < 0) {
                        fprintf(stderr, "dense: entry in dead column\n");
                        exit(6);
                    }
                    v.push_back({cc, (uint8_t)val(x)});
                }
                vector<uint32_t>().swap(rows[r]);
                in.push_back({r, move(v)});
            }
        fprintf(stderr, "  [dense] m=%d input rows=%zu, est. basis memory <= %.2f GB\n", m, in.size(),
                (double)m * m / 1e9);
        long long stagePiv = sparseStagePivots;
        bool first = true;
        while (true) {                                      // one iteration per stage
            vector<uint8_t *> B;
            vector<int> rowOfPiv(m, -1);
            vector<uint8_t> x(m);
            vector<vector<pair<int, uint8_t>>> left;
            auto reduce = [&](const vector<pair<int, uint8_t>> &sp) {
                fill(x.begin(), x.end(), 0);
                for (auto &[c, v] : sp)
                    x[c] = v;
                for (auto &[c, v] : sp) {
                    int br = rowOfPiv[c];
                    if (br >= 0 && x[c])
                        axpyD(x.data(), B[br], x[c], m);
                }
            };
            long long done = 0;
            for (auto &[orig, sp] : in) {
                reduce(sp);
                int j = -1;                                 // first unit entry (skip zero words fast)
                for (int t0 = 0; t0 < m && j < 0; t0 += 8) {
                    uint64_t w = 0;
                    if (t0 + 8 <= m)
                        memcpy(&w, x.data() + t0, 8);
                    else
                        w = 1;
                    if (!w)
                        continue;
                    for (int t = t0; t < min(m, t0 + 8); t++)
                        if (x[t] % p0) {
                            j = t;
                            break;
                        }
                }
                if (j >= 0) {
                    scalD(x.data(), inv(x[j]), m);
                    {
                        long long nb2 = B.size();
                        uint8_t *xd = x.data();
#pragma omp parallel for schedule(static, 64)
                        for (long long b2 = 0; b2 < nb2; b2++)
                            if (B[b2][j])
                                axpyD(B[b2], xd, B[b2][j], m);
                    }
                    uint8_t *nb = (uint8_t *)malloc(m);
                    memcpy(nb, x.data(), m);
                    rowOfPiv[j] = B.size();
                    B.push_back(nb);
                    stagePiv++;
                    pivots++;
                    if (first && orig >= 0)
                        rowdead[orig] = 1;                  // only stage-0 pivots count for column dropping
                } else {
                    vector<pair<int, uint8_t>> rest;
                    for (int t = 0; t < m; t++)
                        if (x[t])
                            rest.push_back({t, x[t]});
                    if (!rest.empty())
                        left.push_back(move(rest));
                }
                vector<pair<int, uint8_t>>().swap(sp);
                if (++done % 200000 == 0)
                    fprintf(stderr, "  [dense] mod %u: %lld rows, basis %zu, leftover %zu\n", mod, done, B.size(),
                            left.size());
            }
            stages.push_back(stagePiv);
            stagePiv = 0;
            first = false;
            fprintf(stderr, "  [dense] stage done: mod %u basis %zu leftover rows %zu\n", mod, B.size(), left.size());
            // Re-reduce leftovers against the FINAL basis (they were reduced against an earlier, smaller
            // one).  Survivors are divisible by p and supported on the non-pivot columns.
            vector<pair<int, vector<pair<int, uint8_t>>>> nxt;
            vector<int> cmap2(m, -1);
            int m2 = 0;
            for (int c = 0; c < m; c++)
                if (rowOfPiv[c] < 0)
                    cmap2[c] = m2++;
            for (auto &sp : left) {
                reduce(sp);
                vector<pair<int, uint8_t>> rest;
                for (int t = 0; t < m; t++)
                    if (x[t]) {
                        if (x[t] % p0) {
                            fprintf(stderr, "dense: unit survived\n");
                            exit(7);
                        }
                        if (rowOfPiv[t] >= 0) {
                            fprintf(stderr, "dense: pivot col nonzero\n");
                            exit(7);
                        }
                        rest.push_back({cmap2[t], (uint8_t)(x[t])});
                    }
                if (!rest.empty())
                    nxt.push_back({-1, move(rest)});
            }
            for (auto p : B)
                free(p);
            if (nxt.empty())
                break;
            if (mod == p0) {                                // no precision left: count what survives
                for (auto &r : nxt)
                    leftover_nonzero += r.second.size();
                break;
            }
            for (auto &r : nxt) {                           // divide by p
                vector<pair<int, uint8_t>> M2;
                for (auto &[c, v] : r.second)
                    if (v / p0)
                        M2.push_back({c, (uint8_t)(v / p0)});
                r.second.swap(M2);
            }
            mod /= p0;
            m = m2;
            in.swap(nxt);
            in.erase(remove_if(in.begin(), in.end(), [](auto &r) { return r.second.empty(); }), in.end());
            if (in.empty())
                break;
        }
    }

    // Full local elimination: sparse stages, dividing by p between stages; possibly a dense finish.
    void run() {
        long long before = 0;
        if (denseMax < 0) {                                 // -D -1: dense from the start (testing)
            dense_finish(0);
            return;
        }
        while (true) {
            while (sweep()) {
                if (goDense)
                    break;
            }
            if (goDense) {
                dense_finish(pivots - before);
                return;
            }
            stages.push_back(pivots - before);
            before = pivots;
            bool live = false;
            for (int r = 0; r < nr && !live; r++)
                if (!rowdead[r] && !rows[r].empty())
                    live = true;
            if (!live)
                break;
            if (mod == p0) {
                for (int r = 0; r < nr; r++)
                    if (!rowdead[r])
                        leftover_nonzero += rows[r].size();
                break;
            }
            for (int r = 0; r < nr; r++)                    // every live entry is divisible by p: divide
                if (!rowdead[r]) {
                    vector<uint32_t> M;
                    for (auto x : rows[r]) {
                        if (val(x) % p0) {
                            fprintf(stderr, "local: nonunit not divisible\n");
                            exit(4);
                        }
                        uint32_t v = val(x) / p0;
                        if (v)
                            M.push_back(col(x) << 5 | v);
                        else
                            colcnt[col(x)]--;
                    }
                    rows[r].swap(M);
                }
            mod /= p0;
        }
    }
};

#include "gr42.h"   // LGR: the GR(4,2) engine (sparse + bit-sliced dense)

// Placeholder "residual file" for engines that have no residual (all field / local engines): an empty
// 0 x 0 matrix, so that downstream tools can treat every mode uniformly.
void empty_residual(const string &fn, int &rr, int &rc, long long &rnnz) {
    FILE *f = fopen(fn.c_str(), "w");
    fprintf(f, "0 0\n");
    fclose(f);
    rr = rc = 0;
    rnnz = 0;
}

// Split "a,b,c" into integers.
static vector<int> parse_ints(const string &s) {
    vector<int> out;
    stringstream ss(s);
    string t;
    while (getline(ss, t, ','))
        out.push_back(stoi(t));
    return out;
}

// Closure of a set of permutations of {0..n-1} under composition (breadth-first search).  The identity
// is always included.  Returned sorted.
static vector<vector<int>> closeup(const vector<vector<int>> &gens, int n) {
    vector<int> id(n);
    iota(id.begin(), id.end(), 0);
    set<vector<int>> G = {id};
    vector<vector<int>> fr = {id};
    while (!fr.empty()) {
        vector<vector<int>> nf;
        for (auto &x : fr)
            for (auto &g : gens) {
                vector<int> y(n);
                for (int t = 0; t < n; t++)
                    y[t] = g[x[t]];
                if (G.insert(y).second)
                    nf.push_back(y);
            }
        fr.swap(nf);
    }
    return vector<vector<int>>(G.begin(), G.end());
}

// Extend values given on generators to a homomorphism on the generated group, by breadth-first search,
// combining with 'op'; aborts (exit 12) if two paths to the same element give different values.
template <class Op>
static map<vector<int>, int> extend_character(const vector<vector<int>> &gens, const vector<int> &vals, int dflt,
                                              int identity, Op op, const char *what) {
    int n = gens.empty() ? V + 2 * E : gens[0].size();
    map<vector<int>, int> val;
    vector<int> id(n);
    iota(id.begin(), id.end(), 0);
    val[id] = identity;
    vector<vector<int>> fr = {id};
    while (!fr.empty()) {
        vector<vector<int>> nf;
        for (auto &x : fr)
            for (size_t gi = 0; gi < gens.size(); gi++) {
                auto &g = gens[gi];
                vector<int> y(n);
                for (int t = 0; t < n; t++)
                    y[t] = g[x[t]];
                int vy = op(val[x], gi < vals.size() ? vals[gi] : dflt);
                auto it = val.find(y);
                if (it == val.end()) {
                    val[y] = vy;
                    nf.push_back(y);
                } else if (it->second != vy) {
                    fprintf(stderr, "%s do not define a character\n", what);
                    exit(12);
                }
            }
        fr.swap(nf);
    }
    return val;
}

int main(int argc, char **argv) {
    // ------------------------------------------------------------------ command line
    ll P = 0;                 // coefficient modulus: 0 = Z; with -q, P = p^K
    string prefix = "res";
    int maxdeg = -1;
    bool cellsonly = false;
    bool genericMatrix = false;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (a == "-p")
            P = atoll(argv[++i]);
        else if (a == "-q") {
            LQ = atoll(argv[++i]);
            LK = atoi(argv[++i]);
            P = 1;
            for (int t = 0; t < LK; t++)
                P *= LQ;
        } else if (a == "-o")
            prefix = argv[++i];
        else if (a == "-d")
            maxdeg = atoi(argv[++i]);
        else if (a == "-D")
            LPE::denseMax = atoll(argv[++i]);
        else if (a == "-M") {
            LPE::denseMin = atoll(argv[++i]);
            LPE::nnzBudget = atoll(argv[++i]);
        } else if (a == "-C")
            LPE::cntThresh = atoi(argv[++i]);
        else if (a == "-G")
            LGR::denseMax = atoll(argv[++i]);
        else if (a == "-F")
            fastEnum = true;
        else if (a == "-W") {
            grMode = true;
            agW = parse_ints(argv[++i]);
        } else if (a == "-Y")
            agSigns = parse_ints(argv[++i]);
        else if (a == "-X") {
            string t = argv[++i];
            twistLab = t.find('l') != string::npos;
            twistEdge = t.find('e') != string::npos;
        } else if (a == "-a") {
            grouped = true;
            stringstream gs(argv[++i]);
            string one;
            while (getline(gs, one, ';'))
                agens.push_back(one);
        } else if (a == "-g") {
            grouped = true;
            stringstream gs(argv[++i]);
            string one;
            while (getline(gs, one, ';'))
                gelts.push_back(parse_ints(one));
        } else if (a == "-k") {
            colored = true;
            string cs = argv[++i];
            for (char ch : cs)
                pcolor.push_back(ch - 'a');
        } else if (a == "-c")
            cellsonly = true;
        else if (a == "-z")
            checksq = true;
        else if (a == "-T")
            transposeTop = true;
        else if (a == "-m")
            genericMatrix = true;
        else {
            fprintf(stderr, "unknown option %s\n", a.c_str());
            return 1;
        }
    }
    if (grMode)
        fastEnum = true;      // the slow enumerator does not know about omega characters
    if (grMode && checksq) {
        fprintf(stderr, "-z is not supported together with -W\n");
        return 1;
    }

    // ------------------------------------------------------------------ generic matrix mode (-m)
    if (genericMatrix) {
        int nr, nc;
        if (scanf("%d %d", &nr, &nc) != 2)
            return 1;
        SparseElim M(nr, nc, P);
        int r, c;
        ll v;
        while (scanf("%d %d %lld", &r, &c, &v) == 3)
            M.rows[r].push_back({c, (int)v});
        M.finalize();
        M.run();
        int rr, rc;
        long long rnnz;
        M.dump_residual(prefix + "_m.txt", rr, rc, rnnz);
        printf("{\"pivots\":%lld,\"res_rows\":%d,\"res_cols\":%d,\"res_file\":\"%s\"}\n", M.pivots, rr, rc,
               (prefix + "_m.txt").c_str());
        return 0;
    }

    // ------------------------------------------------------------------ read the graph
    if (scanf("%d %d %d", &V, &E, &N) != 3)
        return 1;
    etail.resize(E);
    ehead.resize(E);
    hvert.resize(2 * E);
    vhalf.assign(V, {});
    for (int e = 0; e < E; e++) {
        if (scanf("%d %d", &etail[e], &ehead[e]) != 2)
            return 1;
        hvert[2 * e] = etail[e];
        hvert[2 * e + 1] = ehead[e];
        vhalf[etail[e]].push_back(2 * e);
        vhalf[ehead[e]].push_back(2 * e + 1);
    }
    if (3 * E >= 255 || N > 10) {
        fprintf(stderr, "too large\n");
        return 1;
    }

    // ------------------------------------------------------------------ group and base half-edges
    base.assign(V, -1);
    if (!grouped) {
        for (int v = 0; v < V; v++)
            if (!vhalf[v].empty())
                base[v] = vhalf[v][0];   // lowest-numbered half-edge
    } else {
        for (auto &g : gelts)            // extend label generators by the identity
            for (int t = g.size(); t < N; t++)
                g.push_back(t);
        auto Lall = closeup(gelts, N);

        // Automorphism generators.  Each is turned into one permutation 'comb' of V + 2E points:
        // points 0..V-1 are vertices; point V + h is half-edge h.  An edge mapped with reversal sends its
        // tail half-edge to the image's head half-edge and vice versa.
        vector<vector<int>> ag;
        for (auto &str : agens) {
            string vs = str, es;
            size_t bar = str.find('|');
            if (bar != string::npos) {
                vs = str.substr(0, bar);
                es = str.substr(bar + 1);
            }
            vector<int> vm = parse_ints(vs), em, fl;
            if ((int)vm.size() != V) {
                fprintf(stderr, "bad automorphism\n");
                exit(8);
            }
            if (!es.empty()) {           // explicit edge images; "~e" = reversed
                stringstream ss(es);
                string t;
                while (getline(ss, t, ',')) {
                    bool f = t[0] == '~';
                    em.push_back(stoi(f ? t.substr(1) : t));
                    fl.push_back(f);
                }
            } else {                     // match parallel copies in input order, allowing reversal
                map<pair<int, int>, vector<int>> byend;
                for (int e = 0; e < E; e++)
                    byend[{etail[e], ehead[e]}].push_back(e);
                map<pair<int, int>, int> seen;
                em.assign(E, -1);
                fl.assign(E, 0);
                for (int e = 0; e < E; e++) {
                    int ci = seen[{etail[e], ehead[e]}]++;   // e is the ci-th edge with its endpoints
                    pair<int, int> key = {vm[etail[e]], vm[ehead[e]]};
                    if (!byend.count(key)) {
                        key = {key.second, key.first};
                        fl[e] = 1;
                    }
                    auto &L = byend[key];
                    if (ci >= (int)L.size()) {
                        fprintf(stderr, "automorphism does not map edges to edges\n");
                        exit(8);
                    }
                    em[e] = L[ci];
                }
            }
            for (int e = 0; e < E; e++) {   // verify that this really is a graph automorphism
                int a = fl[e] ? ehead[em[e]] : etail[em[e]], b = fl[e] ? etail[em[e]] : ehead[em[e]];
                if (a != vm[etail[e]] || b != vm[ehead[e]]) {
                    fprintf(stderr, "automorphism does not preserve edge %d\n", e);
                    exit(8);
                }
            }
            vector<int> comb(V + 2 * E);
            for (int v = 0; v < V; v++)
                comb[v] = vm[v];
            for (int e = 0; e < E; e++) {
                comb[V + 2 * e] = V + 2 * em[e] + fl[e];
                comb[V + 2 * e + 1] = V + 2 * em[e] + 1 - fl[e];
            }
            ag.push_back(comb);
        }
        auto Aall = closeup(ag, V + 2 * E);
        // characters on the automorphism group: omega exponents (-W, additive mod 3) and signs (-Y)
        auto awx = extend_character(ag, agW, 0, 0, [](int a, int b) { return (a + b) % 3; }, "-W exponents");
        auto achi = extend_character(ag, agSigns, 1, 1, [](int a, int b) { return a * b; }, "-Y signs");

        // Invariant base half-edges.  For each vertex orbit, try the half-edges h of a representative v in
        // turn; h works if every group element A sends (v, h) to a consistent (A(v), A(h)), i.e. the
        // stabiliser of v fixes h.  Then base[A(v)] = A(h) for all A.  If no h works, the reduced complex
        // with these bases is not preserved by the group and we stop (exit 8).
        for (int v = 0; v < V; v++)
            if (base[v] < 0 && !vhalf[v].empty()) {
                bool ok = false;
                for (int h : vhalf[v]) {
                    map<int, int> asg;
                    bool good = true;
                    for (auto &A : Aall) {
                        int w = A[v];
                        int hh = A[V + h] - V;
                        auto it = asg.find(w);
                        if (it != asg.end() && it->second != hh) {
                            good = false;
                            break;
                        }
                        asg[w] = hh;
                    }
                    if (good) {
                        for (auto &[w, hh] : asg)
                            base[w] = hh;
                        ok = true;
                        break;
                    }
                }
                if (!ok) {
                    fprintf(stderr, "no invariant base half-edge at %d\n", v);
                    exit(8);
                }
            }

        // All elements of G = L x A with their character values.
        for (auto &Lp : Lall)
            for (auto &A : Aall) {
                GE g;
                g.lab = Lp;
                g.vm.assign(A.begin(), A.begin() + V);
                g.em.resize(E);
                g.fl.resize(E);
                for (int e = 0; e < E; e++) {
                    int x = A[V + 2 * e] - V;
                    g.em[e] = x / 2;
                    g.fl[e] = x & 1;
                }
                bool id = true;
                for (int t = 0; t < N; t++)
                    if (g.lab[t] != t)
                        id = false;
                for (int v = 0; v < V; v++)
                    if (g.vm[v] != v)
                        id = false;
                for (int e = 0; e < E; e++)
                    if (g.em[e] != e || g.fl[e])
                        id = false;
                for (int h = 0; h < 2 * E; h++)
                    g.hm[h] = 2 * g.em[h / 2] + ((h & 1) ^ g.fl[h / 2]);
                g.chi = (twistLab ? permsign(g.lab) : 1) * (twistEdge ? permsign(g.em) : 1) * achi[A];
                g.w = awx[A];
                if (!id)
                    gall.push_back(g);
            }
        // A fixed pseudo-random order makes the early-abort minimality tests fail fast on average.
        {
            mt19937 rg(12345);
            shuffle(gall.begin(), gall.end(), rg);
        }
        for (auto &g : gall) {
            bool labid = true, autid = true;
            for (int t = 0; t < N; t++)
                if (g.lab[t] != t)
                    labid = false;
            for (int v = 0; v < V; v++)
                if (g.vm[v] != v)
                    autid = false;
            for (int e = 0; e < E; e++)
                if (g.em[e] != e || g.fl[e])
                    autid = false;
            if (autid)
                gLab.push_back(&g);
            else if (labid)
                gAut.push_back(&g);
        }
        {
            GE idg;
            idg.lab.resize(N);
            iota(idg.lab.begin(), idg.lab.end(), 0);
            idg.vm.resize(V);
            iota(idg.vm.begin(), idg.vm.end(), 0);
            idg.em.resize(E);
            iota(idg.em.begin(), idg.em.end(), 0);
            idg.fl.assign(E, 0);
            idg.chi = 1;
            for (int h = 0; h < 2 * E; h++)
                idg.hm[h] = h;
            Aelts.push_back(idg);
            Lelts.push_back(idg);
            for (auto *g : gAut)
                Aelts.push_back(*g);
            for (auto *g : gLab)
                Lelts.push_back(*g);
            if (Aelts.size() * Lelts.size() != gall.size() + 1) {
                fprintf(stderr, "factorisation failed\n");
                exit(11);
            }
        }
        // Coinvariants compute the homology of the quotient only if |G| is invertible.
        long long order = gall.size() + 1;
        if ((LQ && order % LQ == 0) || (grMode && order % 2 == 0)) {
            fprintf(stderr, "group order %lld divisible by p: coinvariants invalid\n", order);
            exit(10);
        }
        fprintf(stderr, "group order %zu (labels %zu x automorphisms %zu); bases:", gall.size() + 1, Lall.size(),
                Aall.size());
        for (int v = 0; v < V; v++)
            fprintf(stderr, " %d", base[v]);
        fprintf(stderr, "\n");
    }

    // ------------------------------------------------------------------ top degree D
    int D = 0;   // number of vertices that can carry a reduced generator (valence >= 2)
    for (int v = 0; v < V; v++)
        if (vhalf[v].size() >= 2)
            D++;
    D = min(D, N);
    if (maxdeg >= 0)
        D = min(D, maxdeg);

    // The boundary of a (decoded) cell c, as a list of (canonical face key, coefficient sign) pairs.
    // Faces that vanish in the quotient are skipped.  The omega exponent of each face is reported via
    // 'ws' (GR mode).
    auto boundary = [&](const Cell &c, vector<pair<Key, int>> &faces, vector<int> &ws) {
        faces.clear();
        ws.clear();
        int idx = 0;
        for (int v = 0; v < V; v++)
            if (c.vocc[v] >= 0) {
                int p = c.vocc[v], h = c.vh[v], h0 = base[v];
                int sgn = (idx % 2 == 0) ? 1 : -1;   // (-1)^(position among occupied vertices)
                idx++;
                for (int t = 0; t < 2; t++) {        // t = 0: push into edge(h);  t = 1: into edge(h0)
                    int hh = t == 0 ? h : h0;
                    Cell f = c;
                    f.vocc[v] = -1;
                    f.vh[v] = -1;
                    int e = hh / 2;
                    if (hh % 2 == 0)
                        f.elist[e].insert(f.elist[e].begin(), p);   // tail end: first in the edge's list
                    else
                        f.elist[e].push_back(p);                    // head end: last in the list
                    int sf = canonicalize(f);
                    if (!sf)
                        continue;
                    faces.push_back({encode(f), (t == 0 ? sgn : -sgn) * sf});
                    ws.push_back(lastW);
                }
            }
    };

    // ------------------------------------------------------------------ -z : check d o d = 0
    if (checksq) {
        vector<vector<Key>> C(D + 1);
        for (int k = 0; k <= D; k++) {
            enumerate(k, C[k]);
            sort(C[k].begin(), C[k].end());
            C[k].erase(unique(C[k].begin(), C[k].end()), C[k].end());
        }
        auto index_of = [&](int k, Key fk) -> long long {
            auto it = lower_bound(C[k].begin(), C[k].end(), fk);
            if (it == C[k].end() || *it != fk) {
                fprintf(stderr, "face not found\n");
                exit(2);
            }
            return it - C[k].begin();
        };
        long long bad = 0, checked = 0;
        vector<pair<Key, int>> f1, f2;
        vector<int> w1, w2;
        for (int k = 2; k <= D; k++)
            for (size_t i = 0; i < C[k].size(); i++) {
                boundary(decode(C[k][i]), f1, w1);
                map<long long, long long> acc;
                for (auto &x : f1) {
                    index_of(k - 1, x.first);
                    boundary(decode(x.first), f2, w2);
                    for (auto &y : f2)
                        acc[index_of(k - 2, y.first)] += (long long)x.second * y.second;
                }
                for (auto &z : acc)
                    if (z.second != 0) {
                        bad++;
                        break;
                    }
                checked++;
            }
        printf("{\"d2_checked\":%lld,\"d2_nonzero\":%lld}\n", checked, bad);
        return 0;
    }

    // ------------------------------------------------------------------ main loop over degrees
    vector<long long> csize(D + 2, 0);
    vector<Key> prevcells;   // sorted (k-1)-cells
    vector<char> dropcol;    // dropcol[i] = 1 if (k-1)-cell i was a pivot row of d_{k-1}
    bool canDrop = true;     // may those columns be dropped from d_k?
    printf("{\"V\":%d,\"E\":%d,\"n\":%d,\"P\":%lld,\"dim\":%d", V, E, N, P, D);
    vector<string> results;
    for (int k = 0; k <= D; k++) {
        vector<Key> cells;
        enumerate(k, cells);
        sort(cells.begin(), cells.end());
        if (colored || grouped)
            cells.erase(unique(cells.begin(), cells.end()), cells.end());
        csize[k] = cells.size();
        if (k >= 1 && !cellsonly) {
            // Store one matrix entry in the format of the engine M.
            auto addent = [&](auto &MM, size_t i, int ci, ll coef) {
                using T = std::decay_t<decltype(MM)>;
                if constexpr (std::is_same_v<T, F2E>) {
                    if (coef & 1)
                        MM.rows[i].push_back(ci);
                } else if constexpr (std::is_same_v<T, LPE>) {
                    ll v = ((coef % P) + P) % P;
                    if (v)
                        MM.rows[i].push_back((uint32_t)ci << 5 | (uint32_t)v);
                } else if constexpr (std::is_same_v<T, LGR>) {
                    if (coef)
                        MM.rows[i].push_back((uint32_t)ci << 5 | (uint32_t)coef);
                } else if constexpr (std::is_same_v<T, SPE>) {
                    ll v = ((coef % P) + P) % P;
                    if (v)
                        MM.rows[i].push_back((uint32_t)ci << 4 | (uint32_t)v);
                } else
                    MM.rows[i].push_back({ci, (int)coef});
            };
            auto dumpres = [&](auto &MM, const string &fn, int &rr, int &rc, long long &rnnz) {
                if constexpr (std::is_same_v<std::decay_t<decltype(MM)>, SparseElim>)
                    MM.dump_residual(fn, rr, rc, rnnz);
                else
                    empty_residual(fn, rr, rc, rnnz);
            };
            // Build the matrix of d_k (rows = k-cells, in parallel), eliminate, record the result.
            auto blockfn = [&](auto &M) {
                using T = std::decay_t<decltype(M)>;
#pragma omp parallel for schedule(dynamic, 4096)
                for (size_t i = 0; i < cells.size(); i++) {
                    vector<pair<Key, int>> faces;
                    vector<int> ws;
                    boundary(decode(cells[i]), faces, ws);
                    for (size_t t = 0; t < faces.size(); t++) {
                        Key fk = faces[t].first;
                        auto it = lower_bound(prevcells.begin(), prevcells.end(), fk);
                        if (it == prevcells.end() || *it != fk) {
                            fprintf(stderr, "face not found\n");
                            exit(2);
                        }
                        ll coef = faces[t].second;
                        if (grMode)
                            coef = grcode((int)coef, ws[t]);     // (+-1) * w^t as a GR(4,2) element
                        else if (P)
                            coef = ((coef % P) + P) % P;
                        int ci = (int)(it - prevcells.begin());
                        if (!(canDrop && !dropcol.empty() && dropcol[ci]))
                            addent(M, i, ci, coef);
                    }
                }
                vector<Key>().swap(prevcells);   // free before elimination
                if constexpr (std::is_same_v<T, F2E>) {
                    if (transposeTop && k == D) {   // rank only, via the transpose (rows = (k-1)-cells)
                        F2E Tm(M.nc, M.nr);
                        for (int i = 0; i < M.nr; i++) {
                            for (int c : M.rows[i])
                                Tm.rows[c].push_back(i);
                            vector<int>().swap(M.rows[i]);
                        }
                        Tm.finalize();
                        Tm.run();
                        M.pivots = Tm.pivots;
                        char b2[256];
                        snprintf(b2, sizeof b2,
                                 "{\"k\":%d,\"rows\":%zu,\"cols\":%d,\"pivots\":%lld,\"res_rows\":0,\"res_cols\":0,"
                                 "\"res_nnz\":0,\"maxabs\":1,\"res_file\":\"none\"}",
                                 k, cells.size(), M.nc, M.pivots);
                        results.push_back(b2);
                        canDrop = false;
                        return;
                    }
                }
                M.finalize();
                fprintf(stderr, "deg %d: matrix built, eliminating\n", k);
                M.run();
                int rr, rc;
                long long rnnz;
                string fn = prefix + "_" + to_string(k) + ".txt";
                dumpres(M, fn, rr, rc, rnnz);
                char buf[512];
                snprintf(buf, sizeof buf,
                         "{\"k\":%d,\"rows\":%zu,\"cols\":%zu,\"pivots\":%lld,\"res_rows\":%d,\"res_cols\":%d,"
                         "\"res_nnz\":%lld,\"maxabs\":%lld,\"res_file\":\"%s\"}",
                         k, cells.size(), prevcells.size(), M.pivots, rr, rc, rnnz, M.maxabs, fn.c_str());
                string rec = buf;
                if constexpr (std::is_same_v<T, LPE> || std::is_same_v<T, LGR>) {
                    rec.pop_back();
                    rec += ",\"stages\":[";
                    for (size_t t = 0; t < M.stages.size(); t++)
                        rec += (t ? "," : "") + to_string(M.stages[t]);
                    rec += "],\"leftover_nnz\":" + to_string(M.leftover_nonzero) + "}";
                    canDrop = (M.stages.size() <= 1 && M.leftover_nonzero == 0);   // only unit pivots used
                } else
                    canDrop = (rr == 0);
                results.push_back(rec);
                dropcol.assign(cells.size(), 0);
                for (size_t i = 0; i < cells.size(); i++)
                    if (M.rowdead[i])
                        dropcol[i] = 1;
            };
            if (grMode) {
                LGR M(cells.size(), prevcells.size());
                blockfn(M);
            } else if (LQ) {
                LPE M(cells.size(), prevcells.size(), (uint32_t)LQ, (uint32_t)P);
                blockfn(M);
            } else if (P == 2) {
                F2E M(cells.size(), prevcells.size());
                blockfn(M);
            } else if (P > 2 && P <= 13) {
                SPE M(cells.size(), prevcells.size(), (uint32_t)P);
                blockfn(M);
            } else {
                SparseElim M(cells.size(), prevcells.size(), P);
                blockfn(M);
            }
        }
        prevcells.swap(cells);
        fprintf(stderr, "deg %d: %lld cells\n", k, csize[k]);
    }
    printf(",\"cells\":[");
    for (int k = 0; k <= D; k++)
        printf("%s%lld", k ? "," : "", csize[k]);
    printf("],\"maps\":[");
    for (size_t i = 0; i < results.size(); i++)
        printf("%s%s", i ? "," : "", results[i].c_str());
    printf("]}\n");
    return 0;
}
