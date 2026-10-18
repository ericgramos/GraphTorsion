// =====================================================================================================
//  prod2.cpp  --  are top classes of a two-vertex graph products of star classes?
// =====================================================================================================
//  SETTING.  G is a graph with exactly two essential (valence >= 2) vertices u and w; every other vertex is
//  a leaf.  Then the reduced Swiatkowski complex of Conf_n(G) lives in degrees 0, 1, 2 and the top homology
//  H_2(Conf_n(G)) = Z_2 = ker d_2 is free.
//
//  QUESTION.  Is H_2 spanned by STAR x STAR TORI, i.e. by the classes
//        c(z_u (x) z_w),   z_u in H_1(Conf_A(Star_u)),  z_w in H_1(Conf_B(Star_w)),  A |_| B = [n] ?
//  Here Star_v is the open star of v (v with its d_v half-edges as arms).  The reduced complex of
//  Conf_k(Star_d) has cells in degrees 0 and 1 only (the centre is occupied or not), so it is a graph and
//  H_1 is its cycle space.  A basis of that cycle space is given by the fundamental cycles of a BFS
//  spanning tree (star_basis below); all coefficients are +-1.  The product of a cycle at u (particles A)
//  and a cycle at w (particles B) is written as a 2-chain of G: the arms of the two stars are glued to the
//  edges of G, an edge between u and w receiving the particles of u's arm (listed from u outward) followed
//  by those of w's arm (in reverse order).
//
//  OUTPUT.  {"n","prime","du","dw","C2" = |C_2|, "b2" = dim ker d_2 mod prime, "gens" = number of torus
//  generators, "rankP" = rank of their span mod prime}.  rankP == b2 means the tori span H_2 (mod prime);
//  doing this for a large prime and for p gives the statement over Z_(p) and Q.
//  A random sample of 200 generators is checked to be cycles (exit 3 otherwise).
//
//  USAGE.   prod2 n prime [k] < graph        (optional k: only tori with |A| = k)
//  USED FOR.  Showing that the torsion of T(a,b,c) lies in the image of
//  star x star tori of Gamma minus a vertex (prod_T222.log, prod_T223.log).
//
//  SHARED HELPER CODE.  Everything above main() in this file (struct Graph, struct Cell, encode/decode of
//  128-bit cell keys, enumerate(), boundary(), lookup(), the F_p eliminator FpElim and, where present,
//  star_basis()/z1_basis()) is the same helper code as in prod2.cpp, where it is documented line by line.
//  The conventions are exactly those of engine/oconf.cpp: half-edge 2e is the tail end of edge e, 2e+1 the
//  head end; the base half-edge of a vertex is its lowest-numbered half-edge; particle p occupies bits
//  [12p, 12p+12) of a key (8 bits location, 4 bits rank on the edge).
//
//  GRAPH INPUT FORMAT for all tools in this directory (note: no particle number on the first line):
//      V E
//      u_0 v_0
//      ...            (edge e oriented u_e -> v_e)
// =====================================================================================================

#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef unsigned long long u64;
typedef unsigned __int128 Key;

// A graph with its half-edge structure.  build() numbers half-edges as in the engine (2e = tail end of
// edge e, 2e+1 = head end), records for each half-edge its vertex (hvert) and for each vertex its sorted
// half-edges (vhalf), and chooses the lowest half-edge at each vertex as base half-edge.
struct Graph {
    int V = 0, E = 0;
    vector<int> tail, head, hvert;
    vector<vector<int>> vhalf;
    vector<int> base;
    void build(int V_, const vector<pair<int, int>> &ed) {
        V = V_;
        E = ed.size();
        tail.resize(E);
        head.resize(E);
        hvert.resize(2 * E);
        vhalf.assign(V, {});
        for (int e = 0; e < E; e++) {
            tail[e] = ed[e].first;
            head[e] = ed[e].second;
            hvert[2 * e] = tail[e];
            hvert[2 * e + 1] = head[e];
            vhalf[tail[e]].push_back(2 * e);
            vhalf[head[e]].push_back(2 * e + 1);
        }
        for (auto &v : vhalf)
            sort(v.begin(), v.end());
        base.assign(V, -1);
        for (int v = 0; v < V; v++)
            if (!vhalf[v].empty())
                base[v] = vhalf[v][0];
    }
};

// Readable cell: vocc[v] = particle at vertex v or -1; vh[v] = half-edge h of its generator h - h0(v);
// el[e] = particles on edge e listed from tail to head.
struct Cell {
    vector<int> vocc, vh;
    vector<vector<int>> el;
};

// 128-bit key of a cell: particle p owns bits [12p, 12p+12): low 8 bits = half-edge (vertex particle)
// or 2E + e (edge particle), high 4 bits = rank on the edge.  Keys are compared as integers.
Key encode(const Graph &G, int N, const Cell &c) {
    Key s = 0;
    for (int v = 0; v < G.V; v++)
        if (c.vocc[v] >= 0)
            s |= (Key)(unsigned)c.vh[v] << (12 * c.vocc[v]);
    for (int e = 0; e < G.E; e++)
        for (int r = 0; r < (int)c.el[e].size(); r++)
            s |= (Key)(unsigned)((2 * G.E + e) | (r << 8)) << (12 * c.el[e][r]);
    return s;
}
// Inverse of encode: rebuild vertex occupation and the edge words (sorting particles by rank).
Cell decode(const Graph &G, int N, Key s) {
    Cell c;
    c.vocc.assign(G.V, -1);
    c.vh.assign(G.V, -1);
    c.el.assign(G.E, {});
    vector<vector<pair<int, int>>> tmp(G.E);
    for (int p = 0; p < N; p++) {
        unsigned w = (unsigned)((s >> (12 * p)) & 0xFFF);
        int loc = w & 0xFF, r = w >> 8;
        if (loc < 2 * G.E) {
            int v = G.hvert[loc];
            c.vocc[v] = p;
            c.vh[v] = loc;
        } else
            tmp[loc - 2 * G.E].push_back({r, p});
    }
    for (int e = 0; e < G.E; e++) {
        sort(tmp[e].begin(), tmp[e].end());
        for (auto &x : tmp[e])
            c.el[e].push_back(x.second);
    }
    return c;
}
// All k-cells for N particles, sorted.  Particles 0..N-1 are placed in turn, each either at a free
// vertex with a non-base half-edge (if fewer than k vertex particles so far) or at any position of any
// edge (if k vertex particles can still be reached).
void enumerate(const Graph &G, int N, int k, vector<Key> &out) {
    Cell c;
    c.vocc.assign(G.V, -1);
    c.vh.assign(G.V, -1);
    c.el.assign(G.E, {});
    int moving = 0;
    function<void(int)> rec = [&](int p) {
        if (p == N) {
            if (moving == k)
                out.push_back(encode(G, N, c));
            return;
        }
        if (moving < k)
            for (int v = 0; v < G.V; v++)
                if (c.vocc[v] < 0)
                    for (int h : G.vhalf[v])
                        if (h != G.base[v]) {
                            c.vocc[v] = p;
                            c.vh[v] = h;
                            moving++;
                            rec(p + 1);
                            moving--;
                            c.vocc[v] = -1;
                            c.vh[v] = -1;
                        }
        if (k - moving <= N - p - 1)
            for (int e = 0; e < G.E; e++) {
                int L = c.el[e].size();
                for (int pos = 0; pos <= L; pos++) {
                    c.el[e].insert(c.el[e].begin() + pos, p);
                    rec(p + 1);
                    c.el[e].erase(c.el[e].begin() + pos);
                }
            }
    };
    rec(0);
    sort(out.begin(), out.end());
}
// Boundary of a cell: for the i-th occupied vertex (increasing order) with particle p and half-edge h,
// sign (-1)^i times [p pushed into edge(h) at its h end] minus [p pushed into edge(h0) at its h0 end].
// Pushing in at a tail end puts p first in the edge word, at a head end last.
void boundary(const Graph &G, int N, const Cell &c, vector<pair<Key, int>> &out) {
    out.clear();
    int idx = 0;
    for (int v = 0; v < G.V; v++)
        if (c.vocc[v] >= 0) {
            int p = c.vocc[v], h = c.vh[v], h0 = G.base[v];
            int sgn = (idx % 2 == 0) ? 1 : -1;
            idx++;
            for (int t = 0; t < 2; t++) {
                int hh = t == 0 ? h : h0;
                Cell f = c;
                f.vocc[v] = -1;
                f.vh[v] = -1;
                int e = hh / 2;
                if (hh % 2 == 0)
                    f.el[e].insert(f.el[e].begin(), p);
                else
                    f.el[e].push_back(p);
                out.push_back({encode(G, N, f), t == 0 ? sgn : -sgn});
            }
        }
}
// Index of key k in the sorted list C, or -1.
int lookup(const vector<Key> &C, Key k) {
    auto it = lower_bound(C.begin(), C.end(), k);
    if (it == C.end() || *it != k)
        return -1;
    return it - C.begin();
}

// ---------------- F_p sparse elimination (rank) ----------------
u64 PR;
// Arithmetic modulo the global prime PR (128-bit intermediate products).
u64 mulm(u64 a, u64 b) {
    return (unsigned __int128)a * b % PR;
}
u64 powm(u64 a, u64 e) {
    u64 r = 1;
    while (e) {
        if (e & 1)
            r = mulm(r, a);
        a = mulm(a, a);
        e >>= 1;
    }
    return r;
}
// Sparse rank computation over F_PR: Markowitz-style pivoting (column with fewest entries, shortest row
// in it), exactly as SPE/SparseElim in engine/oconf.cpp; piv counts the pivots = rank.
struct FpElim {
    int nr, nc;
    vector<vector<pair<int, uint32_t>>> rows;
    vector<vector<int>> colrows;
    vector<int> colcnt;
    vector<char> rowdead, coldead;
    long long piv = 0;
    FpElim(int r, int c) : nr(r), nc(c), rows(r), colrows(c), colcnt(c, 0), rowdead(r, 0), coldead(c, 0) {}
    void finalize() {
        for (int r = 0; r < nr; r++) {
            auto &R = rows[r];
            sort(R.begin(), R.end());
            vector<pair<int, uint32_t>> M;
            for (auto &x : R) {
                if (!M.empty() && M.back().first == x.first)
                    M.back().second = (uint32_t)((M.back().second + (u64)x.second) % PR);
                else
                    M.push_back(x);
            }
            vector<pair<int, uint32_t>> M2;
            for (auto &x : M)
                if (x.second % PR)
                    M2.push_back({x.first, (uint32_t)(x.second % PR)});
            R.swap(M2);
            for (auto &x : R) {
                colrows[x.first].push_back(r);
                colcnt[x.first]++;
            }
        }
    }
    u64 get(int r, int c) const {
        auto &R = rows[r];
        auto it = lower_bound(R.begin(), R.end(), make_pair(c, (uint32_t)0));
        return (it != R.end() && it->first == c) ? it->second : 0;
    }
    void pivot(int r, int c) {
        auto &R = rows[r];
        u64 ui = powm(get(r, c), PR - 2);
        vector<int> tg;
        for (int r2 : colrows[c])
            if (r2 != r && !rowdead[r2])
                tg.push_back(r2);
        sort(tg.begin(), tg.end());
        tg.erase(unique(tg.begin(), tg.end()), tg.end());
        for (int r2 : tg) {
            u64 a = get(r2, c);
            if (!a)
                continue;
            u64 f = mulm(a, ui);
            auto &S = rows[r2];
            vector<pair<int, uint32_t>> out;
            out.reserve(S.size() + R.size());
            size_t i = 0, j = 0;
            while (i < S.size() || j < R.size()) {
                if (j == R.size() || (i < S.size() && S[i].first < R[j].first)) {
                    out.push_back(S[i]);
                    i++;
                } else if (i == S.size() || R[j].first < S[i].first) {
                    u64 v = (PR - mulm(f, R[j].second)) % PR;
                    if (v) {
                        out.push_back({R[j].first, (uint32_t)v});
                        colrows[R[j].first].push_back(r2);
                        colcnt[R[j].first]++;
                    }
                    j++;
                } else {
                    u64 v = (S[i].second + PR - mulm(f, R[j].second)) % PR;
                    if (v)
                        out.push_back({S[i].first, (uint32_t)v});
                    else
                        colcnt[S[i].first]--;
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
        piv++;
        vector<pair<int, uint32_t>>().swap(rows[r]);
        vector<int>().swap(colrows[c]);
    }
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
                for (auto &x : rows[best])
                    touched.push_back(x.first);
                pivot(best, c);
                any = true;
                for (int c2 : touched)
                    if (!coldead[c2] && colcnt[c2] > 0)
                        pq.push({colcnt[c2], c2});
            }
        }
    }
};

// star complex: cycle basis of Conf_k(Star_d) via BFS fundamental cycles
// Cycle basis of H_1(Conf_k(Star_d)).  The reduced complex of a star (one centre, d arms) has 0-cells
// (all particles on arms) and 1-cells (one particle at the centre), so it is a graph: each 1-cell joins
// its two faces.  cycles = fundamental cycles of a BFS spanning forest, as lists of (index of a 1-cell
// in C1, coefficient +-1).
struct StarBasis {
    Graph S;
    int k;
    vector<Key> C1;
    vector<vector<pair<int, int>>> cycles; // list of (index into C1, coefficient)
};
StarBasis star_basis(int d, int k) {
    StarBasis B;
    B.k = k;
    vector<pair<int, int>> ed;
    for (int i = 0; i < d; i++)
        ed.push_back({0, i + 1});
    B.S.build(d + 1, ed);
    vector<Key> C0;
    enumerate(B.S, k, 0, C0);
    enumerate(B.S, k, 1, B.C1);
    int n0 = C0.size(), n1 = B.C1.size();
    vector<int> hd(n1), tl(n1);
    vector<pair<Key, int>> bd;
    vector<vector<pair<int, int>>> adj(n0);
    for (int i = 0; i < n1; i++) {
        boundary(B.S, k, decode(B.S, k, B.C1[i]), bd);
        for (auto &x : bd) {
            int j = lookup(C0, x.first);
            if (x.second > 0)
                hd[i] = j;
            else
                tl[i] = j;
        }
        adj[hd[i]].push_back({i, tl[i]});
        adj[tl[i]].push_back({i, hd[i]});
    }
    vector<int> par(n0, -2), pedge(n0, -1), depth(n0, 0);
    vector<char> tree(n1, 0);
    // BFS over components (Conf of a star with d>=3 is connected)
    for (int s = 0; s < n0; s++)
        if (par[s] == -2) {
            deque<int> q;
            q.push_back(s);
            par[s] = -1;
            while (!q.empty()) {
                int x = q.front();
                q.pop_front();
                for (auto [ei, y] : adj[x])
                    if (par[y] == -2) {
                        par[y] = x;
                        pedge[y] = ei;
                        depth[y] = depth[x] + 1;
                        tree[ei] = 1;
                        q.push_back(y);
                    }
            }
        }
    for (int i = 0; i < n1; i++)
        if (!tree[i]) {
            vector<pair<int, int>> cyc;
            cyc.push_back({i, 1});
            // walk from hd[i] to tl[i] through the tree
            int a = hd[i], b = tl[i];
            vector<pair<int, int>> down; // edges from b side, to be traversed parent->child
            while (a != b) {
                if (depth[a] >= depth[b]) {
                    int e = pedge[a];
                    cyc.push_back({e, hd[e] == a ? -1 : 1});
                    a = par[a];
                } else {
                    int e = pedge[b];
                    down.push_back({e, hd[e] == b ? 1 : -1});
                    b = par[b];
                }
            }
            for (auto &x : down)
                cyc.push_back(x);
            B.cycles.push_back(cyc);
        }
    return B;
}

// main: (1) b2 = dim ker d_2 of Conf_n(G) mod PR; (2) cycle bases of the stars at u and w for every
// split k + (n-k) with 2 <= k <= n-2; (3) for every subset A of k particles and every pair of cycles, the
// product 2-chain (arm words of the two stars glued along the edges of G); (4) the rank of all products.
int main(int argc, char **argv) {
    int n = atoi(argv[1]);
    PR = strtoull(argv[2], 0, 10);
    int onlyA = argc > 3 ? atoi(argv[3]) : -1;
    int V, E;
    if (scanf("%d %d", &V, &E) != 2)
        return 1;
    vector<pair<int, int>> ed(E);
    for (auto &x : ed)
        if (scanf("%d %d", &x.first, &x.second) != 2)
            return 1;
    Graph G;
    G.build(V, ed);
    vector<int> ess;
    for (int v = 0; v < V; v++)
        if (G.vhalf[v].size() >= 2)
            ess.push_back(v);
    if (ess.size() != 2) {
        fprintf(stderr, "need exactly two non-leaf vertices\n");
        return 1;
    }
    int u = ess[0], w = ess[1];
    int du = G.vhalf[u].size(), dw = G.vhalf[w].size();
    vector<Key> C2, C1;
    enumerate(G, n, 2, C2);
    enumerate(G, n, 1, C1);
    // b2 = dim ker d_2 (top degree)
    long long rk2;
    {
        FpElim M(C2.size(), C1.size());
        vector<pair<Key, int>> bd;
        for (size_t i = 0; i < C2.size(); i++) {
            boundary(G, n, decode(G, n, C2[i]), bd);
            for (auto &x : bd)
                M.rows[i].push_back({lookup(C1, x.first), (uint32_t)((x.second + (ll)PR) % PR)});
        }
        M.finalize();
        M.run();
        rk2 = M.piv;
    }
    long long b2 = (long long)C2.size() - rk2;
    fprintf(stderr, "|C2|=%zu b2=%lld\n", C2.size(), b2);
    // star bases
    map<int, StarBasis> SBu, SBw;
    for (int k = 2; k <= n - 2; k++) {
        SBu[k] = star_basis(du, k);
        SBw[k] = star_basis(dw, k);
    }
    // half-edge -> arm index
    map<int, int> armU, armW;
    for (int i = 0; i < du; i++)
        armU[G.vhalf[u][i]] = i;
    for (int i = 0; i < dw; i++)
        armW[G.vhalf[w][i]] = i;
    // generators
    vector<vector<pair<int, uint32_t>>> gens;
    vector<int> gsize;
    long long bad = 0;
    for (int k = 2; k <= n - 2; k++) {
        if (onlyA >= 0 && k != onlyA)
            continue;
        StarBasis &Su = SBu[k];
        StarBasis &Sw = SBw[n - k];
        vector<Cell> cu, cw;
        for (Key x : Su.C1)
            cu.push_back(decode(Su.S, k, x));
        for (Key x : Sw.C1)
            cw.push_back(decode(Sw.S, n - k, x));
        vector<int> sel(n, 0);
        fill(sel.begin(), sel.begin() + k, 1);
        sort(sel.begin(), sel.end(), greater<int>());
        do {
            vector<int> A, Bv;
            for (int i = 0; i < n; i++)
                (sel[i] ? A : Bv).push_back(i);
            for (auto &z1 : Su.cycles)
                for (auto &z2 : Sw.cycles) {
                    vector<pair<int, uint32_t>> row;
                    for (auto [i1, s1] : z1)
                        for (auto [i2, s2] : z2) {
                            const Cell &a = cu[i1];
                            const Cell &b = cw[i2];
                            Cell c;
                            c.vocc.assign(V, -1);
                            c.vh.assign(V, -1);
                            c.el.assign(E, {});
                            c.vocc[u] = A[a.vocc[0]];
                            c.vh[u] = G.vhalf[u][a.vh[0] / 2];
                            c.vocc[w] = Bv[b.vocc[0]];
                            c.vh[w] = G.vhalf[w][b.vh[0] / 2];
                            for (int e = 0; e < E; e++) {
                                vector<int> word;
                                auto part = [&](int h) -> vector<int> {
                                    int vv = G.hvert[h];
                                    if (vv == u) {
                                        vector<int> r;
                                        for (int p : a.el[armU[h]])
                                            r.push_back(A[p]);
                                        return r;
                                    }
                                    if (vv == w) {
                                        vector<int> r;
                                        for (int p : b.el[armW[h]])
                                            r.push_back(Bv[p]);
                                        return r;
                                    }
                                    return {};
                                };
                                vector<int> tp = part(2 * e), hp = part(2 * e + 1);
                                word = tp;
                                for (int t = (int)hp.size() - 1; t >= 0; t--)
                                    word.push_back(hp[t]);
                                c.el[e] = word;
                            }
                            int j = lookup(C2, encode(G, n, c));
                            if (j < 0) {
                                bad++;
                                continue;
                            }
                            row.push_back({j, (uint32_t)(((ll)s1 * s2 + (ll)PR) % PR)});
                        }
                    gens.push_back(row);
                    gsize.push_back(k);
                }
        } while (prev_permutation(sel.begin(), sel.end()));
    }
    if (bad) {
        fprintf(stderr, "BUG: %lld product cells not found\n", bad);
        return 2;
    }
    // sanity: each generator is a cycle (check a sample)
    {
        vector<pair<Key, int>> bd;
        mt19937 rng(1);
        int checks = min<size_t>(200, gens.size());
        for (int t = 0; t < checks; t++) {
            auto &g = gens[rng() % gens.size()];
            map<int, ll> acc;
            for (auto [j, cf] : g) {
                boundary(G, n, decode(G, n, C2[j]), bd);
                ll cc = cf > PR / 2 ? (ll)cf - (ll)PR : cf;
                for (auto &x : bd)
                    acc[lookup(C1, x.first)] += cc * x.second;
            }
            for (auto &kv : acc)
                if (kv.second % (ll)PR) {
                    fprintf(stderr, "BUG: product not a cycle\n");
                    return 3;
                }
        }
    }
    long long nnz = 0;
    for (auto &g : gens)
        nnz += g.size();
    fprintf(stderr, "generators=%zu nnz=%lld\n", gens.size(), nnz);
    FpElim M(gens.size(), C2.size());
    for (size_t i = 0; i < gens.size(); i++)
        M.rows[i].swap(gens[i]);
    vector<vector<pair<int, uint32_t>>>().swap(gens);
    M.finalize();
    M.run();
    printf(
        "{\"n\":%d,\"prime\":%llu,\"du\":%d,\"dw\":%d,\"C2\":%zu,\"b2\":%lld,\"gens\":%zu,\"rankP\":%lld}\n",
        n, PR, du, dw, C2.size(), b2, gsize.size(), M.piv);
    return 0;
}
