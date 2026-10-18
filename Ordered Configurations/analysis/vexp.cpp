// =====================================================================================================
//  vexp.cpp  --  p-torsion of H_k(Conf_n Gamma) via the ordered vertex explosion at one vertex v
// =====================================================================================================
//  MATHEMATICS.  Let Gamma_v be Gamma with v exploded: every half-edge at v becomes the
//  end of a new leaf; edge numbers and orientations are kept.  The reduced complex of Gamma is the mapping
//  cone of
//        Delta = sum_{h != h0} (sigma_h^q - sigma_{h0}^q) :  (+)_{q, h != h0} S_{[n]\q}(Gamma_v) -> S_n(Gamma_v),
//  sigma_h^q inserting particle q at the leaf end of edge(h).  If H_{k-1}(Conf_{n-1} Gamma_v) has no
//  p-torsion, the long exact sequence gives
//        Tors_p H_k(Conf_n Gamma) = Tors_p( C_k(Gamma_v; n) / ( B_k(Gamma_v; n) + Delta(Z_k(Gamma_v; n-1)^(+)) ) ).
//
//  ALGORITHM.
//   (n-1)-particle side:  Pi = pivot columns (k-cells) of d_{k+1} of Gamma_v mod p (a complement of the
//     boundaries); cycle representatives = kernel vectors of d_k restricted to the k-cells outside Pi, one
//     per free column (LPE with keep = true, kernel_unit back-substitution).  The hypothesis "no p-torsion
//     in H_{k-1}(n-1)" is checked: all pivots of that elimination must be units (else exit 3).  If there
//     are no cycle representatives, nothing is added and the torsion is zero.
//   n-particle side:  the rows of d_{k+1}(Gamma_v; n) together with the rows Delta(z) for every cycle
//     representative z, particle q (relabelling the other particles into [n] \ q) and h != h0, all as
//     vectors over the k-cells of Gamma_v at n particles (k-cells that are unit pivot rows of d_k are dropped
//     from the columns, as in the engine).  The local Smith form of this stacked matrix modulo p^K gives
//     the p-primary torsion (stages > 0).
//
//  OUTPUT.  {"n","v","k","p","K","stages":[...],"leftover_nnz",...}.   Stage s > 0 counts Z/p^s summands.
//  USAGE.   vexp n v k p K < graph          (p^K <= 31)
//  VALIDATION.  validation/vexp_validate.py compares with direct computation (vexp_validate*.log).
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

struct Cell {
    vector<int> vocc, vh;
    vector<vector<int>> el;
};

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
int lookup(const vector<Key> &C, Key k) {
    auto it = lower_bound(C.begin(), C.end(), k);
    if (it == C.end() || *it != k)
        return -1;
    return it - C.begin();
}

struct LPE {
    int nr, nc;
    uint32_t p0, mod;
    vector<vector<uint32_t>> rows;
    vector<vector<int>> colrows;
    vector<int> colcnt;
    vector<char> rowdead, coldead;
    long long pivots = 0;
    ll maxabs = 1;
    vector<long long> stages;
    long long leftover_nonzero = 0;
    bool keep = false;
    vector<pair<int, vector<uint32_t>>> piv;
    LPE(int r, int c, uint32_t p, uint32_t m)
        : nr(r), nc(c), p0(p), mod(m), rows(r), colrows(c), colcnt(c, 0), rowdead(r, 0), coldead(c, 0) {}
    static uint32_t col(uint32_t x) {
        return x >> 5;
    }
    static uint32_t val(uint32_t x) {
        return x & 31;
    }
    uint32_t inv(uint32_t a) const {
        for (uint32_t b = 1; b < mod; b++)
            if (a * b % mod == 1)
                return b;
        return 0;
    }
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
    uint32_t get(int r, int c) const {
        auto &R = rows[r];
        auto it = lower_bound(R.begin(), R.end(), (uint32_t)c << 5);
        return (it != R.end() && col(*it) == (uint32_t)c) ? val(*it) : 0;
    }
    bool unit(uint32_t v) const {
        return v % p0 != 0;
    }
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
            S.swap(out);
        }
        for (auto x : R)
            colcnt[col(x)]--;
        rowdead[r] = 1;
        coldead[c] = 1;
        colcnt[c] = 0;
        pivots++;
        if (keep)
            piv.push_back({c, rows[r]});
        vector<uint32_t>().swap(rows[r]);
        vector<int>().swap(colrows[c]);
    }
    bool sweep() {
        bool any = false;
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
            for (int c2 : touched)
                if (!coldead[c2] && colcnt[c2] > 0)
                    pq.push({colcnt[c2], c2});
        }
        return any;
    }
    // kernel vector with free column f set to 1 (requires all pivots units; mod = original modulus)
    vector<pair<int, uint32_t>> kernel_unit(int f, uint32_t M) {
        unordered_map<int, uint32_t> x;
        x[f] = 1;
        for (int i = (int)piv.size() - 1; i >= 0; i--) {
            int c = piv[i].first;
            uint32_t s2 = 0, ac = 0;
            for (auto y : piv[i].second) {
                int cc = col(y);
                if (cc == c)
                    ac = val(y);
                else {
                    auto it = x.find(cc);
                    if (it != x.end())
                        s2 = (s2 + val(y) * it->second) % M;
                }
            }
            if (s2) {
                uint32_t ia = 0;
                for (uint32_t t = 1; t < M; t++)
                    if (ac * t % M == 1)
                        ia = t;
                x[c] = (M - s2) % M * ia % M;
            }
        }
        vector<pair<int, uint32_t>> out;
        for (auto &kv : x)
            if (kv.second)
                out.push_back(kv);
        sort(out.begin(), out.end());
        return out;
    }
    void run() {
        long long before = 0;
        while (true) {
            while (sweep()) {
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
            // all live entries divisible by p0: divide out
            for (int r = 0; r < nr; r++)
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

// vexp: p-torsion of H_k(Conf_n(Gamma)) via the ordered vertex explosion at v.
//   Tors_p H_k(Gamma; n) = Tors_p( C_k(Gamma_v; n) / ( B_k(Gamma_v; n) + Delta( Z_k(Gamma_v; n-1)^{(+) q,h} ) ) )
// valid when H_{k-1}(Conf_{n-1}(Gamma_v)) has no p-torsion (checked: all pivots of d_k(Gamma_v; n-1) are units).
// usage: vexp n v k p K < graph            (works modulo p^K, p^K <= 31)
struct Timer {
    chrono::steady_clock::time_point t = chrono::steady_clock::now();
    double s() {
        return chrono::duration<double>(chrono::steady_clock::now() - t).count();
    }
};
int main(int argc, char **argv) {
    int n = atoi(argv[1]), v = atoi(argv[2]), k = atoi(argv[3]);
    uint32_t p = atoi(argv[4]), K = atoi(argv[5]);
    uint32_t M = 1;
    for (uint32_t i = 0; i < K; i++)
        M *= p;
    if (M > 31) {
        fprintf(stderr, "p^K must be <= 31\n");
        return 1;
    }
    int V, E;
    if (scanf("%d %d", &V, &E) != 2)
        return 1;
    vector<pair<int, int>> ed(E);
    for (auto &x : ed)
        if (scanf("%d %d", &x.first, &x.second) != 2)
            return 1;
    Graph G;
    G.build(V, ed);
    int d = G.vhalf[v].size();
    int h0 = G.base[v];
    vector<int> hv;
    for (int h : G.vhalf[v])
        if (h != h0)
            hv.push_back(h);
    // Gamma_v: same edge indices/orientations, each end at v replaced by a fresh leaf
    vector<pair<int, int>> ex = ed;
    int nv = V;
    for (auto &x : ex) {
        if (x.first == v)
            x.first = nv++;
        if (x.second == v)
            x.second = nv++;
    }
    Graph X;
    X.build(nv, ex);
    Timer T;
    auto bmat = [&](int m, int deg, vector<Key> &Cd, vector<Key> &Cdm1, LPE &L, bool transpose,
                    const vector<char> *dropcol) {
        // rows = deg-cells (or transposed), entries mod M
        vector<pair<Key, int>> bd;
        for (size_t i = 0; i < Cd.size(); i++) {
            boundary(X, m, decode(X, m, Cd[i]), bd);
            for (auto &y : bd) {
                int j = lookup(Cdm1, y.first);
                if (dropcol && (*dropcol)[j])
                    continue;
                uint32_t val = (uint32_t)((y.second % (int)M + (int)M) % M);
                if (transpose)
                    L.rows[j].push_back((uint32_t)i << 5 | val);
                else
                    L.rows[i].push_back((uint32_t)j << 5 | val);
            }
        }
    };
    // ---------- (n-1)-particle side ----------
    int m = n - 1;
    vector<Key> A1, A0, Am;
    enumerate(X, m, k + 1, A1);
    enumerate(X, m, k, A0);
    if (k >= 1)
        enumerate(X, m, k - 1, Am);
    // Pi: pivot columns (k-cells) of d_{k+1} mod p
    vector<char> inPi(A0.size(), 0);
    {
        LPE L(A1.size(), A0.size(), p, p);
        bmat(m, k + 1, A1, A0, L, false, nullptr);
        L.finalize();
        L.run();
        for (size_t c = 0; c < A0.size(); c++)
            if (L.coldead[c])
                inPi[c] = 1;
    }
    // kernel of d_k restricted to columns outside Pi: transpose system rows=(k-1)-cells, cols = all k-cells, Pi columns zeroed
    vector<vector<pair<int, uint32_t>>> Z;
    bool hyp = true;
    {
        LPE L(k >= 1 ? Am.size() : 0, A0.size(), p, M);
        L.keep = true;
        if (k >= 1) {
            vector<pair<Key, int>> bd;
            for (size_t i = 0; i < A0.size(); i++) {
                if (inPi[i])
                    continue;
                boundary(X, m, decode(X, m, A0[i]), bd);
                for (auto &y : bd) {
                    int j = lookup(Am, y.first);
                    L.rows[j].push_back((uint32_t)i << 5 | (uint32_t)((y.second % (int)M + (int)M) % M));
                }
            }
        }
        L.finalize();
        L.run();
        if (L.stages.size() > 1 || L.leftover_nonzero)
            hyp = false;
        for (size_t f = 0; f < A0.size(); f++)
            if (!inPi[f] && !L.coldead[f])
                Z.push_back(L.kernel_unit(f, M));
    }
    long long znnz = 0;
    for (auto &z : Z)
        znnz += z.size();
    fprintf(stderr, "[%.1fs] n-1 side: |C_k|=%zu |Pi|=%ld cycle reps=%zu (nnz %lld) hypothesis %s\n", T.s(),
            A0.size(), (long)count(inPi.begin(), inPi.end(), 1), Z.size(), znnz, hyp ? "ok" : "FAILS");
    if (!hyp) {
        printf("{\"hypothesis\":false}\n");
        return 3;
    }
    if (Z.empty()) { // no relations: the image of H_k(Gamma_v; n) is a subgroup of the free group Z_k, so no torsion
        printf("{\"n\":%d,\"v\":%d,\"k\":%d,\"p\":%u,\"K\":%u,\"stages\":[0],\"leftover_nnz\":0,\"no_"
               "relations\":true,\"secs\":%.1f}\n",
               n, v, k, p, K, T.s());
        return 0;
    }
    vector<Cell> Zc;
    for (Key key : A0)
        Zc.push_back(decode(X, m, key));
    // ---------- n-particle side ----------
    vector<Key> B1, B0, Bm;
    enumerate(X, n, k + 1, B1);
    enumerate(X, n, k, B0);
    if (k >= 1)
        enumerate(X, n, k - 1, Bm);
    // column dropping: k-cells that are unit pivot rows of d_k(Gamma_v; n)
    vector<char> drop(B0.size(), 0);
    bool canDrop = false;
    if (k >= 1) {
        LPE L(B0.size(), Bm.size(), p, M);
        bmat(n, k, B0, Bm, L, false, nullptr);
        L.finalize();
        L.run();
        if (L.stages.size() <= 1 && !L.leftover_nonzero) {
            canDrop = true;
            for (size_t i = 0; i < B0.size(); i++)
                if (L.rowdead[i])
                    drop[i] = 1;
        }
    }
    size_t nrowsM = Z.size() * n * hv.size();
    LPE L(B1.size() + nrowsM, B0.size(), p, M);
    bmat(n, k + 1, B1, B0, L, false, canDrop ? &drop : nullptr);
    size_t r = B1.size();
    for (int q = 0; q < n; q++) {
        auto rel = [&](int j) {
            return j < q ? j : j + 1;
        };
        for (auto &z : Z) {
            for (int h : hv) {
                auto &row = L.rows[r++];
                for (auto [ci, cf] : z) {
                    const Cell &c0 = Zc[ci];
                    for (int t = 0; t < 2; t++) {
                        int hh = t == 0 ? h : h0;
                        Cell c;
                        c.vocc = c0.vocc;
                        c.vh = c0.vh;
                        c.el = c0.el;
                        for (auto &x : c.vocc)
                            if (x >= 0)
                                x = rel(x);
                        for (auto &w : c.el)
                            for (auto &x : w)
                                x = rel(x);
                        int e = hh / 2;
                        if (hh % 2 == 0)
                            c.el[e].insert(c.el[e].begin(), q);
                        else
                            c.el[e].push_back(q);
                        int j = lookup(B0, encode(X, n, c));
                        if (j < 0) {
                            fprintf(stderr, "BUG: stabilized cell missing\n");
                            return 2;
                        }
                        if (canDrop && drop[j])
                            continue;
                        uint32_t val = t == 0 ? cf : (M - cf) % M;
                        row.push_back((uint32_t)j << 5 | val);
                    }
                }
            }
        }
    }
    fprintf(stderr, "[%.1fs] n side: |C_{k+1}|=%zu |C_k|=%zu dropped=%ld Delta rows=%zu\n", T.s(), B1.size(),
            B0.size(), (long)count(drop.begin(), drop.end(), 1), nrowsM);
    L.finalize();
    L.run();
    printf("{\"n\":%d,\"v\":%d,\"k\":%d,\"p\":%u,\"K\":%u,\"stages\":[", n, v, k, p, K);
    for (size_t i = 0; i < L.stages.size(); i++)
        printf("%s%lld", i ? "," : "", L.stages[i]);
    printf("],\"leftover_nnz\":%lld,\"secs\":%.1f}\n", L.leftover_nonzero, T.s());
}
