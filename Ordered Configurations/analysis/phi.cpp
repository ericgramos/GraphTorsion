// =====================================================================================================
//  phi.cpp  --  rank of the "double push" map C_3 -> C_1 for a graph with three essential vertices
// =====================================================================================================
//  For each 3-cell c and each pair of its three occupied vertices x < y, the particles at x and at y are
//  both pushed into an edge (each either along its own half-edge h, or along the base half-edge h0, with a
//  sign -1 for every h0 push), leaving a 1-cell.  The map Phi: C_3 -> C_1 is the sum of these terms; its
//  rank modulo the given prime is printed.  This auxiliary rank enters the bookkeeping of the free part of
//  H_2 for T(2,2,2) (how much of b_2 is spanned by star x star classes): the star x star span has dimension
//  sum_x b_2(Conf_n Gamma_x) - (rank d_3 - rank Phi); see the header of decomp.cpp.
//
//  OUTPUT.  {"C3","C1","rankPhi"}.     USAGE.  phi n prime < graph
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
    // Set up the half-edge tables (half-edge 2e = tail end of edge e, 2e+1 = head end; vhalf[v] sorted) and the base half-edge of each vertex (its lowest half-edge).
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

// Pack a cell into a 128-bit key: 12 bits per particle, the low 8 its location (half-edge if at a vertex, 2E + e if on edge e) and the high 4 its rank in the word of that edge.
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
// Inverse of encode: rebuild the occupants of the vertices, their half-edges and the edge words from a key.
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
// All cells of degree k (k occupied vertices), placing the particles one at a time at a vertex (with a non-base half-edge) or into an edge word.
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
// Signed boundary of a cell: for the j-th occupied vertex v carrying (p, h), (-1)^j ( [p pushed into edge(h)] - [p pushed into the base edge] ).
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
// Index of a key in the sorted list of cells of its degree (binary search); -1 if absent.
int lookup(const vector<Key> &C, Key k) {
    auto it = lower_bound(C.begin(), C.end(), k);
    if (it == C.end() || *it != k)
        return -1;
    return it - C.begin();
}

// ---------------- F_p sparse elimination (rank) ----------------
u64 PR;
u64 mulm(u64 a, u64 b) {
    return (unsigned __int128)a * b % PR;
}
// a^e modulo the prime, by repeated squaring.
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
struct FpElim {
    int nr, nc;
    vector<vector<pair<int, uint32_t>>> rows;
    vector<vector<int>> colrows;
    vector<int> colcnt;
    vector<char> rowdead, coldead;
    long long piv = 0;
    FpElim(int r, int c) : nr(r), nc(c), rows(r), colrows(c), colcnt(c, 0), rowdead(r, 0), coldead(c, 0) {}
    // Sort every sparse row, merge repeated columns, drop zero entries, and build the column structures (rows per column, number of usable pivot entries per column).
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
    // Entry (r, c) of the sparse matrix, found by binary search in row r; 0 if absent.
    u64 get(int r, int c) const {
        auto &R = rows[r];
        auto it = lower_bound(R.begin(), R.end(), make_pair(c, (uint32_t)0));
        return (it != R.end() && it->first == c) ? it->second : 0;
    }
    // Pivot on the usable entry (r, c): subtract the right multiple of row r from every other live row with an entry in column c, then mark row r and column c as dead.
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

// rank of Phi: C_3 -> C_1, w |-> sum over vertex pairs {x<v} of d_v d_x w  (graphs with 3 essential vertices)
// rank(V2) = sum_x b2(G \ x) + rank(Phi) - rank(d_3), V2 = sum_x Im H_2(Conf(G \ x)).
int main(int argc, char **argv) {
    int n = atoi(argv[1]);
    PR = strtoull(argv[2], 0, 10);
    int V, E;
    if (scanf("%d %d", &V, &E) != 2)
        return 1;
    vector<pair<int, int>> ed(E);
    for (auto &x : ed)
        if (scanf("%d %d", &x.first, &x.second) != 2)
            return 1;
    Graph G;
    G.build(V, ed);
    vector<Key> C3, C1;
    enumerate(G, n, 3, C3);
    enumerate(G, n, 1, C1);
    FpElim M(C3.size(), C1.size());
    for (size_t i = 0; i < C3.size(); i++) {
        Cell c = decode(G, n, C3[i]);
        vector<int> occ;
        for (int v = 0; v < V; v++)
            if (c.vocc[v] >= 0)
                occ.push_back(v);
        for (int a = 0; a < 3; a++)
            for (int b = a + 1; b < 3; b++) {
                int x = occ[a], y = occ[b];
                for (int t1 = 0; t1 < 2; t1++)
                    for (int t2 = 0; t2 < 2; t2++) {
                        Cell f = c;
                        int sg = 1;
                        for (int step = 0; step < 2; step++) {
                            int v = step == 0 ? x : y;
                            int t = step == 0 ? t1 : t2;
                            int p = f.vocc[v], hh = t == 0 ? f.vh[v] : G.base[v];
                            f.vocc[v] = -1;
                            f.vh[v] = -1;
                            int e = hh / 2;
                            if (hh % 2 == 0)
                                f.el[e].insert(f.el[e].begin(), p);
                            else
                                f.el[e].push_back(p);
                            if (t == 1)
                                sg = -sg;
                        }
                        int j = lookup(C1, encode(G, n, f));
                        M.rows[i].push_back({j, (uint32_t)((sg + (ll)PR) % PR)});
                    }
            }
    }
    M.finalize();
    M.run();
    printf("{\"C3\":%zu,\"C1\":%zu,\"rankPhi\":%lld}\n", C3.size(), C1.size(), M.piv);
}
