// =====================================================================================================
//  g2tors.cpp  --  the "triangle shortcut": torsion of H_2 for a graph that becomes two-vertex after
//                  deleting one vertex x
// =====================================================================================================
//  MATHEMATICS.  Let Gamma be a graph and x a vertex such that Gamma \ x (edges at x turned into hairs
//  ending at new leaves) has exactly two essential vertices u, w.  By the ordered vertex explosion at x
//  (the complex of Gamma is the mapping cone of Delta, see vexp.cpp) and because H_1(Conf_{n-1}(Gamma \ x))
//  is torsion free,
//        Tors H_2(Conf_n Gamma) = Tors( H_2(Conf_n(Gamma \ x)) / < sigma_h^q y - sigma_{h0}^q y > ),
//  where q runs over particles, h over the non-base half-edges at x, y over the (n-1)-particle star x star
//  tori of Gamma \ x (particles [n] \ q), and sigma_h^q places particle q at the x-end of the hair h.
//  Since H_2(Conf_n(Gamma \ x)) = Z_2 is the top homology (a subgroup of the free group C_2), the torsion
//  of this quotient is the failure of the relation lattice M = span{sigma_h^q y - sigma_{h0}^q y} to be
//  saturated in C_2.  Its p-part has dimension  rank_Q(M) - rank_{F_p}(M).
//
//  WHAT THE PROGRAM DOES.  Builds all relation vectors (one row per (q, torus y, h != h0)) as sparse rows
//  over the cells C_2 of Conf_n(Gamma \ x) and prints their rank modulo the given prime.  Running it with a
//  large prime (1000003) and with p gives the number of Z/p^? summands (the difference of the ranks).
//
//  OUTPUT.  {"n","x","prime","C2","rankM"}.
//  USAGE.   g2tors n x prime < graph
//  RESULTS.  data/survey/triangle_family_shortcut.jsonl (T(a,b,c), n = 5, 6); it reproduced the full-complex
//  torsion of T(1,2,3), T(2,2,2), T(1,3,3) and found 10 summands of 2-power order and 2 of 3-power order
//  for T(1,2,4) (local Smith forms of the full complex show that all of them have order exactly 2 or 3).
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

// ---------------- F_p sparse elimination (rank) ----------------
u64 PR;
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

// g2tors: p-torsion of H_2(Conf_n(Gamma)) for Gamma with 3 essential vertices, computed as
//   Tors( C_2(Gamma \ x) / M ),  M = span{ sigma_h^q(y) - sigma_h0^q(y) },
// y = star-product tori in Conf_{n-1}(Gamma \ x) on the particles other than q, h = half-edges at x.
// Valid when (i) H_1(Conf_{n-1}(Gamma\x)) is torsion-free and (ii) star products span
// H_2(Conf_{n-1}(Gamma\x)) over Z_(p) (checked separately with prod2).
// Also: optional split filter for the n-particle side is not used here.
// usage: g2tors n x prime < Gamma        prints rank of M mod prime
int main(int argc, char **argv) {
    int n = atoi(argv[1]);
    int X = atoi(argv[2]);
    PR = strtoull(argv[3], 0, 10);
    int V, E;
    if (scanf("%d %d", &V, &E) != 2)
        return 1;
    vector<pair<int, int>> ged(E);
    for (auto &x : ged)
        if (scanf("%d %d", &x.first, &x.second) != 2)
            return 1;
    // build G2 = Gamma \ X: keep vertex numbering of non-X vertices, new leaf per edge end at X
    vector<int> vmap(V, -1);
    int nv = 0;
    for (int v = 0; v < V; v++)
        if (v != X)
            vmap[v] = nv++;
    vector<pair<int, int>> ed;
    vector<int> hairEnd(E, -1); // for edges at X: 0 = X is tail, 1 = X is head
    for (int e = 0; e < E; e++) {
        int a = ged[e].first, b = ged[e].second;
        if (a == X && b == X) {
            fprintf(stderr, "loop at x unsupported\n");
            return 1;
        }
        if (a == X) {
            ed.push_back({nv++, vmap[b]});
            hairEnd[e] = 0;
        } else if (b == X) {
            ed.push_back({vmap[a], nv++});
            hairEnd[e] = 1;
        } else
            ed.push_back({vmap[a], vmap[b]});
    }
    Graph G;
    G.build(nv, ed);
    vector<int> ess;
    for (int v = 0; v < G.V; v++)
        if (G.vhalf[v].size() >= 2)
            ess.push_back(v);
    if (ess.size() != 2) {
        fprintf(stderr, "Gamma minus x must have two essential vertices\n");
        return 1;
    }
    int u = ess[0], w = ess[1];
    int du = G.vhalf[u].size(), dw = G.vhalf[w].size();
    // half-edges at X in Gamma (sorted), base first
    vector<int> xedges;
    for (int e = 0; e < E; e++)
        if (hairEnd[e] >= 0)
            xedges.push_back(e); // edge order = half-edge order
    vector<Key> C2;
    enumerate(G, n, 2, C2);
    map<int, StarBasis> SBu, SBw;
    for (int k = 2; k <= n - 3; k++) {
        SBu[k] = star_basis(du, k);
        SBw[k] = star_basis(dw, k);
    }
    map<int, int> armU, armW;
    for (int i = 0; i < du; i++)
        armU[G.vhalf[u][i]] = i;
    for (int i = 0; i < dw; i++)
        armW[G.vhalf[w][i]] = i;
    vector<vector<pair<int, uint32_t>>> rows;
    long long missing = 0;
    for (int q = 0; q < n; q++) {
        vector<int> S;
        for (int i = 0; i < n; i++)
            if (i != q)
                S.push_back(i);
        int m = n - 1;
        for (int k = 2; k <= m - 2; k++) {
            StarBasis &Su = SBu[k];
            StarBasis &Sw = SBw[m - k];
            vector<Cell> cu, cw;
            for (Key t : Su.C1)
                cu.push_back(decode(Su.S, k, t));
            for (Key t : Sw.C1)
                cw.push_back(decode(Sw.S, m - k, t));
            vector<int> sel(m, 0);
            fill(sel.begin(), sel.begin() + k, 1);
            sort(sel.begin(), sel.end(), greater<int>());
            do {
                vector<int> A, B;
                for (int i = 0; i < m; i++)
                    (sel[i] ? A : B).push_back(S[i]);
                for (auto &z1 : Su.cycles)
                    for (auto &z2 : Sw.cycles) {
                        // the torus y as a list of (cell, coef)
                        vector<pair<Cell, int>> y;
                        for (auto [i1, s1] : z1)
                            for (auto [i2, s2] : z2) {
                                const Cell &a = cu[i1];
                                const Cell &b = cw[i2];
                                Cell c;
                                c.vocc.assign(G.V, -1);
                                c.vh.assign(G.V, -1);
                                c.el.assign(G.E, {});
                                c.vocc[u] = A[a.vocc[0]];
                                c.vh[u] = G.vhalf[u][a.vh[0] / 2];
                                c.vocc[w] = B[b.vocc[0]];
                                c.vh[w] = G.vhalf[w][b.vh[0] / 2];
                                for (int e = 0; e < G.E; e++) {
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
                                                r.push_back(B[p]);
                                            return r;
                                        }
                                        return {};
                                    };
                                    vector<int> tp = part(2 * e), hp = part(2 * e + 1), word = tp;
                                    for (int t = (int)hp.size() - 1; t >= 0; t--)
                                        word.push_back(hp[t]);
                                    c.el[e] = word;
                                }
                                y.push_back({c, s1 * s2});
                            }
                        auto stab = [&](const Cell &c, int e) {
                            Cell d = c;
                            if (hairEnd[e] == 0)
                                d.el[e].insert(d.el[e].begin(), q);
                            else
                                d.el[e].push_back(q);
                            return d;
                        };
                        for (size_t hi = 1; hi < xedges.size(); hi++) {
                            vector<pair<int, uint32_t>> row;
                            for (auto &[c, cf] : y) {
                                int j1 = lookup(C2, encode(G, n, stab(c, xedges[hi])));
                                int j0 = lookup(C2, encode(G, n, stab(c, xedges[0])));
                                if (j1 < 0 || j0 < 0) {
                                    missing++;
                                    continue;
                                }
                                row.push_back({j1, (uint32_t)((cf + (ll)PR) % PR)});
                                row.push_back({j0, (uint32_t)((-cf + (ll)PR) % PR)});
                            }
                            rows.push_back(row);
                        }
                    }
            } while (prev_permutation(sel.begin(), sel.end()));
        }
    }
    if (missing) {
        fprintf(stderr, "BUG: %lld cells missing\n", missing);
        return 2;
    }
    long long nnz = 0;
    for (auto &r : rows)
        nnz += r.size();
    fprintf(stderr, "C2=%zu rows=%zu nnz=%lld\n", C2.size(), rows.size(), nnz);
    FpElim M(rows.size(), C2.size());
    for (size_t i = 0; i < rows.size(); i++)
        M.rows[i].swap(rows[i]);
    vector<vector<pair<int, uint32_t>>>().swap(rows);
    M.finalize();
    M.run();
    printf("{\"n\":%d,\"x\":%d,\"prime\":%llu,\"C2\":%zu,\"rankM\":%lld}\n", n, X, PR, C2.size(), M.piv);
}
