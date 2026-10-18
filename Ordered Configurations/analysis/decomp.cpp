// =====================================================================================================
//  decomp.cpp  --  how much of H_2(Conf_n Gamma) is spanned by tori (products of disjoint 1-classes)?
// =====================================================================================================
//  For every essential vertex x of Gamma, every k (particles at x) and every subset A of k particles, and
//  every pair (z_x, z_P) with z_x a cycle of Conf_A(Star_x) and z_P a cycle of Conf_B(Gamma \ x)
//  (B = [n] \ A; z1_basis computes a cycle basis of the 1-dimensional complex of the piece), the product
//  z_x * z_P is written as a 2-chain of Gamma.  The sign of the product includes the Koszul sign of the
//  vertex generator at x relative to the occupied vertex of the piece (kos = +1 if that vertex is smaller
//  than x, else -1).  A sample of generators is checked to be cycles (exit 3 otherwise).
//
//  It prints rankB2 = rank of the boundaries d_3(C_3) and rank_span_mod_B = rank of (tori + boundaries)
//  minus rankB2, i.e. the dimension of the span of the tori inside H_2 (mod prime).
//  mode 1: keep only those BASIS cycles of the piece whose cells all have the same occupied vertex.  This
//          filters a fundamental-cycle basis, so it spans only part of the star x star tori and gives a LOWER
//          BOUND for their span (for T(2,2,2), n = 4 it gives 210, while the star x star span is 448).  The
//          exact star x star span is  sum_x b_2(Conf_n Gamma_x) - (rank d_3 - rank Phi),  with b_2(Gamma_x) from
//          prod2 and rank Phi from phi.  Reason: the 2-cells with x unoccupied are disjoint for different x, so the
//          span is (+)_x Z_2(Gamma_x) modulo its intersection with the boundaries B_2, and that intersection has
//          dimension rank d_3 - rank Phi.
//
//  OUTPUT.  {"n","prime","mode","gens","rankB2","rank_span_mod_B"}.    USAGE.  decomp n prime [mode] < graph
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

// fundamental-cycle basis of the 1-skeleton (cells of degree <= 1) of Conf_k(G): Z_1 spanning set
vector<vector<pair<int, int>>> z1_basis(const Graph &G, int k, vector<Key> &C1) {
    vector<Key> C0;
    enumerate(G, k, 0, C0);
    C1.clear();
    enumerate(G, k, 1, C1);
    int n0 = C0.size(), n1 = C1.size();
    vector<int> hd(n1), tl(n1);
    vector<pair<Key, int>> bd;
    vector<vector<pair<int, int>>> adj(n0);
    for (int i = 0; i < n1; i++) {
        boundary(G, k, decode(G, k, C1[i]), bd);
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
    for (int s0 = 0; s0 < n0; s0++)
        if (par[s0] == -2) {
            deque<int> q;
            q.push_back(s0);
            par[s0] = -1;
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
    vector<vector<pair<int, int>>> cyc;
    for (int i = 0; i < n1; i++)
        if (!tree[i]) {
            vector<pair<int, int>> c;
            c.push_back({i, 1});
            int a = hd[i], b = tl[i];
            vector<pair<int, int>> down;
            while (a != b) {
                if (depth[a] >= depth[b]) {
                    int e = pedge[a];
                    c.push_back({e, hd[e] == a ? -1 : 1});
                    a = par[a];
                } else {
                    int e = pedge[b];
                    down.push_back({e, hd[e] == b ? 1 : -1});
                    b = par[b];
                }
            }
            for (auto &x : down)
                c.push_back(x);
            cyc.push_back(c);
        }
    return cyc;
}

// decomp: rank in H_2(Conf_n(G)) of span{ star_x (x) Z_1(Conf(G cut at the edges of x)) : x essential }
// plus optional star x star only mode.   usage: decomp n prime mode < graph
//   mode 0: x ranges over all essential vertices, second factor = all 1-cycles of the piece
//   mode 1: second factor restricted to basis cycles supported at one vertex (a lower bound, see the header)
int main(int argc, char **argv) {
    int n = atoi(argv[1]);
    PR = strtoull(argv[2], 0, 10);
    int mode = argc > 3 ? atoi(argv[3]) : 0;
    int V, E;
    if (scanf("%d %d", &V, &E) != 2)
        return 1;
    vector<pair<int, int>> ed(E);
    for (auto &x : ed)
        if (scanf("%d %d", &x.first, &x.second) != 2)
            return 1;
    Graph G;
    G.build(V, ed);
    vector<Key> C3, C2;
    enumerate(G, n, 3, C3);
    enumerate(G, n, 2, C2);
    vector<vector<pair<int, uint32_t>>> rows;
    {
        vector<pair<Key, int>> bd;
        for (size_t i = 0; i < C3.size(); i++) {
            boundary(G, n, decode(G, n, C3[i]), bd);
            vector<pair<int, uint32_t>> r;
            for (auto &x : bd)
                r.push_back({lookup(C2, x.first), (uint32_t)((x.second + (ll)PR) % PR)});
            rows.push_back(r);
        }
    }
    size_t nb = rows.size();
    long long ngen = 0;
    vector<int> essv;
    for (int v = 0; v < V; v++)
        if (G.vhalf[v].size() >= 3)
            essv.push_back(v);
    for (int x : essv) {
        int dx = G.vhalf[x].size();
        // piece: G minus x, edges at x become hairs u -> leaf (u = other endpoint)
        vector<pair<int, int>> ped;
        vector<int> pmap(V, -1);
        int pv = 0;
        for (int v = 0; v < V; v++)
            if (v != x)
                pmap[v] = pv++;
        vector<int> pedge_of(E, -1); // Gamma edge -> piece edge
        for (int e = 0; e < E; e++) {
            int a = G.tail[e], b = G.head[e];
            if (a == x && b == x) {
                fprintf(stderr, "loops at x unsupported\n");
                return 1;
            }
            if (a == x) {
                pedge_of[e] = ped.size();
                ped.push_back({pmap[b], pv++});
            } else if (b == x) {
                pedge_of[e] = ped.size();
                ped.push_back({pmap[a], pv++});
            } else {
                pedge_of[e] = ped.size();
                ped.push_back({pmap[a], pmap[b]});
            }
        }
        Graph P;
        P.build(pv, ped);
        // star at x
        vector<pair<int, int>> sed;
        for (int i = 0; i < dx; i++)
            sed.push_back({0, i + 1});
        Graph S;
        S.build(dx + 1, sed);
        map<int, int> armX;
        for (int i = 0; i < dx; i++)
            armX[G.vhalf[x][i]] = i;
        for (int k = 2; k <= n - 1; k++) {
            int m = n - k;
            vector<Key> SC1, PC1;
            auto zs = z1_basis(S, k, SC1);
            auto zp = z1_basis(P, m, PC1);
            if (mode == 1) {
                // keep only piece cycles supported at a single vertex of the piece? (star classes): filter cycles
                // whose cells all have the same occupied vertex and no particle crossing a middle edge
                vector<vector<pair<int, int>>> keep;
                for (auto &z : zp) {
                    int vtx = -1;
                    bool ok = true;
                    for (auto [i, s0] : z) {
                        Cell c = decode(P, m, PC1[i]);
                        int occ = -1;
                        for (int v = 0; v < P.V; v++)
                            if (c.vocc[v] >= 0)
                                occ = v;
                        if (vtx < 0)
                            vtx = occ;
                        else if (vtx != occ)
                            ok = false;
                    }
                    if (ok)
                        keep.push_back(z);
                }
                zp.swap(keep);
            }
            vector<Cell> cs, cp;
            for (Key q : SC1)
                cs.push_back(decode(S, k, q));
            for (Key q : PC1)
                cp.push_back(decode(P, m, q));
            vector<int> sel(n, 0);
            fill(sel.begin(), sel.begin() + k, 1);
            sort(sel.begin(), sel.end(), greater<int>());
            do {
                vector<int> A, B;
                for (int i = 0; i < n; i++)
                    (sel[i] ? A : B).push_back(i);
                for (auto &z1 : zs)
                    for (auto &z2 : zp) {
                        vector<pair<int, uint32_t>> row;
                        for (auto [i1, s1] : z1)
                            for (auto [i2, s2] : z2) {
                                const Cell &a = cs[i1];
                                const Cell &b = cp[i2];
                                Cell c;
                                c.vocc.assign(V, -1);
                                c.vh.assign(V, -1);
                                c.el.assign(E, {});
                                c.vocc[x] = A[a.vocc[0]];
                                c.vh[x] = G.vhalf[x][a.vh[0] / 2];
                                for (int v = 0; v < V; v++)
                                    if (v != x && b.vocc[pmap[v]] >= 0) {
                                        c.vocc[v] = B[b.vocc[pmap[v]]];
                                        int ph = b.vh[pmap[v]];
                                        int pe = ph / 2;
                                        // piece half-edge -> Gamma half-edge at v
                                        int ge = -1;
                                        for (int e = 0; e < E; e++)
                                            if (pedge_of[e] == pe)
                                                ge = e;
                                        c.vh[v] =
                                            (G.tail[ge] == v && (ph % 2 == 0 || P.tail[pe] == P.head[pe]))
                                                ? 2 * ge
                                                : (G.head[ge] == v ? 2 * ge + 1 : 2 * ge);
                                        if (G.tail[ge] != x && G.head[ge] != x)
                                            c.vh[v] = 2 * ge + (ph % 2);
                                    }
                                for (int e = 0; e < E; e++) {
                                    int pe = pedge_of[e];
                                    vector<int> word;
                                    if (G.tail[e] != x && G.head[e] != x) {
                                        for (int p : b.el[pe])
                                            word.push_back(B[p]);
                                    } else {
                                        // hair: piece word listed from u outward; star arm word listed from x outward
                                        vector<int> xs, us;
                                        for (int p : a.el[armX[G.tail[e] == x ? 2 * e : 2 * e + 1]])
                                            xs.push_back(A[p]);
                                        for (int p : b.el[pe])
                                            us.push_back(B[p]);
                                        if (G.tail[e] == x) {
                                            word = xs;
                                            for (int t = (int)us.size() - 1; t >= 0; t--)
                                                word.push_back(us[t]);
                                        } else {
                                            word = us;
                                            for (int t = (int)xs.size() - 1; t >= 0; t--)
                                                word.push_back(xs[t]);
                                        }
                                    }
                                    c.el[e] = word;
                                }
                                int j = lookup(C2, encode(G, n, c));
                                if (j < 0) {
                                    fprintf(stderr, "BUG: product cell missing\n");
                                    return 2;
                                }
                                int vocc = -1;
                                for (int v = 0; v < V; v++)
                                    if (v != x && c.vocc[v] >= 0)
                                        vocc = v;
                                int kos = vocc < x ? 1 : -1; // Koszul sign: sign of d_v in Gamma
                                row.push_back({j, (uint32_t)(((ll)s1 * s2 * kos + (ll)PR) % PR)});
                            }
                        rows.push_back(row);
                        ngen++;
                    }
            } while (prev_permutation(sel.begin(), sel.end()));
        }
    }
    // sanity: generators are cycles
    {
        vector<Key> C1;
        enumerate(G, n, 1, C1);
        vector<pair<Key, int>> bd;
        mt19937 rng(7);
        for (int t = 0; t < 300 && ngen; t++) {
            auto &g = rows[nb + rng() % ngen];
            map<int, ll> acc;
            for (auto [j, cf] : g) {
                boundary(G, n, decode(G, n, C2[j]), bd);
                ll cc = cf > PR / 2 ? (ll)cf - (ll)PR : cf;
                for (auto &y : bd)
                    acc[lookup(C1, y.first)] += cc * y.second;
            }
            for (auto &kv : acc)
                if (kv.second % (ll)PR) {
                    fprintf(stderr, "BUG: generator not a cycle\n");
                    return 3;
                }
        }
    }
    FpElim Mb(nb, C2.size());
    for (size_t i = 0; i < nb; i++)
        Mb.rows[i] = rows[i];
    Mb.finalize();
    Mb.run();
    long long rB = Mb.piv;
    FpElim M(rows.size(), C2.size());
    for (size_t i = 0; i < rows.size(); i++)
        M.rows[i].swap(rows[i]);
    vector<vector<pair<int, uint32_t>>>().swap(rows);
    M.finalize();
    M.run();
    printf("{\"n\":%d,\"prime\":%llu,\"mode\":%d,\"gens\":%lld,\"rankB2\":%lld,\"rank_span_mod_B\":%lld}\n",
           n, PR, mode, ngen, rB, M.piv - rB);
}
