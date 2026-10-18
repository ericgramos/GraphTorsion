// oconf: integral homology of ordered configuration spaces of graphs
// via the reduced ordered Swiatkowski complex (Lutgehetmann, Chettih-Lutgehetmann,
// Hainaut-Knudsen-Wawrykow Def. 3.12).
//
// Input (stdin):  V E n
//                 then E lines "u v" (edge oriented u -> v; loops u==v allowed,
//                 multi-edges allowed).
// Output (stdout): per degree k: |C_k|; per boundary map d_k: rank, number of unit
// pivots, and the residual (non-unit) block written to residual_k.txt for exact SNF.
//
// Cells: each particle sits either
//   * at a vertex v as the reduced generator (h - h0(v)), h a non-base half-edge at v,
//     at most one particle per vertex; or
//   * in the interior of an edge e, edges carrying linear orders (tail -> head).
// Degree = number of particles at vertices.
// d = sum_i (-1)^i [ p_i enters e(h) at h's end  -  p_i enters e(h0) at h0's end ],
// i indexing the occupied vertices in increasing vertex order.

#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

int V, E, N;
bool checksq = false;
vector<int> pcolor; bool colored = false;
vector<int> gperm; bool grouped = false; vector<vector<int>> gelts; bool transposeTop = false; long long LQ = 0; int LK = 0;
vector<int> etail, ehead;
vector<int> hvert;           // half-edge -> vertex; half-edge 2e = tail end, 2e+1 = head end
vector<vector<int>> vhalf;   // vertex -> half-edges
vector<int> base;            // vertex -> base half-edge (-1 if no half-edges)

// cell encoding: 2N bytes; for particle p: loc (half-edge id if at vertex, 2E+e if on edge), rank
struct Cell {
    vector<int> vocc;            // vertex -> particle or -1
    vector<int> vh;              // vertex -> half-edge used
    vector<vector<int>> elist;   // edge -> particles tail..head
};

typedef unsigned __int128 Key;
Key encode(const Cell& c) {
    Key s = 0;
    for (int v = 0; v < V; v++) if (c.vocc[v] >= 0) {
        int p = c.vocc[v];
        s |= (Key)(unsigned)(c.vh[v]) << (12 * p);
    }
    for (int e = 0; e < E; e++) for (int r = 0; r < (int)c.elist[e].size(); r++) {
        int p = c.elist[e][r];
        s |= (Key)(unsigned)((2 * E + e) | (r << 8)) << (12 * p);
    }
    return s;
}

Key encode(const Cell& c);
Cell decode(Key s);
void canonicalize(Cell& c) {
    if (grouped) {
        Key best = encode(c);
        for (auto& g : gelts) {
            Cell cur = c;
            for (int v = 0; v < V; v++) if (cur.vocc[v] >= 0) cur.vocc[v] = g[cur.vocc[v]];
            for (int e = 0; e < E; e++) for (auto& p : cur.elist[e]) p = g[p];
            Key k = encode(cur); if (k < best) best = k;
        }
        c = decode(best); return;
    }
    if (!colored) return;
    // pools of labels per colour, in increasing order
    map<int, vector<int>> pool;
    for (int p = 0; p < N; p++) pool[pcolor[p]].push_back(p);
    map<int, int> used;
    auto relabel = [&](int p) { int col = pcolor[p]; return pool[col][used[col]++]; };
    for (int v = 0; v < V; v++) if (c.vocc[v] >= 0) c.vocc[v] = relabel(c.vocc[v]);
    for (int e = 0; e < E; e++) for (auto& p : c.elist[e]) p = relabel(p);
}
Cell decode(Key s) {
    Cell c; c.vocc.assign(V, -1); c.vh.assign(V, -1); c.elist.assign(E, {});
    vector<vector<pair<int,int>>> tmp(E);
    for (int p = 0; p < N; p++) {
        unsigned w = (unsigned)((s >> (12 * p)) & 0xFFF);
        int loc = w & 0xFF, r = w >> 8;
        if (loc < 2 * E) { int v = hvert[loc]; c.vocc[v] = p; c.vh[v] = loc; }
        else tmp[loc - 2 * E].push_back({r, p});
    }
    for (int e = 0; e < E; e++) { sort(tmp[e].begin(), tmp[e].end()); for (auto& x : tmp[e]) c.elist[e].push_back(x.second); }
    return c;
}

// enumerate all cells of degree k
void enumerate(int k, vector<Key>& out) {
    Cell c; c.vocc.assign(V, -1); c.vh.assign(V, -1); c.elist.assign(E, {});
    int moving = 0;
    function<void(int)> rec = [&](int p) {
        if (p == N) { if (moving == k) { if (grouped) { Cell d = c; canonicalize(d); Key kc = encode(c); if (encode(d) == kc) out.push_back(kc); } else if (colored) { Cell d = c; canonicalize(d); out.push_back(encode(d)); } else out.push_back(encode(c)); } return; }
        int remaining = N - p;
        // place at a vertex
        if (moving < k) {
            for (int v = 0; v < V; v++) if (c.vocc[v] < 0) {
                for (int h : vhalf[v]) if (h != base[v]) {
                    c.vocc[v] = p; c.vh[v] = h; moving++;
                    rec(p + 1);
                    moving--; c.vocc[v] = -1; c.vh[v] = -1;
                }
            }
        }
        // place on an edge (only if we can still reach k moving)
        if (k - moving <= remaining - 1) {
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

// ---------------- sparse elimination ----------------
// Matrix rows = k-cells, cols = (k-1)-cells (i.e. transpose of d_k).
struct SparseElim {
    int nr, nc;
    vector<vector<pair<int, int>>> rows;
    vector<vector<int>> colrows;
    vector<int> colcnt;
    vector<char> rowdead, coldead;
    ll P; // 0 = integers (unit pivots only); else prime modulus
    long long pivots = 0;
    ll maxabs = 1;

    SparseElim(int nr_, int nc_, ll P_) : nr(nr_), nc(nc_), rows(nr_), colrows(nc_), colcnt(nc_, 0), rowdead(nr_, 0), coldead(nc_, 0), P(P_) {}

    void finalize() {
        for (int r = 0; r < nr; r++) {
            auto& R = rows[r];
            sort(R.begin(), R.end());
            vector<pair<int, ll>> M;
            for (auto& x : R) { if (!M.empty() && M.back().first == x.first) M.back().second += x.second; else M.push_back(x); }
            vector<pair<int, int>> M2;
            for (auto& x : M) { ll v = x.second; if (P) { v %= P; if (v < 0) v += P; } if (llabs(v) > INT_MAX) { fprintf(stderr, "overflow\n"); exit(3); } if (v != 0) M2.push_back({x.first, (int)v}); }
            R.swap(M2);
            for (auto& x : R) { colrows[x.first].push_back(r); colcnt[x.first]++; }
        }
    }
    bool isunit(ll v) const { if (P) return v != 0; return v == 1 || v == -1; }
    ll getval(int r, int c) const {
        auto& R = rows[r];
        auto it = lower_bound(R.begin(), R.end(), make_pair(c, INT_MIN));
        if (it != R.end() && it->first == c) return it->second;
        return 0;
    }
    static ll powmod(ll a, ll e, ll m) { ll r = 1; a %= m; while (e) { if (e & 1) r = (__int128)r * a % m; a = (__int128)a * a % m; e >>= 1; } return r; }

    void pivot(int r, int c) {
        ll u = getval(r, c);
        ll uinv = P ? powmod(u, P - 2, P) : u; // u = +-1 over Z
        vector<pair<int, int>>& R = rows[r];
        // collect rows to update
        vector<int> targets;
        for (int r2 : colrows[c]) if (r2 != r && !rowdead[r2]) targets.push_back(r2);
        sort(targets.begin(), targets.end()); targets.erase(unique(targets.begin(), targets.end()), targets.end());
        for (int r2 : targets) {
            ll a = getval(r2, c);
            if (a == 0) continue;
            ll f;
            if (P) f = (__int128)a * uinv % P; else f = a * uinv;
            auto& S = rows[r2];
            vector<pair<int, int>> out; out.reserve(S.size() + R.size());
            size_t i = 0, j = 0;
            while (i < S.size() || j < R.size()) {
                if (j == R.size() || (i < S.size() && S[i].first < R[j].first)) { out.push_back(S[i]); i++; }
                else if (i == S.size() || R[j].first < S[i].first) {
                    ll v;
                    if (P) { v = (P - (__int128)f * R[j].second % P) % P; }
                    else { if (__builtin_mul_overflow(f, R[j].second, &v)) { fprintf(stderr, "overflow\n"); exit(3); } v = -v; }
                    if (llabs(v) > INT_MAX) { fprintf(stderr, "overflow\n"); exit(3); }
                    if (v != 0) { out.push_back({R[j].first, (int)v}); colrows[R[j].first].push_back(r2); colcnt[R[j].first]++; }
                    j++;
                } else {
                    ll v;
                    if (P) { v = ((S[i].second - (__int128)f * R[j].second % P) % P + P) % P; }
                    else { ll t; if (__builtin_mul_overflow(f, R[j].second, &t) || __builtin_sub_overflow(S[i].second, t, &v)) { fprintf(stderr, "overflow\n"); exit(3); } }
                    if (llabs(v) > INT_MAX) { fprintf(stderr, "overflow\n"); exit(3); }
                    if (v != 0) { out.push_back({S[i].first, (int)v}); if (!P) maxabs = max(maxabs, llabs(v)); }
                    else colcnt[S[i].first]--;
                    i++; j++;
                }
            }
            S.swap(out);
        }
        for (auto& x : R) colcnt[x.first]--;
        rowdead[r] = 1; coldead[c] = 1; colcnt[c] = 0;
        vector<pair<int, int>>().swap(rows[r]);
        vector<int>().swap(colrows[c]);
        pivots++;
    }

    void compact_col(int c) {
        auto& L = colrows[c];
        vector<int> M;
        for (int r : L) if (!rowdead[r] && getval(r, c) != 0) M.push_back(r);
        sort(M.begin(), M.end()); M.erase(unique(M.begin(), M.end()), M.end());
        L.swap(M);
    }

    // returns true if some pivot was done
    bool sweep() {
        // priority by column count, Markowitz within column
        typedef pair<int, int> PI;
        priority_queue<PI, vector<PI>, greater<PI>> pq;
        for (int c = 0; c < nc; c++) if (!coldead[c] && colcnt[c] > 0) pq.push({colcnt[c], c});
        bool any = false;
        while (!pq.empty()) {
            auto [cnt, c] = pq.top(); pq.pop();
            if (coldead[c] || colcnt[c] != cnt || cnt == 0) continue;
            if (colrows[c].size() > 2 * (size_t)colcnt[c] + 8) compact_col(c);
            int best = -1; size_t bestlen = SIZE_MAX;
            for (int r : colrows[c]) {
                if (rowdead[r]) continue;
                ll v = getval(r, c);
                if (v != 0 && isunit(v) && rows[r].size() < bestlen) { bestlen = rows[r].size(); best = r; }
            }
            if (best < 0) continue; // no unit pivot in this column right now
            // record neighbours whose counts change
            vector<int> touched;
            for (auto& x : rows[best]) touched.push_back(x.first);
            pivot(best, c);
            any = true;
            for (int c2 : touched) if (!coldead[c2] && colcnt[c2] > 0) pq.push({colcnt[c2], c2});
            // fill may also have changed other columns; they are re-pushed lazily below
            // (columns touched through fill are a superset of 'touched' only via R entries, which we pushed)
        }
        return any;
    }

    void run() { while (sweep()) {} }

    void dump_residual(const string& fn, int& rr, int& rc, long long& rnnz) {
        vector<int> rmap(nr, -1), cmap(nc, -1);
        rr = 0; rc = 0; rnnz = 0;
        for (int r = 0; r < nr; r++) if (!rowdead[r] && !rows[r].empty()) rmap[r] = rr++;
        for (int c = 0; c < nc; c++) if (!coldead[c] && colcnt[c] > 0) cmap[c] = rc++;
        FILE* f = fopen(fn.c_str(), "w");
        fprintf(f, "%d %d\n", rr, rc);
        for (int r = 0; r < nr; r++) if (rmap[r] >= 0) for (auto& x : rows[r]) { fprintf(f, "%d %d %lld\n", rmap[r], cmap[x.first], (long long)x.second); rnnz++; }
        fclose(f);
    }
};

struct F2E {
    int nr, nc;
    vector<vector<int>> rows;
    vector<vector<int>> colrows;
    vector<int> colcnt;
    vector<char> rowdead, coldead;
    vector<pair<int, vector<int>>> piv;   // (pivot column, pivot row at pivot time)
    bool keep; long long pivots = 0; ll maxabs = 1;
    F2E(int r, int c, bool keep_) : nr(r), nc(c), rows(r), colrows(c), colcnt(c, 0), rowdead(r, 0), coldead(c, 0), keep(keep_) {}
    void finalize() {
        for (int r = 0; r < nr; r++) {
            auto& R = rows[r];
            sort(R.begin(), R.end());
            vector<int> M;
            for (size_t i = 0; i < R.size();) { size_t j = i; while (j < R.size() && R[j] == R[i]) j++; if ((j - i) & 1) M.push_back(R[i]); i = j; }
            R.swap(M);
            for (int c : R) { colrows[c].push_back(r); colcnt[c]++; }
        }
    }
    bool has(int r, int c) const { auto& R = rows[r]; return binary_search(R.begin(), R.end(), c); }
    void pivot(int r, int c) {
        vector<int>& R = rows[r];
        vector<int> targets;
        for (int r2 : colrows[c]) if (r2 != r && !rowdead[r2]) targets.push_back(r2);
        sort(targets.begin(), targets.end()); targets.erase(unique(targets.begin(), targets.end()), targets.end());
        for (int r2 : targets) {
            if (!has(r2, c)) continue;
            auto& S = rows[r2];
            vector<int> out; out.reserve(S.size() + R.size());
            size_t i = 0, j = 0;
            while (i < S.size() || j < R.size()) {
                if (j == R.size() || (i < S.size() && S[i] < R[j])) { out.push_back(S[i]); i++; }
                else if (i == S.size() || R[j] < S[i]) { out.push_back(R[j]); colrows[R[j]].push_back(r2); colcnt[R[j]]++; j++; }
                else { colcnt[S[i]]--; i++; j++; }
            }
            S.swap(out);
        }
        for (int x : R) colcnt[x]--;
        rowdead[r] = 1; coldead[c] = 1; colcnt[c] = 0;
        pivots++;
        if (keep) piv.push_back({c, std::move(rows[r])});
        vector<int>().swap(rows[r]);
        vector<int>().swap(colrows[c]);
    }
    void run() {
        bool any = true;
        while (any) {
            any = false;
            typedef pair<int, int> PI;
            priority_queue<PI, vector<PI>, greater<PI>> pq;
            for (int c = 0; c < nc; c++) if (!coldead[c] && colcnt[c] > 0) pq.push({colcnt[c], c});
            while (!pq.empty()) {
                auto [cnt, c] = pq.top(); pq.pop();
                if (coldead[c] || colcnt[c] != cnt || cnt == 0) continue;
                if (colrows[c].size() > 2 * (size_t)colcnt[c] + 8) {
                    vector<int> M; for (int r : colrows[c]) if (!rowdead[r] && has(r, c)) M.push_back(r);
                    sort(M.begin(), M.end()); M.erase(unique(M.begin(), M.end()), M.end()); colrows[c].swap(M);
                }
                int best = -1; size_t bl = SIZE_MAX;
                for (int r : colrows[c]) if (!rowdead[r] && rows[r].size() < bl && has(r, c)) { bl = rows[r].size(); best = r; }
                if (best < 0) continue;
                vector<int> touched(rows[best].begin(), rows[best].end());
                pivot(best, c); any = true;
                for (int c2 : touched) if (!coldead[c2] && colcnt[c2] > 0) pq.push({colcnt[c2], c2});
            }
        }
    }
};

struct SPE {
    int nr, nc; uint32_t p;
    vector<vector<uint32_t>> rows; vector<vector<int>> colrows; vector<int> colcnt; vector<char> rowdead, coldead;
    long long pivots = 0; ll maxabs = 1;
    SPE(int r, int c, uint32_t p_) : nr(r), nc(c), p(p_), rows(r), colrows(c), colcnt(c, 0), rowdead(r, 0), coldead(c, 0) {}
    static uint32_t col(uint32_t x) { return x >> 4; } static uint32_t val(uint32_t x) { return x & 15; }
    uint32_t inv(uint32_t a) const { for (uint32_t b = 1; b < p; b++) if (a * b % p == 1) return b; return 0; }
    void finalize() {
        for (int r = 0; r < nr; r++) {
            auto& R = rows[r]; sort(R.begin(), R.end());
            vector<uint32_t> M;
            for (size_t i = 0; i < R.size();) { uint32_t c = col(R[i]); uint32_t v = 0; while (i < R.size() && col(R[i]) == c) { v = (v + val(R[i])) % p; i++; } if (v) M.push_back(c << 4 | v); }
            R.swap(M); for (auto x : R) { colrows[col(x)].push_back(r); colcnt[col(x)]++; }
        }
    }
    uint32_t get(int r, int c) const { auto& R = rows[r]; auto it = lower_bound(R.begin(), R.end(), (uint32_t)c << 4); return (it != R.end() && col(*it) == (uint32_t)c) ? val(*it) : 0; }
    void pivot(int r, int c) {
        auto& R = rows[r]; uint32_t ui = inv(get(r, c));
        vector<int> tg; for (int r2 : colrows[c]) if (r2 != r && !rowdead[r2]) tg.push_back(r2);
        sort(tg.begin(), tg.end()); tg.erase(unique(tg.begin(), tg.end()), tg.end());
        for (int r2 : tg) {
            uint32_t a = get(r2, c); if (!a) continue; uint32_t f = a * ui % p;
            auto& S = rows[r2]; vector<uint32_t> out; out.reserve(S.size() + R.size()); size_t i = 0, j = 0;
            while (i < S.size() || j < R.size()) {
                if (j == R.size() || (i < S.size() && col(S[i]) < col(R[j]))) { out.push_back(S[i]); i++; }
                else if (i == S.size() || col(R[j]) < col(S[i])) { uint32_t v = (p - f * val(R[j]) % p) % p; if (v) { out.push_back(col(R[j]) << 4 | v); colrows[col(R[j])].push_back(r2); colcnt[col(R[j])]++; } j++; }
                else { uint32_t v = (val(S[i]) + p - f * val(R[j]) % p) % p; if (v) out.push_back(col(S[i]) << 4 | v); else colcnt[col(S[i])]--; i++; j++; }
            }
            S.swap(out);
        }
        for (auto x : R) colcnt[col(x)]--;
        rowdead[r] = 1; coldead[c] = 1; colcnt[c] = 0; pivots++;
        vector<uint32_t>().swap(rows[r]); vector<int>().swap(colrows[c]);
    }
    void run() {
        bool any = true;
        while (any) { any = false;
            typedef pair<int,int> PI; priority_queue<PI, vector<PI>, greater<PI>> pq;
            for (int c = 0; c < nc; c++) if (!coldead[c] && colcnt[c] > 0) pq.push({colcnt[c], c});
            while (!pq.empty()) { auto [cnt, c] = pq.top(); pq.pop();
                if (coldead[c] || colcnt[c] != cnt || cnt == 0) continue;
                if (colrows[c].size() > 2 * (size_t)colcnt[c] + 8) { vector<int> M; for (int r : colrows[c]) if (!rowdead[r] && get(r, c)) M.push_back(r); sort(M.begin(), M.end()); M.erase(unique(M.begin(), M.end()), M.end()); colrows[c].swap(M); }
                int best = -1; size_t bl = SIZE_MAX;
                for (int r : colrows[c]) if (!rowdead[r] && rows[r].size() < bl && get(r, c)) { bl = rows[r].size(); best = r; }
                if (best < 0) continue;
                vector<int> touched; for (auto x : rows[best]) touched.push_back(col(x));
                pivot(best, c); any = true;
                for (int c2 : touched) if (!coldead[c2] && colcnt[c2] > 0) pq.push({colcnt[c2], c2});
            } }
    }
};

// packed local elimination over Z/p^k: entry = col<<5 | val, val < 32. Pivots only on units mod p;
// after each stage the remaining block is divided by p (all live entries are then divisible by p).
struct LPE {
    int nr, nc; uint32_t p0, mod;
    vector<vector<uint32_t>> rows; vector<vector<int>> colrows; vector<int> colcnt; vector<char> rowdead, coldead;
    long long pivots = 0; ll maxabs = 1; vector<long long> stages; long long leftover_nonzero = 0;
    LPE(int r, int c, uint32_t p, uint32_t m) : nr(r), nc(c), p0(p), mod(m), rows(r), colrows(c), colcnt(c, 0), rowdead(r, 0), coldead(c, 0) {}
    static uint32_t col(uint32_t x) { return x >> 5; } static uint32_t val(uint32_t x) { return x & 31; }
    uint32_t inv(uint32_t a) const { for (uint32_t b = 1; b < mod; b++) if (a * b % mod == 1) return b; return 0; }
    void finalize() {
        for (int r = 0; r < nr; r++) {
            auto& R = rows[r]; sort(R.begin(), R.end());
            vector<uint32_t> M;
            for (size_t i = 0; i < R.size();) { uint32_t c = col(R[i]); uint32_t v = 0; while (i < R.size() && col(R[i]) == c) { v = (v + val(R[i])) % mod; i++; } if (v) M.push_back(c << 5 | v); }
            R.swap(M); for (auto x : R) { colrows[col(x)].push_back(r); colcnt[col(x)]++; }
        }
    }
    uint32_t get(int r, int c) const { auto& R = rows[r]; auto it = lower_bound(R.begin(), R.end(), (uint32_t)c << 5); return (it != R.end() && col(*it) == (uint32_t)c) ? val(*it) : 0; }
    bool unit(uint32_t v) const { return v % p0 != 0; }
    void pivot(int r, int c) {
        auto& R = rows[r]; uint32_t ui = inv(get(r, c));
        vector<int> tg; for (int r2 : colrows[c]) if (r2 != r && !rowdead[r2]) tg.push_back(r2);
        sort(tg.begin(), tg.end()); tg.erase(unique(tg.begin(), tg.end()), tg.end());
        for (int r2 : tg) {
            uint32_t a = get(r2, c); if (!a) continue; uint32_t f = a * ui % mod;
            auto& S = rows[r2]; vector<uint32_t> out; out.reserve(S.size() + R.size()); size_t i = 0, j = 0;
            while (i < S.size() || j < R.size()) {
                if (j == R.size() || (i < S.size() && col(S[i]) < col(R[j]))) { out.push_back(S[i]); i++; }
                else if (i == S.size() || col(R[j]) < col(S[i])) { uint32_t v = (mod - f * val(R[j]) % mod) % mod; if (v) { out.push_back(col(R[j]) << 5 | v); colrows[col(R[j])].push_back(r2); colcnt[col(R[j])]++; } j++; }
                else { uint32_t v = (val(S[i]) + mod - f * val(R[j]) % mod) % mod; if (v) out.push_back(col(S[i]) << 5 | v); else colcnt[col(S[i])]--; i++; j++; }
            }
            S.swap(out);
        }
        for (auto x : R) colcnt[col(x)]--;
        rowdead[r] = 1; coldead[c] = 1; colcnt[c] = 0; pivots++;
        vector<uint32_t>().swap(rows[r]); vector<int>().swap(colrows[c]);
    }
    bool sweep() {
        bool any = false;
        typedef pair<int,int> PI; priority_queue<PI, vector<PI>, greater<PI>> pq;
        for (int c = 0; c < nc; c++) if (!coldead[c] && colcnt[c] > 0) pq.push({colcnt[c], c});
        while (!pq.empty()) { auto [cnt, c] = pq.top(); pq.pop();
            if (coldead[c] || colcnt[c] != cnt || cnt == 0) continue;
            if (colrows[c].size() > 2 * (size_t)colcnt[c] + 8) { vector<int> M; for (int r : colrows[c]) if (!rowdead[r] && get(r, c)) M.push_back(r); sort(M.begin(), M.end()); M.erase(unique(M.begin(), M.end()), M.end()); colrows[c].swap(M); }
            int best = -1; size_t bl = SIZE_MAX;
            for (int r : colrows[c]) if (!rowdead[r] && rows[r].size() < bl) { uint32_t v = get(r, c); if (v && unit(v)) { bl = rows[r].size(); best = r; } }
            if (best < 0) continue;
            vector<int> touched; for (auto x : rows[best]) touched.push_back(col(x));
            pivot(best, c); any = true;
            for (int c2 : touched) if (!coldead[c2] && colcnt[c2] > 0) pq.push({colcnt[c2], c2});
        }
        return any;
    }
    void run() {
        long long before = 0;
        while (true) {
            while (sweep()) {}
            stages.push_back(pivots - before); before = pivots;
            bool live = false; for (int r = 0; r < nr && !live; r++) if (!rowdead[r] && !rows[r].empty()) live = true;
            if (!live) break;
            if (mod == p0) { for (int r = 0; r < nr; r++) if (!rowdead[r]) leftover_nonzero += rows[r].size(); break; }
            // all live entries divisible by p0: divide out
            for (int r = 0; r < nr; r++) if (!rowdead[r]) {
                vector<uint32_t> M; for (auto x : rows[r]) { if (val(x) % p0) { fprintf(stderr, "local: nonunit not divisible\n"); exit(4); } uint32_t v = val(x) / p0; if (v) M.push_back(col(x) << 5 | v); else colcnt[col(x)]--; }
                rows[r].swap(M);
            }
            mod /= p0;
        }
    }
};
void f2_dump(const string& fn, int& rr, int& rc, long long& rnnz) { FILE* f = fopen(fn.c_str(), "w"); fprintf(f, "0 0\n"); fclose(f); rr = rc = 0; rnnz = 0; }

int main(int argc, char** argv) {
    ll P = 0; string prefix = "res"; int maxdeg = -1; bool cellsonly = false;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (a == "-p") P = atoll(argv[++i]);
        else if (a == "-q") { LQ = atoll(argv[++i]); LK = atoi(argv[++i]); P = 1; for (int t = 0; t < LK; t++) P *= LQ; }
        else if (a == "-o") prefix = argv[++i];
        else if (a == "-d") maxdeg = atoi(argv[++i]);
        else if (a == "-c") cellsonly = true;
        else if (a == "-z") checksq = true;
        else if (a == "-T") transposeTop = true;
        else if (a == "-g") { grouped = true; string all = argv[++i]; stringstream gs(all); string one;
            while (getline(gs, one, ';')) { vector<int> pm; stringstream ss(one); string t; while (getline(ss, t, ',')) pm.push_back(stoi(t)); gelts.push_back(pm); } }
        else if (a == "-k") { colored = true; string cs = argv[++i]; for (char ch : cs) pcolor.push_back(ch - 'a'); }
    }
    for (int i = 1; i < argc; i++) if (string(argv[i]) == "-m") {
        // generic mode: eliminate a sparse matrix "rows cols\n r c v ..." from stdin
        int nr, nc; if (scanf("%d %d", &nr, &nc) != 2) return 1;
        SparseElim M(nr, nc, P);
        int r, c; ll v;
        while (scanf("%d %d %lld", &r, &c, &v) == 3) M.rows[r].push_back({c, v});
        M.finalize(); M.run();
        int rr, rc; long long rnnz;
        M.dump_residual(prefix + "_m.txt", rr, rc, rnnz);
        printf("{\"pivots\":%lld,\"res_rows\":%d,\"res_cols\":%d,\"res_file\":\"%s\"}\n", M.pivots, rr, rc, (prefix + "_m.txt").c_str());
        return 0;
    }
    if (scanf("%d %d %d", &V, &E, &N) != 3) return 1;
    if (grouped) {   // extend generators to N labels and close under composition
        for (auto& g : gelts) for (int t = g.size(); t < N; t++) g.push_back(t);
        vector<int> id(N); iota(id.begin(), id.end(), 0);
        set<vector<int>> G = {id}; vector<vector<int>> frontier = {id};
        vector<vector<int>> gens = gelts;
        while (!frontier.empty()) { vector<vector<int>> nf;
            for (auto& x : frontier) for (auto& g : gens) { vector<int> y(N); for (int t = 0; t < N; t++) y[t] = g[x[t]]; if (G.insert(y).second) nf.push_back(y); }
            frontier.swap(nf); }
        gelts.clear(); for (auto& x : G) if (x != id) gelts.push_back(x);
        fprintf(stderr, "group order %zu\n", gelts.size() + 1);
    }
    etail.resize(E); ehead.resize(E);
    hvert.resize(2 * E); vhalf.assign(V, {});
    for (int e = 0; e < E; e++) {
        scanf("%d %d", &etail[e], &ehead[e]);
        hvert[2 * e] = etail[e]; hvert[2 * e + 1] = ehead[e];
        vhalf[etail[e]].push_back(2 * e); vhalf[ehead[e]].push_back(2 * e + 1);
    }
    if (3 * E >= 255 || N > 10 || N > 15) { fprintf(stderr, "too large\n"); return 1; }
    base.assign(V, -1);
    for (int v = 0; v < V; v++) if (!vhalf[v].empty()) base[v] = vhalf[v][0];
    int D = 0;
    for (int v = 0; v < V; v++) if (vhalf[v].size() >= 2) D++;
    D = min(D, N);
    if (maxdeg >= 0) D = min(D, maxdeg);

    if (checksq) {
        vector<vector<Key>> C(D + 1);
        for (int k = 0; k <= D; k++) { enumerate(k, C[k]); sort(C[k].begin(), C[k].end()); if (colored) C[k].erase(unique(C[k].begin(), C[k].end()), C[k].end()); }
        auto bd = [&](int k, size_t i, vector<pair<long long,int>>& out) {
            out.clear();
            Cell c = decode(C[k][i]);
            int idx = 0;
            for (int v = 0; v < V; v++) if (c.vocc[v] >= 0) {
                int p = c.vocc[v], h = c.vh[v], h0 = base[v];
                int sgn = (idx % 2 == 0) ? 1 : -1; idx++;
                for (int t = 0; t < 2; t++) {
                    int hh = t == 0 ? h : h0;
                    Cell f = c; f.vocc[v] = -1; f.vh[v] = -1;
                    int e = hh / 2;
                    if (hh % 2 == 0) f.elist[e].insert(f.elist[e].begin(), p); else f.elist[e].push_back(p);
                    canonicalize(f);
                    Key fk = encode(f);
                    auto it = lower_bound(C[k-1].begin(), C[k-1].end(), fk);
                    if (it == C[k-1].end() || *it != fk) { fprintf(stderr, "face not found\n"); exit(2); }
                    out.push_back({(long long)(it - C[k-1].begin()), t == 0 ? sgn : -sgn});
                }
            }
        };
        long long bad = 0, checked = 0;
        vector<pair<long long,int>> b1, b2;
        for (int k = 2; k <= D; k++) {
            for (size_t i = 0; i < C[k].size(); i++) {
                bd(k, i, b1);
                map<long long, long long> acc;
                for (auto& x : b1) { bd(k - 1, x.first, b2); for (auto& y : b2) acc[y.first] += (long long)x.second * y.second; }
                for (auto& z : acc) if (z.second != 0) { bad++; break; }
                checked++;
            }
        }
        printf("{\"d2_checked\":%lld,\"d2_nonzero\":%lld}\n", checked, bad);
        return 0;
    }
    vector<long long> csize(D + 2, 0);
    vector<Key> prevcells;
    vector<char> dropcol;        // (k-1)-cells that were pivot rows of d_{k-1}: dropped from d_k
    bool canDrop = true;
    printf("{\"V\":%d,\"E\":%d,\"n\":%d,\"P\":%lld,\"dim\":%d", V, E, N, P, D);
    vector<string> results;
    for (int k = 0; k <= D; k++) {
        vector<Key> cells;
        enumerate(k, cells);
        sort(cells.begin(), cells.end());
        if (colored || grouped) cells.erase(unique(cells.begin(), cells.end()), cells.end());
        csize[k] = cells.size();
        if (k >= 1 && !cellsonly) {
            auto addent = [&](auto& MM, size_t i, int ci, ll coef) {
                if constexpr (std::is_same_v<std::decay_t<decltype(MM)>, F2E>) { if (coef & 1) MM.rows[i].push_back(ci); }
                else if constexpr (std::is_same_v<std::decay_t<decltype(MM)>, LPE>) { ll v = ((coef % P) + P) % P; if (v) MM.rows[i].push_back((uint32_t)ci << 5 | (uint32_t)v); }
                else if constexpr (std::is_same_v<std::decay_t<decltype(MM)>, SPE>) { ll v = ((coef % P) + P) % P; if (v) MM.rows[i].push_back((uint32_t)ci << 4 | (uint32_t)v); }
                else MM.rows[i].push_back({ci, coef}); };
            auto dumpres = [&](auto& MM, const string& fn, int& rr, int& rc, long long& rnnz) {
                if constexpr (std::is_same_v<std::decay_t<decltype(MM)>, F2E> || std::is_same_v<std::decay_t<decltype(MM)>, SPE> || std::is_same_v<std::decay_t<decltype(MM)>, LPE>) f2_dump(fn, rr, rc, rnnz);
                else MM.dump_residual(fn, rr, rc, rnnz); };
            auto blockfn = [&](auto& M) {
            for (size_t i = 0; i < cells.size(); i++) {
                Cell c = decode(cells[i]);
                int idx = 0;
                for (int v = 0; v < V; v++) if (c.vocc[v] >= 0) {
                    int p = c.vocc[v], h = c.vh[v], h0 = base[v];
                    ll sgn = (idx % 2 == 0) ? 1 : -1;
                    idx++;
                    for (int t = 0; t < 2; t++) {
                        int hh = t == 0 ? h : h0;
                        Cell f = c; f.vocc[v] = -1; f.vh[v] = -1;
                        int e = hh / 2;
                        if (hh % 2 == 0) f.elist[e].insert(f.elist[e].begin(), p);
                        else f.elist[e].push_back(p);
                        canonicalize(f);
                        Key fk = encode(f);
                        auto it = lower_bound(prevcells.begin(), prevcells.end(), fk);
                        if (it == prevcells.end() || *it != fk) { fprintf(stderr, "face not found\n"); exit(2); }
                        ll coef = t == 0 ? sgn : -sgn;
                        if (P) coef = ((coef % P) + P) % P;
                        int ci = (int)(it - prevcells.begin());
                        if (!(canDrop && !dropcol.empty() && dropcol[ci])) addent(M, i, ci, coef);
                    }
                }
            }
            vector<Key>().swap(prevcells);   // no longer needed; frees memory before elimination
            if constexpr (std::is_same_v<std::decay_t<decltype(M)>, F2E>) {
                if (transposeTop && k == D) {
                    // rank-only: eliminate the transpose (rows = (k-1)-cells)
                    F2E Tm(M.nc, M.nr, false);
                    for (int i = 0; i < M.nr; i++) { for (int c : M.rows[i]) Tm.rows[c].push_back(i); vector<int>().swap(M.rows[i]); }
                    Tm.finalize(); Tm.run();
                    M.pivots = Tm.pivots;
                    char b2[256]; snprintf(b2, sizeof b2, "{\"k\":%d,\"rows\":%zu,\"cols\":%d,\"pivots\":%lld,\"res_rows\":0,\"res_cols\":0,\"res_nnz\":0,\"maxabs\":1,\"res_file\":\"none\"}", k, cells.size(), M.nc, M.pivots);
                    results.push_back(b2); canDrop = false; return;
                }
            }
            M.finalize();
            fprintf(stderr, "deg %d: matrix built, eliminating\n", k);
            M.run();
            int rr, rc; long long rnnz;
            string fn = prefix + "_" + to_string(k) + ".txt";
            dumpres(M, fn, rr, rc, rnnz);
            char buf[512];
            snprintf(buf, sizeof buf, "{\"k\":%d,\"rows\":%zu,\"cols\":%zu,\"pivots\":%lld,\"res_rows\":%d,\"res_cols\":%d,\"res_nnz\":%lld,\"maxabs\":%lld,\"res_file\":\"%s\"}",
                     k, cells.size(), prevcells.size(), M.pivots, rr, rc, rnnz, M.maxabs, fn.c_str());
            string rec = buf;
            if constexpr (std::is_same_v<std::decay_t<decltype(M)>, LPE>) {
                rec.pop_back(); rec += ",\"stages\":[";
                for (size_t t = 0; t < M.stages.size(); t++) rec += (t ? "," : "") + to_string(M.stages[t]);
                rec += "],\"leftover_nnz\":" + to_string(M.leftover_nonzero) + "}";
                canDrop = (M.stages.size() <= 1 && M.leftover_nonzero == 0);
            } else canDrop = (rr == 0);
            results.push_back(rec);
            dropcol.assign(cells.size(), 0);
            for (size_t i = 0; i < cells.size(); i++) if (M.rowdead[i]) dropcol[i] = 1;
            };
            if (LQ) { LPE M(cells.size(), prevcells.size(), (uint32_t)LQ, (uint32_t)P); blockfn(M); }
            else if (P == 2) { F2E M(cells.size(), prevcells.size(), false); blockfn(M); }
            else if (P > 2 && P <= 13) { SPE M(cells.size(), prevcells.size(), (uint32_t)P); blockfn(M); }
            else { SparseElim M(cells.size(), prevcells.size(), P); blockfn(M); }
        }
        prevcells.swap(cells);
        fprintf(stderr, "deg %d: %lld cells\n", k, csize[k]);
    }
    printf(",\"cells\":[");
    for (int k = 0; k <= D; k++) printf("%s%lld", k ? "," : "", csize[k]);
    printf("],\"maps\":[");
    for (size_t i = 0; i < results.size(); i++) printf("%s%s", i ? "," : "", results[i].c_str());
    printf("]}\n");
}
