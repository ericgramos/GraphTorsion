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
struct GE { vector<int> lab, vm, em; vector<char> fl; int chi = 1; int hm[256]; };
int twistLab = 0, twistEdge = 0; vector<int> agSigns; bool fastEnum = false;
int permsign(const vector<int>& p) { int n = p.size(), s = 1; vector<char> seen(n, 0); for (int i = 0; i < n; i++) if (!seen[i]) { int len = 0, j = i; while (!seen[j]) { seen[j] = 1; j = p[j]; len++; } if (len % 2 == 0) s = -s; } return s; }
vector<GE> gall;             // all non-identity group elements (label perm x graph automorphism)
vector<string> agens;        // automorphism generators as given

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
static inline Key applyKey(const GE& g, Key s, int& sign);
int applyGE(const GE& g, const Cell& c, Cell& cur) {
    cur.vocc.assign(V, -1); cur.vh.assign(V, -1); cur.elist.assign(E, {});
    int w[16], k = 0;
    for (int v = 0; v < V; v++) if (c.vocc[v] >= 0) { int a = g.vm[v]; cur.vocc[a] = g.lab[c.vocc[v]]; int h = c.vh[v]; cur.vh[a] = 2 * g.em[h / 2] + ((h & 1) ^ g.fl[h / 2]); w[k++] = a; }
    for (int e = 0; e < E; e++) { auto& L = cur.elist[g.em[e]]; L = c.elist[e]; for (auto& p : L) p = g.lab[p]; if (g.fl[e]) reverse(L.begin(), L.end()); }
    int inv = 0; for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) if (w[i] > w[j]) inv++;
    return ((inv & 1) ? -1 : 1) * g.chi;
}
// returns s with [c] = s*[canonical rep]; 0 if c is zero in the coinvariants (stabilised with sign -1)
vector<GE> Aelts, Lelts;   // factorised group: automorphisms (with chi_A) and label perms (with chi_L), identities included
int canonFactored(Key kc, Key& best) {
    best = kc; int bs = 1;
    for (auto& a : Aelts) {
        int sa; Key ka = applyKey(a, kc, sa);
        unsigned f[16]; for (int p = 0; p < N; p++) f[p] = (unsigned)((ka >> (12 * p)) & 0xFFF);
        for (auto& l : Lelts) {
            Key k = 0; for (int p = 0; p < N; p++) k |= (Key)f[p] << (12 * l.lab[p]);
            int sg = sa * l.chi;
            if (k == kc && sg < 0) return 0;
            if (k < best) { best = k; bs = sg; }
        }
    }
    return bs;
}
int canonicalize(Cell& c) {
    if (grouped) {
        Key kc = encode(c), best; int bs = canonFactored(kc, best);
        if (!bs) return 0;
        c = decode(best); return bs;
    }
    if (!colored) return 1;
    // pools of labels per colour, in increasing order
    map<int, vector<int>> pool;
    for (int p = 0; p < N; p++) pool[pcolor[p]].push_back(p);
    map<int, int> used;
    auto relabel = [&](int p) { int col = pcolor[p]; return pool[col][used[col]++]; };
    for (int v = 0; v < V; v++) if (c.vocc[v] >= 0) c.vocc[v] = relabel(c.vocc[v]);
    for (int e = 0; e < E; e++) for (auto& p : c.elist[e]) p = relabel(p);
    return 1;
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

// group action directly on keys (no allocation)
static inline Key applyKey(const GE& g, Key s, int& sign) {
    int len[128]; unsigned w[16];
    for (int e = 0; e < E; e++) len[e] = 0;
    for (int p = 0; p < N; p++) { w[p] = (unsigned)((s >> (12 * p)) & 0xFFF); unsigned loc = w[p] & 0xFF; if ((int)loc >= 2 * E) len[loc - 2 * E]++; }
    Key out = 0; int occ[16], k = 0;
    for (int p = 0; p < N; p++) {
        unsigned loc = w[p] & 0xFF, r = w[p] >> 8, nl, nr = 0;
        if ((int)loc < 2 * E) { nl = g.hm[loc]; occ[k++] = hvert[loc]; }
        else { int e = loc - 2 * E; nl = 2 * E + g.em[e]; nr = g.fl[e] ? len[e] - 1 - r : r; }
        out |= (Key)(nl | (nr << 8)) << (12 * g.lab[p]);
    }
    // occupied vertices in increasing order -> images; parity of inversions
    sort(occ, occ + k); int inv = 0;
    for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) if (g.vm[occ[i]] > g.vm[occ[j]]) inv++;
    sign = ((inv & 1) ? -1 : 1) * g.chi;
    return out;
}
bool isCanonicalKey(Key kc) {
    for (auto& g : gall) { int sg; Key k = applyKey(g, kc, sg); if (k < kc) return false; if (k == kc && sg < 0) return false; }
    return true;
}
bool isCanonical(const Cell& c) {
    Key kc = encode(c); static thread_local Cell cur;
    for (auto& g : gall) { int sg = applyGE(g, c, cur); Key k = encode(cur); if (k < kc) return false; if (k == kc && sg < 0) return false; }
    return true;
}

// fast orbit enumeration: A-canonical shapes x L-canonical labelings, then full canonicalization
vector<const GE*> gLab, gAut;   // label-only and automorphism-only elements (non-identity)
int canonKey(Key kc, Key& best) {
    best = kc; int bs = 1;
    for (auto& g : gall) { int sg; Key k = applyKey(g, kc, sg); if (k == kc && sg < 0) return 0; if (k < best) { best = k; bs = sg; } }
    return bs;
}
void enumerate2(int k, vector<Key>& out) {
    vector<int> occh(V, -1); vector<int> lens(E, 0);
    vector<unsigned> fields; fields.reserve(N);
    long long shapes = 0, canonShapes = 0;
    function<void(int,int)> recV;   // vertex choices
    function<void(int,int)> recE;   // edge lengths
    auto process = [&]() {
        fields.clear();
        for (int v = 0; v < V; v++) if (occh[v] >= 0) fields.push_back((unsigned)occh[v]);
        for (int e = 0; e < E; e++) for (int r = 0; r < lens[e]; r++) fields.push_back((unsigned)((2 * E + e) | (r << 8)));
        shapes++;
        // A-canonical shape: compare sorted field multisets
        vector<unsigned> sig(fields); sort(sig.begin(), sig.end());
        Key tk = 0; for (int i = 0; i < N; i++) tk |= (Key)fields[i] << (12 * i);
        for (auto* g : gAut) { int sg; Key k = applyKey(*g, tk, sg); vector<unsigned> s2(N); for (int i = 0; i < N; i++) s2[i] = (unsigned)((k >> (12 * i)) & 0xFFF); sort(s2.begin(), s2.end()); if (s2 < sig) return; }
        canonShapes++;
        vector<int> perm(N); iota(perm.begin(), perm.end(), 0);
        do {
            Key kc = 0; for (int i = 0; i < N; i++) kc |= (Key)fields[i] << (12 * perm[i]);
            bool ok = true; for (auto* g : gLab) { int sg; Key kk = applyKey(*g, kc, sg); if (kk < kc) { ok = false; break; } }
            if (!ok) continue;
            Key best; int sg = canonFactored(kc, best); if (sg) out.push_back(best);
        } while (next_permutation(perm.begin(), perm.end()));
    };
    recE = [&](int e, int rem) {
        if (e == E - 1) { lens[e] = rem; process(); lens[e] = 0; return; }
        for (int l = 0; l <= rem; l++) { lens[e] = l; recE(e + 1, rem - l); }
        lens[e] = 0;
    };
    recV = [&](int v, int left) {
        if (v == V) { if (left == 0) { if (E > 0) recE(0, N - k); else if (N == k) process(); } return; }
        if (V - v > left) recV(v + 1, left);       // leave v empty
        if (left > 0) for (int h : vhalf[v]) if (h != base[v]) { occh[v] = h; recV(v + 1, left - 1); occh[v] = -1; }
    };
    recV(0, k);
    sort(out.begin(), out.end()); out.erase(unique(out.begin(), out.end()), out.end());
    fprintf(stderr, "  [enum2] deg %d: shapes %lld, A-canonical %lld, orbit cells %zu\n", k, shapes, canonShapes, out.size());
}
// enumerate all cells of degree k
void enumerate(int k, vector<Key>& out) {
    if (grouped && fastEnum) { enumerate2(k, out); return; }
    Cell c; c.vocc.assign(V, -1); c.vh.assign(V, -1); c.elist.assign(E, {});
    int moving = 0;
    function<void(int)> rec = [&](int p) {
        if (p == N) { if (moving == k) { if (grouped) { Key kc = encode(c); if (isCanonicalKey(kc)) out.push_back(kc); } else if (colored) { Cell d = c; canonicalize(d); out.push_back(encode(d)); } else out.push_back(encode(c)); } return; }
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
    long long pivots = 0, lastCheck = 0; bool goDense = false; static inline long long denseMax = 0; static inline int cntThresh = 48; static inline long long denseMin = 30000, nnzBudget = 120000000; ll maxabs = 1; vector<long long> stages; long long leftover_nonzero = 0;
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
            S.assign(out.begin(), out.end());
        }
        for (auto x : R) colcnt[col(x)]--;
        rowdead[r] = 1; coldead[c] = 1; colcnt[c] = 0; pivots++;
        vector<uint32_t>().swap(rows[r]); vector<int>().swap(colrows[c]);
    }
    bool sweep() {
        bool any = false;
        typedef pair<int,int> PI; priority_queue<PI, vector<PI>, greater<PI>> pq;
        for (int c = 0; c < nc; c++) if (!coldead[c] && colcnt[c] > 0) pq.push({colcnt[c], c});
        long long sinceDiag = 0;
        while (!pq.empty()) { auto [cnt, c] = pq.top(); pq.pop();
            if (pq.size() > 4 * (size_t)nc + 1000000) {   // drop stale queue entries
                vector<PI> keep; keep.reserve(nc); while (!pq.empty()) { auto t = pq.top(); pq.pop(); if (!coldead[t.second] && colcnt[t.second] == t.first && t.first > 0) keep.push_back(t); }
                sort(keep.begin(), keep.end()); keep.erase(unique(keep.begin(), keep.end()), keep.end());
                pq = priority_queue<PI, vector<PI>, greater<PI>>(keep.begin(), keep.end()); }
            if (++sinceDiag >= 200000) { sinceDiag = 0;
                long long nnz = 0, cre = 0, lc = 0; for (int r2 = 0; r2 < nr; r2++) if (!rowdead[r2]) nnz += rows[r2].size();
                for (int c2 = 0; c2 < nc; c2++) { cre += colrows[c2].size(); if (!coldead[c2] && colcnt[c2] > 0) lc++; }
                fprintf(stderr, "  [elim] pivots=%lld rowNNZ=%lld colrowsEntries=%lld liveCols=%lld pq=%zu cnt=%d\n", pivots, nnz, cre, lc, pq.size(), cnt);
                if (cre > nnz + nnz / 2 + 1000000) {   // compact stale column lists
                    for (int c2 = 0; c2 < nc; c2++) if (!coldead[c2] && colrows[c2].size() > (size_t)colcnt[c2]) {
                        vector<int> M; M.reserve(colcnt[c2]); for (int r : colrows[c2]) if (!rowdead[r] && get(r, c2)) M.push_back(r);
                        sort(M.begin(), M.end()); M.erase(unique(M.begin(), M.end()), M.end()); colrows[c2].assign(M.begin(), M.end()); colcnt[c2] = M.size(); } } }
            if (coldead[c] || colcnt[c] != cnt || cnt == 0) continue;
            if (colrows[c].size() > 2 * (size_t)colcnt[c] + 8) { vector<int> M; for (int r : colrows[c]) if (!rowdead[r] && get(r, c)) M.push_back(r); sort(M.begin(), M.end()); M.erase(unique(M.begin(), M.end()), M.end()); colrows[c].swap(M); }
            int best = -1; size_t bl = SIZE_MAX;
            for (int r : colrows[c]) if (!rowdead[r] && rows[r].size() < bl) { uint32_t v = get(r, c); if (v && unit(v)) { bl = rows[r].size(); best = r; } }
            if (best < 0) continue;
            vector<int> touched; for (auto x : rows[best]) touched.push_back(col(x));
            pivot(best, c); any = true;
            if (denseMax > 0 && cnt >= cntThresh && pivots - lastCheck >= 2000) { lastCheck = pivots;
                long long lc = 0; for (int c2 = 0; c2 < nc; c2++) if (!coldead[c2] && colcnt[c2] > 0) lc++;
                long long nz = 0; if (lc <= denseMax && lc > denseMin) for (int r2 = 0; r2 < nr; r2++) if (!rowdead[r2]) nz += rows[r2].size();
                if (lc <= denseMax && (lc <= denseMin || nz > nnzBudget)) { goDense = true; fprintf(stderr, "  [dense] switching: pivots=%lld liveCols=%lld colcnt=%d\n", pivots, lc, cnt); return true; } }
            for (int c2 : touched) if (!coldead[c2] && colcnt[c2] > 0) pq.push({colcnt[c2], c2});
        }
        return any;
    }

    // ---------- dense phase: streaming RREF over Z/mod with unit pivots, stage by stage ----------
    template <uint32_t MD> static void axpy(uint8_t* __restrict x, const uint8_t* __restrict b, uint32_t f, int m) {
        // x <- x - f*b  (mod MD)
        uint32_t g = (MD - f) % MD;
        for (int j = 0; j < m; j++) { uint16_t t = (uint16_t)x[j] + (uint16_t)g * (uint16_t)b[j]; x[j] = (uint8_t)(t % MD); }
    }
    template <uint32_t MD> static void scal(uint8_t* x, uint32_t f, int m) { for (int j = 0; j < m; j++) x[j] = (uint8_t)(((uint16_t)x[j] * f) % MD); }
    void axpyD(uint8_t* x, const uint8_t* b, uint32_t f, int m) {
        switch (mod) {
#define AX(V) case V: axpy<V>(x, b, f, m); break;
            AX(2) AX(4) AX(8) AX(16) AX(3) AX(9) AX(27) AX(5) AX(25) AX(7) AX(11) AX(13) AX(17) AX(19) AX(23) AX(29) AX(31)
#undef AX
            default: fprintf(stderr, "dense: unsupported modulus %u\n", mod); exit(5);
        }
    }
    void scalD(uint8_t* x, uint32_t f, int m) {
        switch (mod) {
#define SC(V) case V: scal<V>(x, f, m); break;
            SC(2) SC(4) SC(8) SC(16) SC(3) SC(9) SC(27) SC(5) SC(25) SC(7) SC(11) SC(13) SC(17) SC(19) SC(23) SC(29) SC(31)
#undef SC
            default: exit(5);
        }
    }

    // ---------- bit-sliced Z/4 dense phase (p0 = 2, mod = 4): two bit planes per row ----------
    static inline int getv(const uint64_t* lo, const uint64_t* hi, int c) { return (int)((lo[c >> 6] >> (c & 63)) & 1) | (int)(((hi[c >> 6] >> (c & 63)) & 1) << 1); }
    // x <- x - f*b  over Z/4, bit-sliced
    static inline void z4axpy(uint64_t* __restrict xl, uint64_t* __restrict xh, const uint64_t* __restrict bl, const uint64_t* __restrict bh, int f, int W) {
        if (f == 2) { for (int w = 0; w < W; w++) xh[w] ^= bl[w]; return; }
        if (f == 1) {   // add -b = (b0, b1 ^ b0)
            for (int w = 0; w < W; w++) { uint64_t a0 = xl[w], a1 = xh[w], c0 = bl[w], c1 = bh[w] ^ bl[w]; xl[w] = a0 ^ c0; xh[w] = a1 ^ c1 ^ (a0 & c0); } return; }
        for (int w = 0; w < W; w++) { uint64_t a0 = xl[w], a1 = xh[w], c0 = bl[w], c1 = bh[w]; xl[w] = a0 ^ c0; xh[w] = a1 ^ c1 ^ (a0 & c0); }   // f == 3: add b
    }
    void dense_z4(long long sparseStagePivots) {
        vector<int> cmap(nc, -1); int m = 0;
        for (int c = 0; c < nc; c++) if (!coldead[c] && colcnt[c] > 0) cmap[c] = m++;
        vector<vector<int>>().swap(colrows);
        int W = (m + 63) / 64;
        fprintf(stderr, "  [z4] m=%d, basis memory <= %.2f GB\n", m, 2.0 * m * W * 8 / 1e9);
        vector<uint64_t*> BL, BH; vector<int> rowOfPiv(m, -1);
        vector<uint64_t> xl(W), xh(W);
        vector<vector<uint64_t>> Yl; vector<int> yp;   // online F2 RREF of halved leftovers (y = x/2 mod 2)
        auto insY = [&](vector<uint64_t>& y) {
            for (size_t i = 0; i < Yl.size(); i++) { int c = yp[i]; if ((y[c >> 6] >> (c & 63)) & 1) for (int w = 0; w < W; w++) y[w] ^= Yl[i][w]; }
            int c = -1; for (int w = 0; w < W && c < 0; w++) if (y[w]) c = w * 64 + __builtin_ctzll(y[w]);
            if (c < 0) return;
            for (size_t i = 0; i < Yl.size(); i++) if ((Yl[i][c >> 6] >> (c & 63)) & 1) for (int w = 0; w < W; w++) Yl[i][w] ^= y[w];
            Yl.push_back(y); yp.push_back(c);
        };
        long long done = 0, stagePiv = sparseStagePivots;
        for (int r = 0; r < nr; r++) {
            if (rowdead[r] || rows[r].empty()) continue;
            fill(xl.begin(), xl.end(), 0); fill(xh.begin(), xh.end(), 0);
            for (auto e : rows[r]) { int c = cmap[col(e)]; if (c < 0) { fprintf(stderr, "z4: dead col\n"); exit(6); } int v = val(e) & 3; if (v & 1) xl[c >> 6] |= 1ULL << (c & 63); if (v & 2) xh[c >> 6] |= 1ULL << (c & 63); }
            for (auto e : rows[r]) { int c = cmap[col(e)]; int br = rowOfPiv[c]; if (br >= 0) { int f = getv(xl.data(), xh.data(), c); if (f) z4axpy(xl.data(), xh.data(), BL[br], BH[br], f, W); } }
            vector<uint32_t>().swap(rows[r]);
            int j = -1; for (int w = 0; w < W; w++) if (xl[w]) { j = w * 64 + __builtin_ctzll(xl[w]); break; }
            if (j >= 0) {
                if (getv(xl.data(), xh.data(), j) == 3) for (int w = 0; w < W; w++) { xh[w] ^= xl[w]; }   // negate: (x0, x1 ^ x0)
                long long nb = BL.size(); const uint64_t* pl = xl.data(); const uint64_t* ph = xh.data();
                #pragma omp parallel for schedule(static, 64)
                for (long long b = 0; b < nb; b++) { int f = getv(BL[b], BH[b], j); if (f) z4axpy(BL[b], BH[b], pl, ph, f, W); }
                uint64_t* nl = (uint64_t*)malloc(W * 8); uint64_t* nh = (uint64_t*)malloc(W * 8);
                if (!nl || !nh) { fprintf(stderr, "z4: out of memory\n"); exit(9); }
                memcpy(nl, xl.data(), W * 8); memcpy(nh, xh.data(), W * 8);
                rowOfPiv[j] = BL.size(); BL.push_back(nl); BH.push_back(nh); stagePiv++; pivots++; rowdead[r] = 1;
            } else {
                bool nz = false; for (int w = 0; w < W; w++) if (xh[w]) { nz = true; break; }
                if (nz) { vector<uint64_t> y(xh.begin(), xh.end()); insY(y); }
            }
            if (++done % 500000 == 0) fprintf(stderr, "  [z4] %lld rows, basis %zu, leftover span %zu\n", done, BL.size(), Yl.size());
        }
        stages.push_back(stagePiv);
        fprintf(stderr, "  [z4] stage 0 done: basis %zu, leftover span %zu\n", BL.size(), Yl.size());
        // final reduction of halved leftovers against the stage-0 basis mod 2 (lo planes), then F2 rank
        vector<vector<uint64_t>> F; int W2 = W; int m2 = m;
        for (auto& y : Yl) {
            vector<int> hits; for (int w = 0; w < W; w++) { uint64_t t = y[w]; while (t) { int c = w * 64 + __builtin_ctzll(t); t &= t - 1; if (rowOfPiv[c] >= 0) hits.push_back(c); } }
            for (int c : hits) { const uint64_t* bl = BL[rowOfPiv[c]]; for (int w = 0; w < W; w++) y[w] ^= bl[w]; }
            bool nz = false; for (int w = 0; w < W; w++) { uint64_t t = y[w]; while (t) { int c = w * 64 + __builtin_ctzll(t); t &= t - 1; if (rowOfPiv[c] >= 0) { fprintf(stderr, "z4: pivot col nonzero\n"); exit(7); } nz = true; } }
            if (nz) F.push_back(move(y));
        }
        vector<vector<uint64_t>>().swap(Yl);
        for (auto p : BL) free(p); for (auto p : BH) free(p);
        // F2 rank of the halved leftovers
        long long r1 = 0;
        for (int c = 0; c < m2 && r1 < (long long)F.size(); c++) {
            int w = c >> 6; uint64_t bit = 1ULL << (c & 63); long long piv = -1;
            for (long long i = r1; i < (long long)F.size(); i++) if (F[i][w] & bit) { piv = i; break; }
            if (piv < 0) continue;
            swap(F[piv], F[r1]);
            for (long long i = 0; i < (long long)F.size(); i++) if (i != r1 && (F[i][w] & bit)) for (int t = 0; t < W2; t++) F[i][t] ^= F[r1][t];
            r1++;
        }
        if (r1) { stages.push_back(r1); pivots += r1; }
        fprintf(stderr, "  [z4] stage 1 (Z/2 summands): %lld\n", r1);
    }

    // ---------- bit-sliced Z/9 dense phase (p0 = 3, mod = 9): four bit planes, values 0..8 ----------
    static void f3axpy(vector<uint64_t>& A1, vector<uint64_t>& A2, const vector<uint64_t>& B1v, const vector<uint64_t>& B2v, int f, int W) {
        for (int w = 0; w < W; w++) {
            uint64_t a1 = A1[w], a2 = A2[w], b1 = B1v[w], b2 = B2v[w]; if (f == 2) swap(b1, b2);
            A1[w] = (a1 & ~(b1 | b2)) | (b1 & ~(a1 | a2)) | (a2 & b2);
            A2[w] = (a2 & ~(b1 | b2)) | (b2 & ~(a1 | a2)) | (a1 & b1);
        }
    }
    struct R9 { uint64_t* p[4]; };
    static inline void add9w(uint64_t& x0, uint64_t& x1, uint64_t& x2, uint64_t& x3, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3) {
        uint64_t s0 = x0 ^ y0, c = x0 & y0;
        uint64_t t = x1 ^ y1; uint64_t s1 = t ^ c; c = (x1 & y1) | (c & t);
        t = x2 ^ y2; uint64_t s2 = t ^ c; c = (x2 & y2) | (c & t);
        t = x3 ^ y3; uint64_t s3 = t ^ c; uint64_t s4 = (x3 & y3) | (c & t);
        uint64_t m = s4 | (s3 & (s0 | s1 | s2));
        // add 23 = 10111b under mask m (i.e. subtract 9 mod 32)
        uint64_t r0 = s0 ^ m; c = s0 & m;
        t = s1 ^ m; uint64_t r1 = t ^ c; c = (s1 & m) | (c & t);
        t = s2 ^ m; uint64_t r2 = t ^ c; c = (s2 & m) | (c & t);
        uint64_t r3 = s3 ^ c;
        x0 = r0; x1 = r1; x2 = r2; x3 = r3;
    }
    static inline void add9row(uint64_t* const* X, uint64_t* const* Y, int W) {
        uint64_t *a = X[0], *b = X[1], *c = X[2], *d = X[3]; const uint64_t *e = Y[0], *f = Y[1], *g = Y[2], *h = Y[3];
        for (int w = 0; w < W; w++) add9w(a[w], b[w], c[w], d[w], e[w], f[w], g[w], h[w]);
    }
    static inline int get9(uint64_t* const* X, int c) { int w = c >> 6, b = c & 63; return (int)((X[0][w] >> b) & 1) | (int)(((X[1][w] >> b) & 1) << 1) | (int)(((X[2][w] >> b) & 1) << 2) | (int)(((X[3][w] >> b) & 1) << 3); }
    static inline void set9(uint64_t* const* X, int c, int v) { int w = c >> 6; uint64_t bit = 1ULL << (c & 63); for (int i = 0; i < 4; i++) { if ((v >> i) & 1) X[i][w] |= bit; else X[i][w] &= ~bit; } }
    // multiples g*Y for g = 0..8 into M (9 rows of 4 planes); M[0] = 0
    static void multiples9(uint64_t* const* Y, vector<vector<uint64_t>>& M, int W) {
        M.assign(36, vector<uint64_t>(W, 0));
        auto P = [&](int g, int i) { return M[4 * g + i].data(); };
        for (int i = 0; i < 4; i++) memcpy(P(1, i), Y[i], W * 8);
        for (int g = 2; g <= 8; g++) {
            for (int i = 0; i < 4; i++) memcpy(P(g, i), P(g - 1, i), W * 8);
            uint64_t* X[4] = {P(g, 0), P(g, 1), P(g, 2), P(g, 3)}; uint64_t* Yy[4] = {P(1, 0), P(1, 1), P(1, 2), P(1, 3)};
            add9row(X, Yy, W);
        }
    }
    void dense_z9(long long sparseStagePivots) {
        vector<int> cmap(nc, -1); int m = 0;
        for (int c = 0; c < nc; c++) if (!coldead[c] && colcnt[c] > 0) cmap[c] = m++;
        vector<vector<int>>().swap(colrows);
        int W = (m + 63) / 64;
        fprintf(stderr, "  [z9] m=%d, basis memory <= %.2f GB\n", m, 4.0 * m * W * 8 / 1e9);
        vector<R9> B; vector<int> rowOfPiv(m, -1);
        vector<uint64_t> xs(4 * (size_t)W); uint64_t* X[4] = {xs.data(), xs.data() + W, xs.data() + 2 * W, xs.data() + 3 * W};
        vector<uint64_t> accs(36 * (size_t)W);   // per-coefficient accumulators S_g, g = 1..8
        vector<vector<uint64_t>> Mx;
        vector<vector<uint64_t>> Ya, Yb; vector<int> ypiv;
        auto insertY = [&](vector<uint64_t>& a, vector<uint64_t>& b) {
            for (size_t i = 0; i < Ya.size(); i++) { int c = ypiv[i]; int w = c >> 6; uint64_t bit = 1ULL << (c & 63); int f = (a[w] & bit) ? 1 : ((b[w] & bit) ? 2 : 0); if (f) f3axpy(a, b, Ya[i], Yb[i], 3 - f, W); }
            int c = -1; for (int w = 0; w < W && c < 0; w++) { uint64_t q = a[w] | b[w]; if (q) c = w * 64 + __builtin_ctzll(q); }
            if (c < 0) return;
            if ((b[c >> 6] >> (c & 63)) & 1) swap(a, b);   // normalise pivot to 1 (multiply by 2 = swap planes)
            for (size_t i = 0; i < Ya.size(); i++) { int w = c >> 6; uint64_t bit = 1ULL << (c & 63); int f = (Ya[i][w] & bit) ? 1 : ((Yb[i][w] & bit) ? 2 : 0); if (f) f3axpy(Ya[i], Yb[i], a, b, 3 - f, W); }
            Ya.push_back(a); Yb.push_back(b); ypiv.push_back(c);
        };
        vector<vector<pair<int,int>>> left;
        long long done = 0, stagePiv = sparseStagePivots;
        auto reduce = [&](const vector<pair<int,int>>& ents) {
            fill(xs.begin(), xs.end(), 0);
            for (auto& [c, v] : ents) set9(X, c, v);
            // group basis rows by the coefficient (9 - x[c]) and accumulate
            bool used[9] = {false};
            for (auto& [c, v] : ents) { int br = rowOfPiv[c]; if (br < 0) continue; int f = get9(X, c); if (!f) continue; int g = (9 - f) % 9;
                uint64_t* A[4] = {&accs[(4 * g) * (size_t)W], &accs[(4 * g + 1) * (size_t)W], &accs[(4 * g + 2) * (size_t)W], &accs[(4 * g + 3) * (size_t)W]};
                if (!used[g]) { for (int i = 0; i < 4; i++) memcpy(A[i], B[br].p[i], W * 8); used[g] = true; } else add9row(A, B[br].p, W); }
            // x += sum_g g * S_g   (multiply via repeated doubling-free additions: g copies)
            for (int g = 1; g <= 8; g++) if (used[g]) {
                uint64_t* A[4] = {&accs[(4 * g) * (size_t)W], &accs[(4 * g + 1) * (size_t)W], &accs[(4 * g + 2) * (size_t)W], &accs[(4 * g + 3) * (size_t)W]};
                for (int t = 0; t < g; t++) add9row(X, A, W);
            }
        };
        auto unitcol = [&]() { for (int w = 0; w < W; w++) { uint64_t nz = X[0][w] | X[1][w] | X[2][w] | X[3][w]; while (nz) { int c = w * 64 + __builtin_ctzll(nz); nz &= nz - 1; if (get9(X, c) % 3) return c; } } return -1; };
        static const int inv9[9] = {0, 1, 5, 0, 7, 2, 0, 4, 8};
        for (int r = 0; r < nr; r++) {
            if (rowdead[r] || rows[r].empty()) continue;
            vector<pair<int,int>> ents; ents.reserve(rows[r].size());
            for (auto e : rows[r]) { int c = cmap[col(e)]; if (c < 0) { fprintf(stderr, "z9: dead col\n"); exit(6); } ents.push_back({c, (int)(val(e) % 9)}); }
            vector<uint32_t>().swap(rows[r]);
            reduce(ents);
            int j = unitcol();
            if (j >= 0) {
                int u = inv9[get9(X, j)];
                if (u != 1) { multiples9(X, Mx, W); for (int i = 0; i < 4; i++) memcpy(X[i], Mx[4 * u + i].data(), W * 8); }
                multiples9(X, Mx, W);
                long long nb = B.size();
                #pragma omp parallel for schedule(static, 64)
                for (long long b = 0; b < nb; b++) { int f = get9(B[b].p, j); if (f) { int g = (9 - f) % 9; uint64_t* Y[4] = {Mx[4 * g].data(), Mx[4 * g + 1].data(), Mx[4 * g + 2].data(), Mx[4 * g + 3].data()}; add9row(B[b].p, Y, W); } }
                R9 nr9; for (int i = 0; i < 4; i++) { nr9.p[i] = (uint64_t*)malloc(W * 8); if (!nr9.p[i]) { fprintf(stderr, "z9: out of memory\n"); exit(9); } memcpy(nr9.p[i], X[i], W * 8); }
                rowOfPiv[j] = B.size(); B.push_back(nr9); stagePiv++; pivots++; rowdead[r] = 1;
            } else {
                // leftover: all entries divisible by 3; keep y = x/3 mod 3 in an online F3 RREF (span is all that matters)
                vector<uint64_t> ya(W, 0), yb(W, 0); bool any = false;
                for (int w = 0; w < W; w++) { uint64_t nz = X[0][w] | X[1][w] | X[2][w] | X[3][w]; while (nz) { int c = w * 64 + __builtin_ctzll(nz); nz &= nz - 1; int v = get9(X, c);
                    if (v % 3) { fprintf(stderr, "z9: unit in leftover\n"); exit(7); } if (v == 3) ya[w] |= 1ULL << (c & 63); else yb[w] |= 1ULL << (c & 63); any = true; } }
                if (any) insertY(ya, yb);
            }
            if (++done % 200000 == 0) fprintf(stderr, "  [z9] %lld rows, basis %zu, leftover span %zu\n", done, B.size(), Ya.size());
        }
        stages.push_back(stagePiv);
        fprintf(stderr, "  [z9] stage 0 done: basis %zu, leftover span %zu\n", B.size(), Ya.size());
        // final reduction of each leftover against the stage-0 basis (mod 3), then F3 rank
        map<int, pair<vector<uint64_t>, vector<uint64_t>>> bcache;
        auto bmod3 = [&](int br) -> pair<vector<uint64_t>, vector<uint64_t>>& {
            auto it = bcache.find(br); if (it != bcache.end()) return it->second;
            vector<uint64_t> a(W, 0), b(W, 0);
            for (int w = 0; w < W; w++) { uint64_t nz = B[br].p[0][w] | B[br].p[1][w] | B[br].p[2][w] | B[br].p[3][w]; while (nz) { int c = w * 64 + __builtin_ctzll(nz); nz &= nz - 1; int v = get9(B[br].p, c) % 3; if (v == 1) a[w] |= 1ULL << (c & 63); else if (v == 2) b[w] |= 1ULL << (c & 63); } }
            return bcache[br] = {move(a), move(b)}; };
        vector<vector<uint64_t>> F0, F1;
        for (size_t i = 0; i < Ya.size(); i++) {
            vector<uint64_t> a = Ya[i], b = Yb[i];
            vector<pair<int,int>> hits;   // (pivot column, value) - RREF basis: values at pivot columns are unaffected by other subtractions
            for (int w = 0; w < W; w++) { uint64_t q = a[w] | b[w]; while (q) { int c = w * 64 + __builtin_ctzll(q); q &= q - 1; if (rowOfPiv[c] >= 0) hits.push_back({c, ((a[w] >> (c & 63)) & 1) ? 1 : 2}); } }
            for (auto& [c, f] : hits) { auto& bb = bmod3(rowOfPiv[c]); f3axpy(a, b, bb.first, bb.second, 3 - f, W); }
            for (int w = 0; w < W; w++) { uint64_t q = a[w] | b[w]; while (q) { int c = w * 64 + __builtin_ctzll(q); q &= q - 1; if (rowOfPiv[c] >= 0) { fprintf(stderr, "z9: pivot col nonzero after final reduction\n"); exit(7); } } }
            bool nz = false; for (int w = 0; w < W && !nz; w++) if (a[w] | b[w]) nz = true;
            if (nz) { F0.push_back(move(a)); F1.push_back(move(b)); }
        }
        bcache.clear();
        for (auto& r9 : B) for (int i = 0; i < 4; i++) free(r9.p[i]);
        long long nrw = F0.size(), r1 = 0;
        auto getf = [&](long long i, int c) { int w = c >> 6; uint64_t bit = 1ULL << (c & 63); return (F0[i][w] & bit) ? 1 : ((F1[i][w] & bit) ? 2 : 0); };
        for (int c = 0; c < m && r1 < nrw; c++) {
            long long piv = -1; for (long long i = r1; i < nrw; i++) if (getf(i, c)) { piv = i; break; }
            if (piv < 0) continue;
            swap(F0[piv], F0[r1]); swap(F1[piv], F1[r1]);
            int pv = getf(r1, c);
            for (long long i = 0; i < nrw; i++) if (i != r1) { int v = getf(i, c); if (v) f3axpy(F0[i], F1[i], F0[r1], F1[r1], (v == pv) ? 2 : 1, W); }
            r1++;
        }
        if (r1) { stages.push_back(r1); pivots += r1; }
        fprintf(stderr, "  [z9] stage 1 (Z/3 summands): %lld\n", r1);
    }
    vector<int> denseRowOrigin;   // which original row each stage-0 dense pivot came from
    void dense_finish(long long sparseStagePivots) {
        if (p0 == 2 && mod == 4) { dense_z4(sparseStagePivots); return; }
        if (p0 == 3 && mod == 9) { dense_z9(sparseStagePivots); return; }
        vector<int> cmap(nc, -1); int m = 0;
        for (int c = 0; c < nc; c++) if (!coldead[c] && colcnt[c] > 0) cmap[c] = m++;
        vector<vector<int>>().swap(colrows);
        // input rows (sparse, original-row index kept for stage-0 bookkeeping)
        vector<pair<int, vector<pair<int,uint8_t>>>> in;
        for (int r = 0; r < nr; r++) if (!rowdead[r] && !rows[r].empty()) {
            vector<pair<int,uint8_t>> v; v.reserve(rows[r].size());
            for (auto x : rows[r]) { int cc = cmap[col(x)]; if (cc < 0) { fprintf(stderr, "dense: entry in dead column\n"); exit(6); } v.push_back({cc, (uint8_t)val(x)}); }
            vector<uint32_t>().swap(rows[r]); in.push_back({r, move(v)});
        }
        fprintf(stderr, "  [dense] m=%d input rows=%zu, est. basis memory <= %.2f GB\n", m, in.size(), (double)m * m / 1e9);
        long long stagePiv = sparseStagePivots; bool first = true;
        while (true) {
            vector<uint8_t*> B; vector<int> rowOfPiv(m, -1); vector<int> pivOfRow;
            vector<uint8_t> x(m);
            vector<vector<pair<int,uint8_t>>> left;
            auto reduce = [&](const vector<pair<int,uint8_t>>& sp) {
                fill(x.begin(), x.end(), 0);
                for (auto& [c, v] : sp) x[c] = v;
                for (auto& [c, v] : sp) { int br = rowOfPiv[c]; if (br >= 0 && x[c]) axpyD(x.data(), B[br], x[c], m); }
            };
            long long done = 0;
            for (auto& [orig, sp] : in) {
                reduce(sp);
                int j = -1;
                for (int t0 = 0; t0 < m && j < 0; t0 += 8) {
                    uint64_t w = 0; if (t0 + 8 <= m) memcpy(&w, x.data() + t0, 8); else w = 1;
                    if (!w) continue;
                    for (int t = t0; t < min(m, t0 + 8); t++) if (x[t] % p0) { j = t; break; }
                }
                if (j >= 0) {
                    scalD(x.data(), inv(x[j]), m);
                    { long long nb2 = B.size(); uint8_t* xd = x.data();
                    #pragma omp parallel for schedule(static, 64)
                    for (long long b2 = 0; b2 < nb2; b2++) if (B[b2][j]) axpyD(B[b2], xd, B[b2][j], m); }
                    uint8_t* nb = (uint8_t*)malloc(m); memcpy(nb, x.data(), m);
                    rowOfPiv[j] = B.size(); pivOfRow.push_back(j); B.push_back(nb);
                    stagePiv++; pivots++;
                    if (first) { if (orig >= 0) rowdead[orig] = 1; }
                } else {
                    vector<pair<int,uint8_t>> rest; for (int t = 0; t < m; t++) if (x[t]) rest.push_back({t, x[t]});
                    if (!rest.empty()) left.push_back(move(rest));
                }
                vector<pair<int,uint8_t>>().swap(sp);
                if (++done % 200000 == 0) fprintf(stderr, "  [dense] mod %u: %lld rows, basis %zu, leftover %zu\n", mod, done, B.size(), left.size());
            }
            stages.push_back(stagePiv); stagePiv = 0; first = false;
            fprintf(stderr, "  [dense] stage done: mod %u basis %zu leftover rows %zu\n", mod, B.size(), left.size());
            // re-reduce leftovers against the final basis; survivors are divisible by p0
            vector<pair<int, vector<pair<int,uint8_t>>>> nxt;
            vector<int> cmap2(m, -1); int m2 = 0; for (int c = 0; c < m; c++) if (rowOfPiv[c] < 0) cmap2[c] = m2++;
            for (auto& sp : left) {
                reduce(sp);
                vector<pair<int,uint8_t>> rest;
                for (int t = 0; t < m; t++) if (x[t]) { if (x[t] % p0) { fprintf(stderr, "dense: unit survived\n"); exit(7); } if (rowOfPiv[t] >= 0) { fprintf(stderr, "dense: pivot col nonzero\n"); exit(7); } rest.push_back({cmap2[t], (uint8_t)(x[t])}); }
                if (!rest.empty()) nxt.push_back({-1, move(rest)});
            }
            for (auto p : B) free(p);
            if (nxt.empty()) break;
            if (mod == p0) { for (auto& r : nxt) leftover_nonzero += r.second.size(); break; }
            for (auto& r : nxt) { vector<pair<int,uint8_t>> M2; for (auto& [c, v] : r.second) if (v / p0) M2.push_back({c, (uint8_t)(v / p0)}); r.second.swap(M2); }
            mod /= p0; m = m2; in.swap(nxt);
            in.erase(remove_if(in.begin(), in.end(), [](auto& r) { return r.second.empty(); }), in.end());
            if (in.empty()) break;
        }
    }
    void run() {
        long long before = 0;
        if (denseMax < 0) { dense_finish(0); return; }
        while (true) {
            while (sweep()) { if (goDense) break; }
            if (goDense) { dense_finish(pivots - before); return; }
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
        else if (a == "-D") LPE::denseMax = atoll(argv[++i]);
        else if (a == "-F") fastEnum = true;
        else if (a == "-Y") { string t = argv[++i]; stringstream ss(t); string x; while (getline(ss, x, ',')) agSigns.push_back(stoi(x)); }
        else if (a == "-M") { LPE::denseMin = atoll(argv[++i]); LPE::nnzBudget = atoll(argv[++i]); }
        else if (a == "-X") { string t = argv[++i]; twistLab = t.find('l') != string::npos; twistEdge = t.find('e') != string::npos; }
        else if (a == "-a") { grouped = true; string all = argv[++i]; stringstream gs(all); string one; while (getline(gs, one, ';')) agens.push_back(one); }
        else if (a == "-C") LPE::cntThresh = atoi(argv[++i]);
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
    etail.resize(E); ehead.resize(E);
    hvert.resize(2 * E); vhalf.assign(V, {});
    for (int e = 0; e < E; e++) {
        scanf("%d %d", &etail[e], &ehead[e]);
        hvert[2 * e] = etail[e]; hvert[2 * e + 1] = ehead[e];
        vhalf[etail[e]].push_back(2 * e); vhalf[ehead[e]].push_back(2 * e + 1);
    }
    if (3 * E >= 255 || N > 10 || N > 15) { fprintf(stderr, "too large\n"); return 1; }
    base.assign(V, -1);
    if (!grouped) { for (int v = 0; v < V; v++) if (!vhalf[v].empty()) base[v] = vhalf[v][0]; }
    else {
        auto closeup = [](vector<vector<int>> gens, int n) {   // closure of perms of {0..n-1}
            vector<int> id(n); iota(id.begin(), id.end(), 0);
            set<vector<int>> G = {id}; vector<vector<int>> fr = {id};
            while (!fr.empty()) { vector<vector<int>> nf;
                for (auto& x : fr) for (auto& g : gens) { vector<int> y(n); for (int t = 0; t < n; t++) y[t] = g[x[t]]; if (G.insert(y).second) nf.push_back(y); }
                fr.swap(nf); }
            return vector<vector<int>>(G.begin(), G.end()); };
        for (auto& g : gelts) for (int t = g.size(); t < N; t++) g.push_back(t);
        auto Lall = closeup(gelts, N);
        // automorphisms: "v0,v1,...[|e0,e1,...]" ; combined vertex+edge perm as one perm of V+E points
        vector<vector<int>> ag;
        for (auto& str : agens) {
            string vs = str, es; size_t bar = str.find('|'); if (bar != string::npos) { vs = str.substr(0, bar); es = str.substr(bar + 1); }
            vector<int> vm, em; vector<int> fl; { stringstream ss(vs); string t; while (getline(ss, t, ',')) vm.push_back(stoi(t)); }
            if ((int)vm.size() != V) { fprintf(stderr, "bad automorphism\n"); exit(8); }
            if (!es.empty()) { stringstream ss(es); string t; while (getline(ss, t, ',')) { bool f = t[0] == '~'; em.push_back(stoi(f ? t.substr(1) : t)); fl.push_back(f); } }
            else {   // match parallel copies in order, allowing reversal
                map<pair<int,int>, vector<int>> byend; for (int e = 0; e < E; e++) byend[{etail[e], ehead[e]}].push_back(e);
                map<pair<int,int>, int> seen; em.assign(E, -1); fl.assign(E, 0);
                for (int e = 0; e < E; e++) { int ci = seen[{etail[e], ehead[e]}]++; pair<int,int> key = {vm[etail[e]], vm[ehead[e]]};
                    if (!byend.count(key)) { key = {key.second, key.first}; fl[e] = 1; }
                    auto& L = byend[key]; if (ci >= (int)L.size()) { fprintf(stderr, "automorphism does not map edges to edges\n"); exit(8); } em[e] = L[ci]; }
            }
            for (int e = 0; e < E; e++) { int a = fl[e] ? ehead[em[e]] : etail[em[e]], b = fl[e] ? etail[em[e]] : ehead[em[e]]; if (a != vm[etail[e]] || b != vm[ehead[e]]) { fprintf(stderr, "automorphism does not preserve edge %d\n", e); exit(8); } }
            vector<int> comb(V + 2 * E); for (int v = 0; v < V; v++) comb[v] = vm[v];
            for (int e = 0; e < E; e++) { comb[V + 2 * e] = V + 2 * em[e] + fl[e]; comb[V + 2 * e + 1] = V + 2 * em[e] + 1 - fl[e]; }
            ag.push_back(comb);
        }
        auto Aall = closeup(ag, V + 2 * E);
        // optional character on the automorphism group from generator signs (-Y); checked for consistency
        map<vector<int>, int> achi;
        { vector<int> id(V + 2 * E); iota(id.begin(), id.end(), 0); achi[id] = 1; vector<vector<int>> fr = {id};
          while (!fr.empty()) { vector<vector<int>> nf;
            for (auto& x : fr) for (size_t gi = 0; gi < ag.size(); gi++) { auto& g = ag[gi]; vector<int> y(V + 2 * E); for (int t = 0; t < V + 2 * E; t++) y[t] = g[x[t]];
                int sy = achi[x] * (gi < agSigns.size() ? agSigns[gi] : 1);
                auto it = achi.find(y); if (it == achi.end()) { achi[y] = sy; nf.push_back(y); } else if (it->second != sy) { fprintf(stderr, "-Y signs do not define a character\n"); exit(12); } }
            fr.swap(nf); } }
        // invariant base half-edges
        for (int v = 0; v < V; v++) if (base[v] < 0 && !vhalf[v].empty()) {
            bool ok = false;
            for (int h : vhalf[v]) {
                map<int,int> asg; bool good = true;
                for (auto& A : Aall) { int w = A[v]; int hh = A[V + h] - V; auto it = asg.find(w); if (it != asg.end() && it->second != hh) { good = false; break; } asg[w] = hh; }
                if (good) { for (auto& [w, hh] : asg) base[w] = hh; ok = true; break; }
            }
            if (!ok) { fprintf(stderr, "no invariant base half-edge at %d\n", v); exit(8); }
        }
        for (auto& L : Lall) for (auto& A : Aall) {
            GE g; g.lab = L; g.vm.assign(A.begin(), A.begin() + V); g.em.resize(E); g.fl.resize(E);
            for (int e = 0; e < E; e++) { int x = A[V + 2 * e] - V; g.em[e] = x / 2; g.fl[e] = x & 1; }
            bool id = true; for (int t = 0; t < N; t++) if (g.lab[t] != t) id = false; for (int v = 0; v < V; v++) if (g.vm[v] != v) id = false; for (int e = 0; e < E; e++) if (g.em[e] != e || g.fl[e]) id = false;
            for (int h = 0; h < 2 * E; h++) g.hm[h] = 2 * g.em[h / 2] + ((h & 1) ^ g.fl[h / 2]);
            g.chi = (twistLab ? permsign(g.lab) : 1) * (twistEdge ? permsign(g.em) : 1) * achi[A];
            if (!id) gall.push_back(g);
        }
        { mt19937 rg(12345); shuffle(gall.begin(), gall.end(), rg); }
        for (auto& g : gall) { bool labid = true, autid = true; for (int t = 0; t < N; t++) if (g.lab[t] != t) labid = false;
            for (int v = 0; v < V; v++) if (g.vm[v] != v) autid = false; for (int e = 0; e < E; e++) if (g.em[e] != e || g.fl[e]) autid = false;
            if (autid) gLab.push_back(&g); else if (labid) gAut.push_back(&g); }
        { GE idg; idg.lab.resize(N); iota(idg.lab.begin(), idg.lab.end(), 0); idg.vm.resize(V); iota(idg.vm.begin(), idg.vm.end(), 0); idg.em.resize(E); iota(idg.em.begin(), idg.em.end(), 0); idg.fl.assign(E, 0); idg.chi = 1;
          for (int h = 0; h < 2 * E; h++) idg.hm[h] = h;
          Aelts.push_back(idg); Lelts.push_back(idg); for (auto* g : gAut) Aelts.push_back(*g); for (auto* g : gLab) Lelts.push_back(*g);
          if (Aelts.size() * Lelts.size() != gall.size() + 1) { fprintf(stderr, "factorisation failed\n"); exit(11); } }
        if (LQ && (long long)(gall.size() + 1) % LQ == 0) { fprintf(stderr, "group order %zu divisible by p: coinvariants invalid\n", gall.size() + 1); exit(10); }
        fprintf(stderr, "group order %zu (labels %zu x automorphisms %zu); bases:", gall.size() + 1, Lall.size(), Aall.size());
        for (int v = 0; v < V; v++) fprintf(stderr, " %d", base[v]); fprintf(stderr, "\n");
    }
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
                    int sf = canonicalize(f); if (!sf) continue;
                    Key fk = encode(f);
                    auto it = lower_bound(C[k-1].begin(), C[k-1].end(), fk);
                    if (it == C[k-1].end() || *it != fk) { fprintf(stderr, "face not found\n"); exit(2); }
                    out.push_back({(long long)(it - C[k-1].begin()), (t == 0 ? sgn : -sgn) * sf});
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
            #pragma omp parallel for schedule(dynamic, 4096)
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
                        int sf = canonicalize(f); if (!sf) continue;
                        Key fk = encode(f);
                        auto it = lower_bound(prevcells.begin(), prevcells.end(), fk);
                        if (it == prevcells.end() || *it != fk) { fprintf(stderr, "face not found\n"); exit(2); }
                        ll coef = (t == 0 ? sgn : -sgn) * sf;
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
