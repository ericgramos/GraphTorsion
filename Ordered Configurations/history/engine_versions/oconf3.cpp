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
        if (p == N) { if (moving == k) out.push_back(encode(c)); return; }
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

int main(int argc, char** argv) {
    ll P = 0; string prefix = "res"; int maxdeg = -1; bool cellsonly = false;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (a == "-p") P = atoll(argv[++i]);
        else if (a == "-o") prefix = argv[++i];
        else if (a == "-d") maxdeg = atoi(argv[++i]);
        else if (a == "-c") cellsonly = true;
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
    for (int v = 0; v < V; v++) if (!vhalf[v].empty()) base[v] = vhalf[v][0];
    int D = 0;
    for (int v = 0; v < V; v++) if (vhalf[v].size() >= 2) D++;
    D = min(D, N);
    if (maxdeg >= 0) D = min(D, maxdeg);

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
        csize[k] = cells.size();
        if (k >= 1 && !cellsonly) {
            SparseElim M(cells.size(), prevcells.size(), P);
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
                        Key fk = encode(f);
                        auto it = lower_bound(prevcells.begin(), prevcells.end(), fk);
                        if (it == prevcells.end() || *it != fk) { fprintf(stderr, "face not found\n"); return 2; }
                        ll coef = t == 0 ? sgn : -sgn;
                        if (P) coef = ((coef % P) + P) % P;
                        int ci = (int)(it - prevcells.begin());
                        if (!(canDrop && !dropcol.empty() && dropcol[ci])) M.rows[i].push_back({ci, coef});
                    }
                }
            }
            M.finalize();
            M.run();
            int rr, rc; long long rnnz;
            string fn = prefix + "_" + to_string(k) + ".txt";
            M.dump_residual(fn, rr, rc, rnnz);
            char buf[512];
            snprintf(buf, sizeof buf, "{\"k\":%d,\"rows\":%zu,\"cols\":%zu,\"pivots\":%lld,\"res_rows\":%d,\"res_cols\":%d,\"res_nnz\":%lld,\"maxabs\":%lld,\"res_file\":\"%s\"}",
                     k, cells.size(), prevcells.size(), M.pivots, rr, rc, rnnz, M.maxabs, fn.c_str());
            results.push_back(buf);
            // pivot rows of d_k (k-cells) are dropped from d_{k+1} if d_k reduced by unit pivots only
            canDrop = (rr == 0);
            dropcol.assign(cells.size(), 0);
            for (size_t i = 0; i < cells.size(); i++) if (M.rowdead[i]) dropcol[i] = 1;
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
