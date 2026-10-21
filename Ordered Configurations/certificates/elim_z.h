// =====================================================================================================
//  elim_z.h  --  exact sparse elimination over Z with pivots +-1 only; used by exactsolve and reducecyc
// =====================================================================================================
//  Pivoting only on entries +-1 keeps everything integral: every row operation row_t <- row_t - m * row_s has
//  an integer multiplier m, and the square submatrix (pivot rows x pivot columns) of the ORIGINAL matrix has
//  determinant +-1.  Matrix entries are stored as 32-bit integers (to save memory); every update is computed in
//  128 bits and the program stops if an entry would leave [-2^30, 2^30] (it never came close in our runs).  Vectors that
//  are combined with the matrix (targets, multipliers) are kept in 128 bits, bounded by 2^100.
//  The same records as in elim_modq.h are kept: the operation log, and the order of the pivots (column, row).
//  Facts (1)-(4) listed there hold here as well, exactly over Z.
// =====================================================================================================
#pragma once
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef __int128 I;

const I LIM = (I)1 << 100;
void chk(I v) { if (v > LIM || v < -LIM) { fprintf(stderr, "overflow beyond 2^100\n"); exit(7); } }
void chk64(I v) { if (v > ((I)1 << 30) || v < -((I)1 << 30)) { fprintf(stderr, "entry overflow in elimination (|entry| > 2^30)\n"); exit(7); } }
typedef int32_t ent;                           // storage type of the matrix entries
string i2s(I v) {
    if (v == 0) return "0";
    bool neg = v < 0; if (neg) v = -v;
    string s; while (v) { s += char('0' + (int)(v % 10)); v /= 10; }
    if (neg) s += '-';
    reverse(s.begin(), s.end());
    return s;
}
// Parse a decimal integer (possibly negative, possibly beyond 64 bits) into 128 bits.
I s2i(const string &s) {
    I v = 0; size_t i = 0; bool neg = false;
    if (s[0] == '-') { neg = true; i = 1; }
    for (; i < s.size(); i++) v = v * 10 + (s[i] - '0');
    return neg ? -v : v;
}

struct OpZ { int t, s; ent m; };               // row_t <- row_t - m * row_s

struct ZElim {
    int nr, nc;
    vector<vector<pair<int, ent>>> rows;
    vector<vector<int>> colrows;
    vector<int> colcnt;                       // number of entries +-1 per live column
    vector<char> rowdead, coldead;
    vector<OpZ> ops;
    vector<pair<int, int>> piv;               // (column, row) in pivot order; the row is frozen from then on
    ZElim(int r, int c) : nr(r), nc(c), rows(r), colrows(c), colcnt(c, 0), rowdead(r, 0), coldead(c, 0) {}
    static bool unit(ll v) { return v == 1 || v == -1; }
    void finalize() {
        for (int r = 0; r < nr; r++) {
            auto &R = rows[r];
            sort(R.begin(), R.end(), [](auto &a, auto &b) { return a.first < b.first; });
            vector<pair<int, ent>> M;
            for (auto &x : R) { if (!M.empty() && M.back().first == x.first) M.back().second += x.second; else M.push_back(x); }
            R.clear();
            for (auto &x : M) if (x.second) R.push_back(x);
            for (auto &x : R) { colrows[x.first].push_back(r); if (unit(x.second)) colcnt[x.first]++; }
        }
    }
    // Entry (r, c) of the sparse matrix, found by binary search in row r; 0 if absent.
    ll get(int r, int c) const {
        auto &R = rows[r];
        auto it = lower_bound(R.begin(), R.end(), c, [](const pair<int, ent> &a, int c) { return a.first < c; });
        return (it != R.end() && it->first == c) ? it->second : 0;
    }
    // Pivot on the usable entry (r, c): subtract the right multiple of row r from every other live row with an entry in column c, then mark row r and column c as dead.
    void pivot(int r, int c) {
        auto &R = rows[r];
        ll pv = get(r, c);                    // +-1, so pv is its own inverse
        vector<int> targets;
        for (int r2 : colrows[c]) if (r2 != r && !rowdead[r2]) targets.push_back(r2);
        sort(targets.begin(), targets.end());
        targets.erase(unique(targets.begin(), targets.end()), targets.end());
        for (int r2 : targets) {
            ll a = get(r2, c);
            if (!a) continue;
            ll f = a * pv;
            chk64(f); ops.push_back({r2, r, (ent)f});
            auto &S = rows[r2];
            vector<pair<int, ent>> out;
            out.reserve(S.size() + R.size());
            size_t i = 0, j = 0;
            while (i < S.size() || j < R.size()) {
                if (j == R.size() || (i < S.size() && S[i].first < R[j].first)) { out.push_back(S[i]); i++; }
                else if (i == S.size() || R[j].first < S[i].first) {
                    I v = -(I)f * R[j].second; chk64(v);
                    out.push_back({R[j].first, (ent)v}); colrows[R[j].first].push_back(r2); if (unit((ll)v)) colcnt[R[j].first]++;
                    j++;
                } else {
                    ll old = S[i].second; I v = (I)old - (I)f * R[j].second; chk64(v);
                    if (unit(old)) colcnt[S[i].first]--;
                    if (v) { out.push_back({S[i].first, (ent)v}); if (unit((ll)v)) colcnt[S[i].first]++; }
                    i++; j++;
                }
            }
            S.swap(out);
        }
        for (auto &x : R) if (unit(x.second)) colcnt[x.first]--;
        rowdead[r] = 1; coldead[c] = 1; colcnt[c] = 0;
        piv.push_back({c, r});
        vector<int>().swap(colrows[c]);       // the pivot row itself is kept: it is the snapshot
    }
    void run() {                              // same Markowitz order as QElim::run, with +-1 entries
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
                    vector<int> M;
                    for (int r : colrows[c]) if (!rowdead[r] && get(r, c)) M.push_back(r);
                    sort(M.begin(), M.end()); M.erase(unique(M.begin(), M.end()), M.end());
                    colrows[c].swap(M);
                }
                int best = -1; size_t bl = SIZE_MAX;
                for (int r : colrows[c]) if (!rowdead[r] && rows[r].size() < bl && unit(get(r, c))) { bl = rows[r].size(); best = r; }
                if (best < 0) continue;
                vector<int> touched;
                for (auto &x : rows[best]) touched.push_back(x.first);
                pivot(best, c);
                any = true;
                for (int c2 : touched) if (!coldead[c2] && colcnt[c2] > 0) pq.push({colcnt[c2], c2});
            }
        }
    }
    // Reduce a vector x (over the columns) against the pivot rows, in pivot order.  Afterwards x vanishes on
    // every pivot column (fact (3)).  g[r] receives the multiple of final row r that was subtracted.
    void reduce(vector<I> &x, vector<I> &g) const {
        for (auto &pr : piv) {
            int c = pr.first, r = pr.second;
            I a = x[c];
            if (!a) continue;
            I m = a * (I)get(r, c);           // get(r, c) = +-1
            g[r] += m;
            for (auto &e : rows[r]) { x[e.first] -= m * (I)e.second; chk(x[e.first]); }
        }
    }
    // Replay the log backwards: turns coefficients g on the final rows into coefficients on the original rows.
    void replay_back(vector<I> &g) const {
        for (size_t o = ops.size(); o-- > 0;) {
            auto &op = ops[o];
            if (g[op.t]) { g[op.s] -= (I)op.m * g[op.t]; chk(g[op.s]); }
        }
    }
};
