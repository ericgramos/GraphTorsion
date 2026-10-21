// =====================================================================================================
//  elim_modq.h  --  sparse elimination over Z/Q, Q = p^K, with unit pivots; used by certq and orderq
// =====================================================================================================
//  Over the local ring Z/p^K an entry is a unit exactly when it is prime to p.  Pivoting only on units keeps
//  every operation a reduction of an operation over the localisation Z_(p), so what is computed mod p^K is
//  the reduction of an honest Z_(p) computation.  After all unit pivots are exhausted, the rows and columns
//  that are left alive form a residual block whose entries are all divisible by p; the local Smith form of
//  that block (dense, see local_smith below) gives the p-primary elementary divisors below p^K.
//
//  Two records can be kept for later back-substitution (keep = true):
//    ops  every row operation  row_t <- row_t - m * row_s, in order;
//    piv  for every pivot, its column and a snapshot of the pivot row at the moment it was used.
//  Facts used by the callers.  (1) A dead (pivot) row is never modified again, so its snapshot is also its final
//  value.  (2) A row that is still alive at the end has zero entries in every pivot column.  (3) The snapshot
//  of a pivot has no entries in the columns of EARLIER pivots.  (4) The final matrix is L * D where L is the
//  product of the logged operations, so for a vector g over the final rows, g L (a vector over the original
//  rows) is obtained by replaying the log backwards: for each op (t, s, m), last first, g[s] -= m * g[t].
// =====================================================================================================
#pragma once
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long u64;

u64 Q, PP;                                     // the modulus p^K and the prime p
u64 mulq(u64 a, u64 b) { return (unsigned __int128)a * b % Q; }
bool isunit(u64 a) { return a % PP != 0; }
u64 ppow(int e) { u64 r = 1; while (e--) r *= PP; return r; }
// p-adic valuation of a residue mod Q; 1000 stands for "zero mod Q".
int val(u64 x) { if (x % Q == 0) return 1000; int v = 0; while (x % PP == 0) { x /= PP; v++; } return v; }
// Inverse of a unit mod Q by the extended Euclidean algorithm.
u64 egcd_inv(u64 a) {
    __int128 t = 0, nt = 1, r = Q, nr = a % Q;
    while (nr) {
        __int128 q = r / nr;
        tie(t, nt) = make_pair(nt, t - q * nt);
        tie(r, nr) = make_pair(nr, r - q * nr);
    }
    if (r != 1) { fprintf(stderr, "egcd_inv: not a unit\n"); exit(5); }
    if (t < 0) t += Q;
    return (u64)t;
}
uint32_t modq_ll(long long x) { long long q = (long long)Q; return (uint32_t)(((x % q) + q) % q); }

struct Op { int t, s; uint32_t m; };          // row_t <- row_t - m * row_s

struct QElim {
    int nr, nc;
    vector<vector<pair<int, uint32_t>>> rows;  // sparse rows, sorted by column after finalize()
    vector<vector<int>> colrows;               // column -> rows that may have an entry there (a superset)
    vector<int> colcnt;                        // column -> number of UNIT entries in live rows
    vector<char> rowdead, coldead;
    bool keep;                                 // record ops and pivot snapshots?
    vector<Op> ops;
    // (pivot column, snapshot of the pivot row) and the row of each pivot, recorded only when keep = true.
    // Read them through snapshots() and pivot_rows(), which stop the program if they were not recorded:
    // using an empty snapshot list silently gives wrong answers (this happened once, in an early spanq).
    vector<pair<int, vector<pair<int, uint32_t>>>> piv_;
    vector<int> pivrow_;
    const vector<pair<int, vector<pair<int, uint32_t>>>> &snapshots() const {
        if (!keep) { fprintf(stderr, "QElim: snapshots requested but keep = false\n"); exit(9); }
        return piv_;
    }
    const vector<int> &pivot_rows() const {
        if (!keep) { fprintf(stderr, "QElim: pivot rows requested but keep = false\n"); exit(9); }
        return pivrow_;
    }
    long long pivots = 0;
    QElim(int r, int c, bool keep_) : nr(r), nc(c), rows(r), colrows(c), colcnt(c, 0), rowdead(r, 0), coldead(c, 0), keep(keep_) {}

    // Append a row (used for extra relation rows); call before finalize().
    void add_row(vector<pair<int, uint32_t>> r) { rows.push_back(std::move(r)); rowdead.push_back(0); nr++; }

    // Sort each row, merge repeated columns, drop zeros, and build the column structures.
    void finalize() {
        for (int r = 0; r < nr; r++) {
            auto &R = rows[r];
            sort(R.begin(), R.end());
            vector<pair<int, uint32_t>> M;
            for (auto &x : R) {
                if (!M.empty() && M.back().first == x.first) M.back().second = (uint32_t)((M.back().second + (u64)x.second) % Q);
                else M.push_back({x.first, (uint32_t)(x.second % Q)});
            }
            R.clear();
            for (auto &x : M) if (x.second) R.push_back(x);
            for (auto &x : R) { colrows[x.first].push_back(r); if (isunit(x.second)) colcnt[x.first]++; }
        }
    }
    // Entry (r, c) of the sparse matrix, found by binary search in row r; 0 if absent.
    u64 get(int r, int c) const {
        auto &R = rows[r];
        auto it = lower_bound(R.begin(), R.end(), make_pair(c, (uint32_t)0));
        return (it != R.end() && it->first == c) ? it->second : 0;
    }
    // Pivot on the unit entry (r, c): clear column c from every other live row.
    void pivot(int r, int c) {
        auto &R = rows[r];
        u64 ui = egcd_inv(get(r, c));
        vector<int> targets;
        for (int r2 : colrows[c]) if (r2 != r && !rowdead[r2]) targets.push_back(r2);
        sort(targets.begin(), targets.end());
        targets.erase(unique(targets.begin(), targets.end()), targets.end());
        for (int r2 : targets) {
            u64 a = get(r2, c);
            if (!a) continue;
            u64 f = mulq(a, ui);                      // row_r2 <- row_r2 - f * row_r
            if (keep) ops.push_back({r2, r, (uint32_t)f});
            auto &S = rows[r2];
            vector<pair<int, uint32_t>> out;
            out.reserve(S.size() + R.size());
            size_t i = 0, j = 0;
            while (i < S.size() || j < R.size()) {    // merge of two sorted sparse rows
                if (j == R.size() || (i < S.size() && S[i].first < R[j].first)) { out.push_back(S[i]); i++; }
                else if (i == S.size() || R[j].first < S[i].first) {
                    u64 v = (Q - mulq(f, R[j].second)) % Q;
                    if (v) { out.push_back({R[j].first, (uint32_t)v}); colrows[R[j].first].push_back(r2); if (isunit(v)) colcnt[R[j].first]++; }
                    j++;
                } else {
                    u64 old = S[i].second, v = (old + Q - mulq(f, R[j].second)) % Q;
                    if (isunit(old)) colcnt[S[i].first]--;
                    if (v) { out.push_back({S[i].first, (uint32_t)v}); if (isunit(v)) colcnt[S[i].first]++; }
                    i++; j++;
                }
            }
            S.swap(out);
        }
        for (auto &x : R) if (isunit(x.second)) colcnt[x.first]--;
        rowdead[r] = 1; coldead[c] = 1; colcnt[c] = 0; pivots++;
        if (keep) { piv_.push_back({c, R}); pivrow_.push_back(r); }
        vector<pair<int, uint32_t>>().swap(rows[r]);
        vector<int>().swap(colrows[c]);
    }
    // Markowitz-style order: repeatedly take the column with the fewest unit entries and, in it, the shortest
    // row whose entry is a unit.  Stops when no live column has a unit entry.
    void run() {
        bool any = true;
        while (any) {
            any = false;
            typedef pair<int, int> PI;
            priority_queue<PI, vector<PI>, greater<PI>> pq;
            for (int c = 0; c < nc; c++) if (!coldead[c] && colcnt[c] > 0) pq.push({colcnt[c], c});
            while (!pq.empty()) {
                auto [cnt, c] = pq.top(); pq.pop();
                if (coldead[c] || colcnt[c] != cnt || cnt == 0) continue;   // stale queue entry
                if (colrows[c].size() > 2 * (size_t)colcnt[c] + 8) {        // compact the row list
                    vector<int> M;
                    for (int r : colrows[c]) if (!rowdead[r] && get(r, c)) M.push_back(r);
                    sort(M.begin(), M.end()); M.erase(unique(M.begin(), M.end()), M.end());
                    colrows[c].swap(M);
                }
                int best = -1; size_t bl = SIZE_MAX;
                for (int r : colrows[c]) if (!rowdead[r] && rows[r].size() < bl && isunit(get(r, c))) { bl = rows[r].size(); best = r; }
                if (best < 0) continue;
                vector<int> touched;
                for (auto &x : rows[best]) touched.push_back(x.first);
                pivot(best, c);
                any = true;
                for (int c2 : touched) if (!coldead[c2] && colcnt[c2] > 0) pq.push({colcnt[c2], c2});
            }
        }
    }
    // Rows still alive and nonzero, and the columns they touch: the residual block.
    void residual(vector<int> &rr, vector<int> &rc, vector<int> &cpos) const {
        rr.clear(); rc.clear(); cpos.assign(nc, -1);
        for (int r = 0; r < nr; r++)
            if (!rowdead[r] && !rows[r].empty()) {
                rr.push_back(r);
                for (auto &x : rows[r]) if (cpos[x.first] < 0) { cpos[x.first] = 0; rc.push_back(x.first); }
            }
        sort(rc.begin(), rc.end());
        for (size_t j = 0; j < rc.size(); j++) cpos[rc[j]] = j;
    }
};

// Dense local Smith form mod Q:  U R V = diag(p^dv[0], p^dv[1], ...) with dv nondecreasing.
// R is nr x nc and is overwritten.  U (nr x nr) is computed only if wantU.  At each step the entry of least
// valuation in the remaining block is moved to the corner (rows and columns swapped, also in U and V), the
// column below it is cleared by row operations and the row to its right by column operations.  The quotient
// q = (R_is / p^v) * unit^{-1} is legitimate because every remaining entry has valuation >= v.  Stops when the
// remaining block is zero mod Q.  Returns the valuations dv.
vector<int> local_smith(vector<vector<u64>> &R, int nr, int nc, vector<vector<u64>> &U, vector<vector<u64>> &Vm, bool wantU, int KK) {
    if (wantU) { U.assign(nr, vector<u64>(nr, 0)); for (int i = 0; i < nr; i++) U[i][i] = 1; }
    Vm.assign(nc, vector<u64>(nc, 0));
    for (int j = 0; j < nc; j++) Vm[j][j] = 1;
    vector<int> dv;
    int r = min(nr, nc);
    for (int s = 0; s < r; s++) {
        int bi = -1, bj = -1, bv = 1000;
        for (int i = s; i < nr && bv > 0; i++)
            for (int j = s; j < nc; j++) { int v = val(R[i][j]); if (v < bv) { bv = v; bi = i; bj = j; if (!v) break; } }
        if (bv >= KK) break;
        swap(R[s], R[bi]);
        if (wantU) swap(U[s], U[bi]);
        if (bj != s) { for (int i = 0; i < nr; i++) swap(R[i][s], R[i][bj]); for (int i = 0; i < nc; i++) swap(Vm[i][s], Vm[i][bj]); }
        u64 pv = ppow(bv), ui = egcd_inv(R[s][s] / pv);
        for (int i = s + 1; i < nr; i++) if (R[i][s] % Q) {               // clear the column below the corner
            u64 q = mulq(R[i][s] / pv, ui);
            for (int j = s; j < nc; j++) if (R[s][j]) R[i][j] = (R[i][j] + Q - mulq(q, R[s][j])) % Q;
            if (wantU) for (int j = 0; j < nr; j++) if (U[s][j]) U[i][j] = (U[i][j] + Q - mulq(q, U[s][j])) % Q;
        }
        for (int j = s + 1; j < nc; j++) if (R[s][j] % Q) {               // clear the row to its right
            u64 q = mulq(R[s][j] / pv, ui);
            R[s][j] = 0;                                                  // column s is zero below s now
            for (int i = 0; i < nc; i++) if (Vm[i][s]) Vm[i][j] = (Vm[i][j] + Q - mulq(q, Vm[i][s])) % Q;
        }
        dv.push_back(bv);
    }
    return dv;
}
