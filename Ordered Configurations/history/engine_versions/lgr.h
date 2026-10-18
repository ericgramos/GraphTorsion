// ---------- GR(4,2) = Z/4[w]/(w^2+w+1) elimination (the omega-isotypic part of a Z/3-symmetric complex, p = 2) ----------
// element code: a + 4*b  (a, b in Z/4 at level 4; in {0,1} at level 2, i.e. F4)
struct GRT {
    uint8_t add[2][16][16], mul[2][16][16], neg[2][16], inv[2][16];
    GRT() {
        for (int L = 0; L < 2; L++) { int M = L == 0 ? 4 : 2;
            for (int x = 0; x < 16; x++) for (int y = 0; y < 16; y++) {
                int a = x & 3, b = x >> 2, c = y & 3, d = y >> 2;
                if (a >= M || b >= M || c >= M || d >= M) { add[L][x][y] = mul[L][x][y] = 0; continue; }
                add[L][x][y] = ((a + c) % M) | (((b + d) % M) << 2);
                int r = ((a * c - b * d) % M + M) % M, s = ((a * d + b * c - b * d) % M + M) % M;   // w^2 = -1 - w
                mul[L][x][y] = r | (s << 2);
            }
            for (int x = 0; x < 16; x++) { int a = x & 3, b = x >> 2; neg[L][x] = ((M - a) % M) | (((M - b) % M) << 2); inv[L][x] = 0;
                for (int y = 0; y < 16; y++) if (mul[L][x][y] == 1) { inv[L][x] = y; break; } }
        }
    }
    static bool unit(int x) { return (x & 5) != 0; }
};
static GRT grt;
// omega^t * (+-1) as a code
static inline int grcode(int sign, int t) { static const int wp[3] = {1, 4, 3 | (3 << 2)}; int c = wp[((t % 3) + 3) % 3]; return sign > 0 ? c : grt.neg[0][c]; }

struct LGR {
    int nr, nc, L = 0;   // L = 0: GR(4,2); L = 1: F4
    vector<vector<uint32_t>> rows; vector<vector<int>> colrows; vector<int> colcnt; vector<char> rowdead, coldead;
    long long pivots = 0, lastCheck = 0; bool goDense = false; ll maxabs = 1; vector<long long> stages; long long leftover_nonzero = 0;
    static inline long long denseMax = 70000, denseMin = 30000, nnzBudget = 120000000; static inline int cntThresh = 48;
    LGR(int r, int c) : nr(r), nc(c), rows(r), colrows(c), colcnt(c, 0), rowdead(r, 0), coldead(c, 0) {}
    static uint32_t col(uint32_t x) { return x >> 5; } static uint32_t val(uint32_t x) { return x & 31; }
    void finalize() {
        for (int r = 0; r < nr; r++) {
            auto& R = rows[r]; sort(R.begin(), R.end());
            vector<uint32_t> M;
            for (size_t i = 0; i < R.size();) { uint32_t c = col(R[i]); int v = 0; while (i < R.size() && col(R[i]) == c) { v = grt.add[L][v][val(R[i])]; i++; } if (v) M.push_back(c << 5 | v); }
            R.swap(M); for (auto x : R) { colrows[col(x)].push_back(r); colcnt[col(x)]++; }
        }
    }
    int get(int r, int c) const { auto& R = rows[r]; auto it = lower_bound(R.begin(), R.end(), (uint32_t)c << 5); return (it != R.end() && col(*it) == (uint32_t)c) ? (int)val(*it) : 0; }
    void pivot(int r, int c) {
        auto& R = rows[r]; int ui = grt.inv[L][get(r, c)];
        vector<int> tg; for (int r2 : colrows[c]) if (r2 != r && !rowdead[r2]) tg.push_back(r2);
        sort(tg.begin(), tg.end()); tg.erase(unique(tg.begin(), tg.end()), tg.end());
        for (int r2 : tg) {
            int a = get(r2, c); if (!a) continue; int f = grt.neg[L][grt.mul[L][a][ui]];   // row2 += f * row
            auto& S = rows[r2]; vector<uint32_t> out; out.reserve(S.size() + R.size()); size_t i = 0, j = 0;
            while (i < S.size() || j < R.size()) {
                if (j == R.size() || (i < S.size() && col(S[i]) < col(R[j]))) { out.push_back(S[i]); i++; }
                else if (i == S.size() || col(R[j]) < col(S[i])) { int v = grt.mul[L][f][val(R[j])]; if (v) { out.push_back(col(R[j]) << 5 | v); colrows[col(R[j])].push_back(r2); colcnt[col(R[j])]++; } j++; }
                else { int v = grt.add[L][val(S[i])][grt.mul[L][f][val(R[j])]]; if (v) out.push_back(col(S[i]) << 5 | v); else colcnt[col(S[i])]--; i++; j++; }
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
        while (!pq.empty()) { auto [cnt, c] = pq.top(); pq.pop();
            if (pq.size() > 4 * (size_t)nc + 1000000) {
                vector<PI> keep; keep.reserve(nc); while (!pq.empty()) { auto t = pq.top(); pq.pop(); if (!coldead[t.second] && colcnt[t.second] == t.first && t.first > 0) keep.push_back(t); }
                sort(keep.begin(), keep.end()); keep.erase(unique(keep.begin(), keep.end()), keep.end());
                pq = priority_queue<PI, vector<PI>, greater<PI>>(keep.begin(), keep.end()); }
            if (coldead[c] || colcnt[c] != cnt || cnt == 0) continue;
            if (colrows[c].size() > 2 * (size_t)colcnt[c] + 8) { vector<int> M; for (int r : colrows[c]) if (!rowdead[r] && get(r, c)) M.push_back(r); sort(M.begin(), M.end()); M.erase(unique(M.begin(), M.end()), M.end()); colrows[c].swap(M); }
            int best = -1; size_t bl = SIZE_MAX;
            for (int r : colrows[c]) if (!rowdead[r] && rows[r].size() < bl) { int v = get(r, c); if (v && GRT::unit(v)) { bl = rows[r].size(); best = r; } }
            if (best < 0) continue;
            vector<int> touched; for (auto x : rows[best]) touched.push_back(col(x));
            pivot(best, c); any = true;
            if (L == 0 && denseMax > 0 && cnt >= cntThresh && pivots - lastCheck >= 2000) { lastCheck = pivots;
                long long lc = 0; for (int c2 = 0; c2 < nc; c2++) if (!coldead[c2] && colcnt[c2] > 0) lc++;
                long long nz = 0; if (lc <= denseMax && lc > denseMin) for (int r2 = 0; r2 < nr; r2++) if (!rowdead[r2]) nz += rows[r2].size();
                if (lc <= denseMax && (lc <= denseMin || nz > nnzBudget)) { goDense = true; fprintf(stderr, "  [gr dense] switching: pivots=%lld liveCols=%lld colcnt=%d\n", pivots, lc, cnt); return true; } }
            for (int c2 : touched) if (!coldead[c2] && colcnt[c2] > 0) pq.push({colcnt[c2], c2});
        }
        return any;
    }
#include "grdense.h"
    void run() {
        long long before = 0;
        while (true) {
            while (sweep()) { if (goDense) break; }
            if (goDense) { dense_gr(pivots - before); return; }
            stages.push_back(pivots - before); before = pivots;
            bool live = false; for (int r = 0; r < nr && !live; r++) if (!rowdead[r] && !rows[r].empty()) live = true;
            if (!live) break;
            if (L == 1) { for (int r = 0; r < nr; r++) if (!rowdead[r]) leftover_nonzero += rows[r].size(); break; }
            for (int r = 0; r < nr; r++) if (!rowdead[r]) {
                vector<uint32_t> M; for (auto x : rows[r]) { int v = val(x); if (GRT::unit(v)) { fprintf(stderr, "gr: nonunit not divisible\n"); exit(4); } int h = ((v & 3) >> 1) | (((v >> 2) >> 1) << 2); if (h) M.push_back(col(x) << 5 | h); else colcnt[col(x)]--; }
                rows[r].swap(M);
            }
            L = 1;
        }
    }
};
