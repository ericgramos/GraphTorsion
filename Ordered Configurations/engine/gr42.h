// =====================================================================================================
//  gr42.h  --  elimination over the Galois ring GR(4,2) = Z/4[w] / (w^2 + w + 1)
//
//  WHY THIS EXISTS
//  ---------------
//  Suppose a group G of ODD order acts on the cell complex C, and G has a quotient Z/3 (for example a
//  rotation of a triangle).  The characters of Z/3 take values in the cube roots of unity, so the parts of
//  H_*(C) on which the rotation acts nontrivially are invisible in the plain coinvariants C_G.  They are
//  seen by the TWISTED coinvariants  C (x)_{Z[G]} R_chi  with  R = Z_(2)[w],  w^2 + w + 1 = 0, where G acts
//  on R through chi(g) = w^{t(g)}.  The prime 2 stays prime in R (2 is inert in Q(w)), R/2 = F_4, and
//  R/4 = GR(4,2) is the Galois ring with 16 elements (Wan, "Lectures on finite fields and Galois rings").
//  Computing the local Smith form over R/4 gives the number of R/2 summands (each one is (Z/2)^2 as an
//  abelian group) in the omega-isotypic part of H_*(Conf_n; Z_(2)).  The trivial part comes from an
//  ordinary run without -W; trivial part + omega part together recover the coinvariants for the
//  subgroup of G that acts trivially on the rotation (see README, Section 5).
//
//  ENCODING
//  --------
//  An element a + b w (a, b in Z/4) is coded as the integer  a + 4 b  in 0..15.  After the first stage
//  everything is divisible by 2 and is halved; the ring becomes R/2 = F_4, with the same encoding but
//  a, b in {0, 1}.  "Level" L = 0 means GR(4,2), L = 1 means F_4.  An element is a unit iff it is nonzero
//  modulo 2, i.e. iff (a mod 2, b mod 2) != (0, 0), i.e. iff (code & 5) != 0.
//
//  Multiplication uses w^2 = -1 - w:
//     (a + b w)(c + d w) = (a c - b d) + (a d + b c - b d) w.
// =====================================================================================================
#pragma once

struct GRT {
    uint8_t add[2][16][16], mul[2][16][16], neg[2][16], inv[2][16];
    GRT() {
        for (int L = 0; L < 2; L++) {
            int M = L == 0 ? 4 : 2;
            for (int x = 0; x < 16; x++)
                for (int y = 0; y < 16; y++) {
                    int a = x & 3, b = x >> 2, c = y & 3, d = y >> 2;
                    if (a >= M || b >= M || c >= M || d >= M) {
                        add[L][x][y] = mul[L][x][y] = 0;
                        continue;
                    }
                    add[L][x][y] = ((a + c) % M) | (((b + d) % M) << 2);
                    int r = ((a * c - b * d) % M + M) % M, s = ((a * d + b * c - b * d) % M + M) % M;
                    mul[L][x][y] = r | (s << 2);
                }
            for (int x = 0; x < 16; x++) {
                int a = x & 3, b = x >> 2;
                neg[L][x] = ((M - a) % M) | (((M - b) % M) << 2);
                inv[L][x] = 0;
                for (int y = 0; y < 16; y++)
                    if (mul[L][x][y] == 1) {
                        inv[L][x] = y;
                        break;
                    }
            }
        }
    }
    static bool unit(int x) {
        return (x & 5) != 0;
    }
};
static GRT grt;

// The coefficient (+-1) * w^t as a GR(4,2) code.  w^0 = 1 -> 1; w^1 -> code 4; w^2 = -1 - w -> 3 + 3*4.
static inline int grcode(int sign, int t) {
    static const int wp[3] = {1, 4, 3 | (3 << 2)};
    int c = wp[((t % 3) + 3) % 3];
    return sign > 0 ? c : grt.neg[0][c];
}

// -----------------------------------------------------------------------------------------------------
//  LGR: sparse local elimination over GR(4,2), then F_4.  Same structure as LPE in oconf.cpp; entries are
//  packed as (column << 5) | code.  Only units are used as pivots; between the two stages every entry is
//  halved.  A dense phase (dense_gr below) takes over when at most denseMax columns are alive.
// -----------------------------------------------------------------------------------------------------
struct LGR {
    int nr, nc, L = 0;
    vector<vector<uint32_t>> rows;
    vector<vector<int>> colrows;
    vector<int> colcnt;
    vector<char> rowdead, coldead;
    long long pivots = 0, lastCheck = 0;
    bool goDense = false;
    ll maxabs = 1;
    vector<long long> stages;
    long long leftover_nonzero = 0;
    static inline long long denseMax = 70000, denseMin = 30000, nnzBudget = 120000000;
    static inline int cntThresh = 48;
    LGR(int r, int c) : nr(r), nc(c), rows(r), colrows(c), colcnt(c, 0), rowdead(r, 0), coldead(c, 0) {}
    static uint32_t col(uint32_t x) {
        return x >> 5;
    }
    static uint32_t val(uint32_t x) {
        return x & 31;
    }
    void finalize() {
        for (int r = 0; r < nr; r++) {
            auto &R = rows[r];
            sort(R.begin(), R.end());
            vector<uint32_t> M;
            for (size_t i = 0; i < R.size();) {
                uint32_t c = col(R[i]);
                int v = 0;
                while (i < R.size() && col(R[i]) == c) {
                    v = grt.add[L][v][val(R[i])];
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
    int get(int r, int c) const {
        auto &R = rows[r];
        auto it = lower_bound(R.begin(), R.end(), (uint32_t)c << 5);
        return (it != R.end() && col(*it) == (uint32_t)c) ? (int)val(*it) : 0;
    }
    void pivot(int r, int c) {
        auto &R = rows[r];
        int ui = grt.inv[L][get(r, c)];
        vector<int> tg;
        for (int r2 : colrows[c])
            if (r2 != r && !rowdead[r2])
                tg.push_back(r2);
        sort(tg.begin(), tg.end());
        tg.erase(unique(tg.begin(), tg.end()), tg.end());
        for (int r2 : tg) {
            int a = get(r2, c);
            if (!a)
                continue;
            int f = grt.neg[L][grt.mul[L][a][ui]];   // row2 <- row2 + f * row, f = -a/u
            auto &S = rows[r2];
            vector<uint32_t> out;
            out.reserve(S.size() + R.size());
            size_t i = 0, j = 0;
            while (i < S.size() || j < R.size()) {
                if (j == R.size() || (i < S.size() && col(S[i]) < col(R[j]))) {
                    out.push_back(S[i]);
                    i++;
                } else if (i == S.size() || col(R[j]) < col(S[i])) {
                    int v = grt.mul[L][f][val(R[j])];
                    if (v) {
                        out.push_back(col(R[j]) << 5 | v);
                        colrows[col(R[j])].push_back(r2);
                        colcnt[col(R[j])]++;
                    }
                    j++;
                } else {
                    int v = grt.add[L][val(S[i])][grt.mul[L][f][val(R[j])]];
                    if (v)
                        out.push_back(col(S[i]) << 5 | v);
                    else
                        colcnt[col(S[i])]--;
                    i++;
                    j++;
                }
            }
            S.assign(out.begin(), out.end());
        }
        for (auto x : R)
            colcnt[col(x)]--;
        rowdead[r] = 1;
        coldead[c] = 1;
        colcnt[c] = 0;
        pivots++;
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
            if (pq.size() > 4 * (size_t)nc + 1000000) {   // drop stale queue entries
                vector<PI> keep;
                keep.reserve(nc);
                while (!pq.empty()) {
                    auto t = pq.top();
                    pq.pop();
                    if (!coldead[t.second] && colcnt[t.second] == t.first && t.first > 0)
                        keep.push_back(t);
                }
                sort(keep.begin(), keep.end());
                keep.erase(unique(keep.begin(), keep.end()), keep.end());
                pq = priority_queue<PI, vector<PI>, greater<PI>>(keep.begin(), keep.end());
            }
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
                    int v = get(r, c);
                    if (v && GRT::unit(v)) {
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
            if (L == 0 && denseMax > 0 && cnt >= cntThresh && pivots - lastCheck >= 2000) {
                lastCheck = pivots;
                long long lc = 0;
                for (int c2 = 0; c2 < nc; c2++)
                    if (!coldead[c2] && colcnt[c2] > 0)
                        lc++;
                long long nz = 0;
                if (lc <= denseMax && lc > denseMin)
                    for (int r2 = 0; r2 < nr; r2++)
                        if (!rowdead[r2])
                            nz += rows[r2].size();
                if (lc <= denseMax && (lc <= denseMin || nz > nnzBudget)) {
                    goDense = true;
                    fprintf(stderr, "  [gr dense] switching: pivots=%lld liveCols=%lld colcnt=%d\n", pivots, lc, cnt);
                    return true;
                }
            }
            for (int c2 : touched)
                if (!coldead[c2] && colcnt[c2] > 0)
                    pq.push({colcnt[c2], c2});
        }
        return any;
    }

    // ---------------------------------------------------------------------------------------------------
    //  Dense phase over GR(4,2), bit-sliced.  Planes: 0 = a0, 1 = a1 (the Z/4 "real" part a = a0 + 2 a1),
    //  2 = b0, 3 = b1 (the coefficient b of w).  Same streaming reduced-echelon algorithm as LPE's dense
    //  phases; leftovers (all entries even) are halved to F_4 vectors (planes a1, b1) and kept as an
    //  online F_4 row-echelon basis.
    // ---------------------------------------------------------------------------------------------------
    struct R4 {
        uint64_t *p[4];
    };
    static inline void z4add(uint64_t &x0, uint64_t &x1, uint64_t y0, uint64_t y1) {
        uint64_t c = x0 & y0;
        x0 ^= y0;
        x1 ^= y1 ^ c;
    }
    static void gradd(uint64_t *const *X, uint64_t *const *Y, int W) {   // X += Y componentwise in Z/4
        uint64_t *a0 = X[0], *a1 = X[1], *b0 = X[2], *b1 = X[3];
        const uint64_t *c0 = Y[0], *c1 = Y[1], *d0 = Y[2], *d1 = Y[3];
        for (int w = 0; w < W; w++) {
            z4add(a0[w], a1[w], c0[w], c1[w]);
            z4add(b0[w], b1[w], d0[w], d1[w]);
        }
    }
    static inline int grget(uint64_t *const *X, int c) {
        int w = c >> 6, s = c & 63;
        return (int)((X[0][w] >> s) & 1) | (int)(((X[1][w] >> s) & 1) << 1) | (int)(((X[2][w] >> s) & 1) << 2) |
               (int)(((X[3][w] >> s) & 1) << 3);
    }
    static inline void grset(uint64_t *const *X, int c, int v) {
        int w = c >> 6;
        uint64_t bit = 1ULL << (c & 63);
        for (int i = 0; i < 4; i++) {
            if ((v >> i) & 1)
                X[i][w] |= bit;
            else
                X[i][w] &= ~bit;
        }
    }
    // O = f * Y for a GR code f = c + d w:  f*y = c*y + d*(w*y),  w*(a + b w) = -b + (a - b) w.
    static void grscale(uint64_t *const *O, uint64_t *const *Y, int f, int W) {
        int c = f & 3, d = f >> 2;
        for (int w = 0; w < W; w++) {
            uint64_t a0 = Y[0][w], a1 = Y[1][w], b0 = Y[2][w], b1 = Y[3][w];
            uint64_t nb0 = b0, nb1 = b1 ^ b0;      // -b
            uint64_t wa0 = nb0, wa1 = nb1;         // real part of w*y
            uint64_t wb0 = a0, wb1 = a1;
            z4add(wb0, wb1, nb0, nb1);             // w-part of w*y: a - b
            uint64_t r0 = 0, r1 = 0, s0 = 0, s1 = 0;
            auto addc = [&](int k, uint64_t x0, uint64_t x1, uint64_t &o0, uint64_t &o1) {   // o += k*x in Z/4
                if (k == 1)
                    z4add(o0, o1, x0, x1);
                else if (k == 2)
                    o1 ^= x0;
                else if (k == 3)
                    z4add(o0, o1, x0, x1 ^ x0);
            };
            addc(c, a0, a1, r0, r1);
            addc(c, b0, b1, s0, s1);
            addc(d, wa0, wa1, r0, r1);
            addc(d, wb0, wb1, s0, s1);
            O[0][w] = r0;
            O[1][w] = r1;
            O[2][w] = s0;
            O[3][w] = s1;
        }
    }
    // F_4 helpers on two planes (a, b) meaning a + b w:  (A, B) += (c + d w) * (X, Y).
    static void f4axpy(vector<uint64_t> &A, vector<uint64_t> &Bp, const uint64_t *X, const uint64_t *Y, int c, int d,
                       int W) {
        for (int w = 0; w < W; w++) {
            uint64_t x = X[w], y = Y[w];
            uint64_t na = (c ? x : 0) ^ (d ? y : 0), nb = (c ? y : 0) ^ (d ? x : 0) ^ (d ? y : 0);
            A[w] ^= na;
            Bp[w] ^= nb;
        }
    }
    void dense_gr(long long sparseStagePivots) {
        vector<int> cmap(nc, -1);
        int m = 0;
        for (int c = 0; c < nc; c++)
            if (!coldead[c] && colcnt[c] > 0)
                cmap[c] = m++;
        vector<vector<int>>().swap(colrows);
        int W = (m + 63) / 64;
        fprintf(stderr, "  [gr] m=%d, basis memory <= %.2f GB\n", m, 4.0 * m * W * 8 / 1e9);
        vector<R4> B;
        vector<int> rowOfPiv(m, -1);
        vector<uint64_t> xs(4 * (size_t)W);
        uint64_t *X[4] = {xs.data(), xs.data() + W, xs.data() + 2 * W, xs.data() + 3 * W};
        vector<uint64_t> accs(64 * (size_t)W), tmps(4 * (size_t)W);
        uint64_t *T[4] = {tmps.data(), tmps.data() + W, tmps.data() + 2 * W, tmps.data() + 3 * W};
        vector<uint64_t> Ms(64 * (size_t)W);   // the 16 multiples of the newest pivot row
        auto MP = [&](int f, int i) {
            return Ms.data() + (4 * f + i) * (size_t)W;
        };
        auto AP = [&](int f, int i) {
            return accs.data() + (4 * f + i) * (size_t)W;
        };
        vector<vector<uint64_t>> Ya, Yb;       // online F_4 echelon basis of halved leftovers
        vector<int> ypiv;
        auto f4get = [&](const vector<uint64_t> &a, const vector<uint64_t> &b, int c) {
            int w = c >> 6, s = c & 63;
            return (int)((a[w] >> s) & 1) | (int)(((b[w] >> s) & 1) << 1);
        };
        auto f4inv = [](int v) {   // codes: 1 = 1, 2 = w, 3 = 1 + w = w^2;  w^{-1} = w^2
            return v == 1 ? 1 : (v == 2 ? 3 : 2);
        };
        auto insertY = [&](vector<uint64_t> &a, vector<uint64_t> &b) {
            for (size_t i = 0; i < Ya.size(); i++) {
                int v = f4get(a, b, ypiv[i]);
                if (v)
                    f4axpy(a, b, Ya[i].data(), Yb[i].data(), v & 1, v >> 1, W);   // char 2: subtract = add
            }
            int c = -1;
            for (int w = 0; w < W && c < 0; w++) {
                uint64_t q = a[w] | b[w];
                if (q)
                    c = w * 64 + __builtin_ctzll(q);
            }
            if (c < 0)
                return;
            int iv = f4inv(f4get(a, b, c));
            if (iv != 1) {
                vector<uint64_t> na(W, 0), nb(W, 0);
                f4axpy(na, nb, a.data(), b.data(), iv & 1, iv >> 1, W);
                a.swap(na);
                b.swap(nb);
            }
            for (size_t i = 0; i < Ya.size(); i++) {
                int v = f4get(Ya[i], Yb[i], c);
                if (v)
                    f4axpy(Ya[i], Yb[i], a.data(), b.data(), v & 1, v >> 1, W);
            }
            Ya.push_back(a);
            Yb.push_back(b);
            ypiv.push_back(c);
        };
        long long done = 0, stagePiv = sparseStagePivots;
        for (int r = 0; r < nr; r++) {
            if (rowdead[r] || rows[r].empty())
                continue;
            vector<pair<int, int>> ents;
            for (auto e : rows[r]) {
                int c = cmap[col(e)];
                if (c < 0) {
                    fprintf(stderr, "gr: dead col\n");
                    exit(6);
                }
                ents.push_back({c, (int)val(e)});
            }
            vector<uint32_t>().swap(rows[r]);
            fill(xs.begin(), xs.end(), 0);
            for (auto &[c, v] : ents)
                grset(X, c, v);
            bool used[16] = {false};   // group basis rows by coefficient g = -x[c], then add g * (sum)
            for (auto &[c, v] : ents) {
                int br = rowOfPiv[c];
                if (br < 0)
                    continue;
                int f = grget(X, c);
                if (!f)
                    continue;
                int g = grt.neg[0][f];
                uint64_t *A[4] = {AP(g, 0), AP(g, 1), AP(g, 2), AP(g, 3)};
                if (!used[g]) {
                    for (int i = 0; i < 4; i++)
                        memcpy(A[i], B[br].p[i], W * 8);
                    used[g] = true;
                } else
                    gradd(A, B[br].p, W);
            }
            for (int g = 1; g < 16; g++)
                if (used[g]) {
                    uint64_t *A[4] = {AP(g, 0), AP(g, 1), AP(g, 2), AP(g, 3)};
                    grscale(T, A, g, W);
                    gradd(X, T, W);
                }
            int j = -1;                // first unit: a0 or b0 bit set
            for (int w = 0; w < W && j < 0; w++) {
                uint64_t q = X[0][w] | X[2][w];
                if (q)
                    j = w * 64 + __builtin_ctzll(q);
            }
            if (j >= 0) {
                int u = grt.inv[0][grget(X, j)];
                if (u != 1) {
                    grscale(T, X, u, W);
                    for (int i = 0; i < 4; i++)
                        memcpy(X[i], T[i], W * 8);
                }
                for (int f = 0; f < 16; f++) {
                    uint64_t *O[4] = {MP(f, 0), MP(f, 1), MP(f, 2), MP(f, 3)};
                    grscale(O, X, f, W);
                }
                long long nb = B.size();
#pragma omp parallel for schedule(static, 64)
                for (long long b = 0; b < nb; b++) {
                    int f = grget(B[b].p, j);
                    if (f) {
                        int g = grt.neg[0][f];
                        uint64_t *Y[4] = {MP(g, 0), MP(g, 1), MP(g, 2), MP(g, 3)};
                        gradd(B[b].p, Y, W);
                    }
                }
                R4 nr4;
                for (int i = 0; i < 4; i++) {
                    nr4.p[i] = (uint64_t *)malloc(W * 8);
                    if (!nr4.p[i]) {
                        fprintf(stderr, "gr: out of memory\n");
                        exit(9);
                    }
                    memcpy(nr4.p[i], X[i], W * 8);
                }
                rowOfPiv[j] = B.size();
                B.push_back(nr4);
                stagePiv++;
                pivots++;
                rowdead[r] = 1;
            } else {
                bool any = false;
                for (int w = 0; w < W; w++)
                    if (X[1][w] | X[3][w]) {
                        any = true;
                        break;
                    }
                if (any) {
                    vector<uint64_t> a(X[1], X[1] + W), b(X[3], X[3] + W);   // halve: (a1, b1)
                    insertY(a, b);
                }
            }
            if (++done % 200000 == 0)
                fprintf(stderr, "  [gr] %lld rows, basis %zu, leftover span %zu\n", done, B.size(), Ya.size());
        }
        stages.push_back(stagePiv);
        fprintf(stderr, "  [gr] stage 0 done: basis %zu, leftover span %zu\n", B.size(), Ya.size());
        vector<vector<uint64_t>> F0, F1;
        for (size_t i = 0; i < Ya.size(); i++) {
            vector<uint64_t> a = Ya[i], b = Yb[i];
            vector<pair<int, int>> hits;
            for (int w = 0; w < W; w++) {
                uint64_t q = a[w] | b[w];
                while (q) {
                    int c = w * 64 + __builtin_ctzll(q);
                    q &= q - 1;
                    if (rowOfPiv[c] >= 0)
                        hits.push_back({c, f4get(a, b, c)});
                }
            }
            for (auto &[c, v] : hits) {   // basis row mod 2 = planes (a0, b0)
                auto &bb = B[rowOfPiv[c]];
                f4axpy(a, b, bb.p[0], bb.p[2], v & 1, v >> 1, W);
            }
            bool nz = false;
            for (int w = 0; w < W; w++) {
                uint64_t q = a[w] | b[w];
                while (q) {
                    int c = w * 64 + __builtin_ctzll(q);
                    q &= q - 1;
                    if (rowOfPiv[c] >= 0) {
                        fprintf(stderr, "gr: pivot col nonzero\n");
                        exit(7);
                    }
                    nz = true;
                }
            }
            if (nz) {
                F0.push_back(move(a));
                F1.push_back(move(b));
            }
        }
        for (auto &r4 : B)
            for (int i = 0; i < 4; i++)
                free(r4.p[i]);
        long long nrw = F0.size(), r1 = 0;   // F_4 rank of the final leftovers = number of R/2 summands
        for (int c = 0; c < m && r1 < nrw; c++) {
            long long piv = -1;
            for (long long i = r1; i < nrw; i++)
                if (f4get(F0[i], F1[i], c)) {
                    piv = i;
                    break;
                }
            if (piv < 0)
                continue;
            swap(F0[piv], F0[r1]);
            swap(F1[piv], F1[r1]);
            int iv = f4inv(f4get(F0[r1], F1[r1], c));
            if (iv != 1) {
                vector<uint64_t> na(W, 0), nb(W, 0);
                f4axpy(na, nb, F0[r1].data(), F1[r1].data(), iv & 1, iv >> 1, W);
                F0[r1].swap(na);
                F1[r1].swap(nb);
            }
            for (long long i = 0; i < nrw; i++)
                if (i != r1) {
                    int v = f4get(F0[i], F1[i], c);
                    if (v)
                        f4axpy(F0[i], F1[i], F0[r1].data(), F1[r1].data(), v & 1, v >> 1, W);
                }
            r1++;
        }
        if (r1) {
            stages.push_back(r1);
            pivots += r1;
        }
        fprintf(stderr, "  [gr] stage 1 (R/2 summands): %lld\n", r1);
    }

    void run() {
        long long before = 0;
        while (true) {
            while (sweep()) {
                if (goDense)
                    break;
            }
            if (goDense) {
                dense_gr(pivots - before);
                return;
            }
            stages.push_back(pivots - before);
            before = pivots;
            bool live = false;
            for (int r = 0; r < nr && !live; r++)
                if (!rowdead[r] && !rows[r].empty())
                    live = true;
            if (!live)
                break;
            if (L == 1) {
                for (int r = 0; r < nr; r++)
                    if (!rowdead[r])
                        leftover_nonzero += rows[r].size();
                break;
            }
            for (int r = 0; r < nr; r++)   // halve everything: GR(4,2) -> F_4
                if (!rowdead[r]) {
                    vector<uint32_t> M;
                    for (auto x : rows[r]) {
                        int v = val(x);
                        if (GRT::unit(v)) {
                            fprintf(stderr, "gr: nonunit not divisible\n");
                            exit(4);
                        }
                        int h = ((v & 3) >> 1) | (((v >> 2) >> 1) << 2);
                        if (h)
                            M.push_back(col(x) << 5 | h);
                        else
                            colcnt[col(x)]--;
                    }
                    rows[r].swap(M);
                }
            L = 1;
        }
    }
};
