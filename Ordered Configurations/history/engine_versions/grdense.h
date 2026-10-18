    // ---- dense phase over GR(4,2), bit-sliced planes (a0, a1, b0, b1) ----
    struct R4 { uint64_t* p[4]; };
    static inline void z4add(uint64_t& x0, uint64_t& x1, uint64_t y0, uint64_t y1) { uint64_t c = x0 & y0; x0 ^= y0; x1 ^= y1 ^ c; }
    static void gradd(uint64_t* const* X, uint64_t* const* Y, int W) {
        uint64_t *a0 = X[0], *a1 = X[1], *b0 = X[2], *b1 = X[3]; const uint64_t *c0 = Y[0], *c1 = Y[1], *d0 = Y[2], *d1 = Y[3];
        for (int w = 0; w < W; w++) { z4add(a0[w], a1[w], c0[w], c1[w]); z4add(b0[w], b1[w], d0[w], d1[w]); }
    }
    static inline int grget(uint64_t* const* X, int c) { int w = c >> 6, s = c & 63; return (int)((X[0][w] >> s) & 1) | (int)(((X[1][w] >> s) & 1) << 1) | (int)(((X[2][w] >> s) & 1) << 2) | (int)(((X[3][w] >> s) & 1) << 3); }
    static inline void grset(uint64_t* const* X, int c, int v) { int w = c >> 6; uint64_t bit = 1ULL << (c & 63); for (int i = 0; i < 4; i++) { if ((v >> i) & 1) X[i][w] |= bit; else X[i][w] &= ~bit; } }
    // out = f * y  (f a GR code), all rows W words, 4 planes each
    static void grscale(uint64_t* const* O, uint64_t* const* Y, int f, int W) {
        int c = f & 3, d = f >> 2;
        for (int w = 0; w < W; w++) {
            uint64_t a0 = Y[0][w], a1 = Y[1][w], b0 = Y[2][w], b1 = Y[3][w];
            // w*y = (-b, a - b)
            uint64_t nb0 = b0, nb1 = b1 ^ b0;            // -b
            uint64_t wa0 = nb0, wa1 = nb1;                // real part of w*y
            uint64_t wb0 = a0, wb1 = a1; z4add(wb0, wb1, nb0, nb1);   // a - b
            uint64_t r0 = 0, r1 = 0, s0 = 0, s1 = 0;
            auto addc = [&](int k, uint64_t x0, uint64_t x1, uint64_t& o0, uint64_t& o1) {   // o += k * x  (k in Z/4)
                if (k == 1) z4add(o0, o1, x0, x1); else if (k == 2) { o1 ^= x0; } else if (k == 3) z4add(o0, o1, x0, x1 ^ x0); };
            addc(c, a0, a1, r0, r1); addc(c, b0, b1, s0, s1);
            addc(d, wa0, wa1, r0, r1); addc(d, wb0, wb1, s0, s1);
            O[0][w] = r0; O[1][w] = r1; O[2][w] = s0; O[3][w] = s1;
        }
    }
    // F4 helpers (planes a, b; element a + b w)
    static void f4axpy(vector<uint64_t>& A, vector<uint64_t>& Bp, const uint64_t* X, const uint64_t* Y, int c, int d, int W) {   // (A,B) += (c + d w) * (X,Y)
        for (int w = 0; w < W; w++) { uint64_t x = X[w], y = Y[w]; uint64_t na = (c ? x : 0) ^ (d ? y : 0), nb = (c ? y : 0) ^ (d ? x : 0) ^ (d ? y : 0); A[w] ^= na; Bp[w] ^= nb; }
    }
    void dense_gr(long long sparseStagePivots) {
        vector<int> cmap(nc, -1); int m = 0;
        for (int c = 0; c < nc; c++) if (!coldead[c] && colcnt[c] > 0) cmap[c] = m++;
        vector<vector<int>>().swap(colrows);
        int W = (m + 63) / 64;
        fprintf(stderr, "  [gr] m=%d, basis memory <= %.2f GB\n", m, 4.0 * m * W * 8 / 1e9);
        vector<R4> B; vector<int> rowOfPiv(m, -1);
        vector<uint64_t> xs(4 * (size_t)W); uint64_t* X[4] = {xs.data(), xs.data() + W, xs.data() + 2 * W, xs.data() + 3 * W};
        vector<uint64_t> accs(64 * (size_t)W), tmps(4 * (size_t)W); uint64_t* T[4] = {tmps.data(), tmps.data() + W, tmps.data() + 2 * W, tmps.data() + 3 * W};
        vector<uint64_t> Ms(64 * (size_t)W);   // 16 multiples of the new pivot row
        auto MP = [&](int f, int i) { return Ms.data() + (4 * f + i) * (size_t)W; };
        auto AP = [&](int f, int i) { return accs.data() + (4 * f + i) * (size_t)W; };
        vector<vector<uint64_t>> Ya, Yb; vector<int> ypiv;
        auto f4get = [&](const vector<uint64_t>& a, const vector<uint64_t>& b, int c) { int w = c >> 6, s = c & 63; return (int)((a[w] >> s) & 1) | (int)(((b[w] >> s) & 1) << 1); };
        auto f4inv = [](int v) { return v == 1 ? 1 : (v == 2 ? 3 : 2); };   // codes: 1 = 1, 2 = w, 3 = 1 + w = w^2
        auto insertY = [&](vector<uint64_t>& a, vector<uint64_t>& b) {
            for (size_t i = 0; i < Ya.size(); i++) { int v = f4get(a, b, ypiv[i]); if (v) f4axpy(a, b, Ya[i].data(), Yb[i].data(), v & 1, v >> 1, W); }
            int c = -1; for (int w = 0; w < W && c < 0; w++) { uint64_t q = a[w] | b[w]; if (q) c = w * 64 + __builtin_ctzll(q); }
            if (c < 0) return;
            int iv = f4inv(f4get(a, b, c)); if (iv != 1) { vector<uint64_t> na(W, 0), nb(W, 0); f4axpy(na, nb, a.data(), b.data(), iv & 1, iv >> 1, W); a.swap(na); b.swap(nb); }
            for (size_t i = 0; i < Ya.size(); i++) { int v = f4get(Ya[i], Yb[i], c); if (v) f4axpy(Ya[i], Yb[i], a.data(), b.data(), v & 1, v >> 1, W); }
            Ya.push_back(a); Yb.push_back(b); ypiv.push_back(c);
        };
        long long done = 0, stagePiv = sparseStagePivots;
        for (int r = 0; r < nr; r++) {
            if (rowdead[r] || rows[r].empty()) continue;
            vector<pair<int,int>> ents; for (auto e : rows[r]) { int c = cmap[col(e)]; if (c < 0) { fprintf(stderr, "gr: dead col\n"); exit(6); } ents.push_back({c, (int)val(e)}); }
            vector<uint32_t>().swap(rows[r]);
            fill(xs.begin(), xs.end(), 0); for (auto& [c, v] : ents) grset(X, c, v);
            bool used[16] = {false};
            for (auto& [c, v] : ents) { int br = rowOfPiv[c]; if (br < 0) continue; int f = grget(X, c); if (!f) continue; int g = grt.neg[0][f];
                uint64_t* A[4] = {AP(g, 0), AP(g, 1), AP(g, 2), AP(g, 3)};
                if (!used[g]) { for (int i = 0; i < 4; i++) memcpy(A[i], B[br].p[i], W * 8); used[g] = true; } else gradd(A, B[br].p, W); }
            for (int g = 1; g < 16; g++) if (used[g]) { uint64_t* A[4] = {AP(g, 0), AP(g, 1), AP(g, 2), AP(g, 3)}; grscale(T, A, g, W); gradd(X, T, W); }
            int j = -1; for (int w = 0; w < W && j < 0; w++) { uint64_t q = X[0][w] | X[2][w]; if (q) j = w * 64 + __builtin_ctzll(q); }
            if (j >= 0) {
                int u = grt.inv[0][grget(X, j)];
                if (u != 1) { grscale(T, X, u, W); for (int i = 0; i < 4; i++) memcpy(X[i], T[i], W * 8); }
                for (int f = 0; f < 16; f++) { uint64_t* O[4] = {MP(f, 0), MP(f, 1), MP(f, 2), MP(f, 3)}; grscale(O, X, f, W); }
                long long nb = B.size();
                #pragma omp parallel for schedule(static, 64)
                for (long long b = 0; b < nb; b++) { int f = grget(B[b].p, j); if (f) { int g = grt.neg[0][f]; uint64_t* Y[4] = {MP(g, 0), MP(g, 1), MP(g, 2), MP(g, 3)}; gradd(B[b].p, Y, W); } }
                R4 nr4; for (int i = 0; i < 4; i++) { nr4.p[i] = (uint64_t*)malloc(W * 8); if (!nr4.p[i]) { fprintf(stderr, "gr: out of memory\n"); exit(9); } memcpy(nr4.p[i], X[i], W * 8); }
                rowOfPiv[j] = B.size(); B.push_back(nr4); stagePiv++; pivots++; rowdead[r] = 1;
            } else {
                bool any = false; for (int w = 0; w < W; w++) if (X[1][w] | X[3][w]) { any = true; break; }
                if (any) { vector<uint64_t> a(X[1], X[1] + W), b(X[3], X[3] + W); insertY(a, b); }
            }
            if (++done % 200000 == 0) fprintf(stderr, "  [gr] %lld rows, basis %zu, leftover span %zu\n", done, B.size(), Ya.size());
        }
        stages.push_back(stagePiv);
        fprintf(stderr, "  [gr] stage 0 done: basis %zu, leftover span %zu\n", B.size(), Ya.size());
        vector<vector<uint64_t>> F0, F1;
        for (size_t i = 0; i < Ya.size(); i++) {
            vector<uint64_t> a = Ya[i], b = Yb[i]; vector<pair<int,int>> hits;
            for (int w = 0; w < W; w++) { uint64_t q = a[w] | b[w]; while (q) { int c = w * 64 + __builtin_ctzll(q); q &= q - 1; if (rowOfPiv[c] >= 0) hits.push_back({c, f4get(a, b, c)}); } }
            for (auto& [c, v] : hits) { auto& bb = B[rowOfPiv[c]]; f4axpy(a, b, bb.p[0], bb.p[2], v & 1, v >> 1, W); }   // b mod 2 = planes (a0, b0)
            bool nz = false; for (int w = 0; w < W; w++) { uint64_t q = a[w] | b[w]; while (q) { int c = w * 64 + __builtin_ctzll(q); q &= q - 1; if (rowOfPiv[c] >= 0) { fprintf(stderr, "gr: pivot col nonzero\n"); exit(7); } nz = true; } }
            if (nz) { F0.push_back(move(a)); F1.push_back(move(b)); }
        }
        for (auto& r4 : B) for (int i = 0; i < 4; i++) free(r4.p[i]);
        long long nrw = F0.size(), r1 = 0;
        for (int c = 0; c < m && r1 < nrw; c++) {
            long long piv = -1; for (long long i = r1; i < nrw; i++) if (f4get(F0[i], F1[i], c)) { piv = i; break; }
            if (piv < 0) continue;
            swap(F0[piv], F0[r1]); swap(F1[piv], F1[r1]);
            int iv = f4inv(f4get(F0[r1], F1[r1], c)); if (iv != 1) { vector<uint64_t> na(W, 0), nb(W, 0); f4axpy(na, nb, F0[r1].data(), F1[r1].data(), iv & 1, iv >> 1, W); F0[r1].swap(na); F1[r1].swap(nb); }
            for (long long i = 0; i < nrw; i++) if (i != r1) { int v = f4get(F0[i], F1[i], c); if (v) f4axpy(F0[i], F1[i], F0[r1].data(), F1[r1].data(), v & 1, v >> 1, W); }
            r1++;
        }
        if (r1) { stages.push_back(r1); pivots += r1; }
        fprintf(stderr, "  [gr] stage 1 (R/2 summands): %lld\n", r1);
    }
