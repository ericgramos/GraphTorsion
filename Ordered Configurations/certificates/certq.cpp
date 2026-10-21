// =====================================================================================================
//  certq.cpp  --  explicit generators and certificates for torsion of prime-power order p^k
// =====================================================================================================
//  WHAT IT PRODUCES.  For a prime p, an exponent k and a degree deg, an integral (deg+1)-chain w and a
//  Z/p^k-valued deg-cochain f such that
//    (1) d w = 0 (mod p^k), so z = d w / p^k is an integral deg-cycle and p^k z is a boundary;
//    (2) f(d s) = 0 (mod p^k) for every (deg+1)-cell s, i.e. f is a cocycle mod p^k;
//    (3) f(z) is prime to p.
//  Then [z] has order EXACTLY p^k in H_deg(Conf_n Gamma; Z): if p^(k-1) z were a boundary d u, then
//  p^(k-1) f(z) = f(d u) = 0 (mod p^k), so p would divide f(z).  The certificate file has the format of
//  certp.cpp with "PRIME p^k" and is checked by certificates/verify_cert.py, which shares no code with this
//  program.
//
//  ALGORITHM.  Work over Q = p^K with K = k + 1 (or larger).
//   1. Reduce d_1 .. d_deg mod Q with unit pivots, dropping the pivot rows of d_j from the columns of d_{j+1}
//      exactly as the engine does.  These maps should reduce completely (no p-torsion below degree deg);
//      a warning is printed otherwise.  Dropping is valid for cycles: a deg-cycle is determined by its values
//      on the kept cells, because the dropped cells carry a minor of d_deg that is a unit mod p.
//   2. Build D = d_{deg+1} restricted to the kept deg-cells (rows = (deg+1)-cells, columns = kept deg-cells) and
//      eliminate it mod Q with unit pivots, logging every row operation and snapshotting every pivot row.
//   3. The rows and columns left alive form a residual block R whose entries are all divisible by p.  Its local
//      Smith form U R V = diag(p^dv_t * unit) is computed densely with both transformations.
//   4. Pick t with dv_t = k.  u = U_t (a combination of residual rows) satisfies u R = p^k * unit * (V^-1)_t.
//      Since live rows vanish on pivot columns, the same combination of FINAL rows of D has boundary supported
//      on residual columns; replaying the log backwards turns it into a chain w on the original (deg+1)-cells
//      with d w = p^k * (a cycle).  This z = d w / p^k is the generator of the t-th cyclic summand.
//   5. g = V_t on residual columns satisfies R g = 0 mod p^k.  It is extended to the pivot columns by
//      back-substitution through the pivot snapshots in reverse pivot order (snapshot o has entries only in
//      its own column, in later pivot columns and in residual columns), making D' f = 0 mod Q for every final
//      row, hence D f = 0.  On dropped cells f = 0.  By construction f(z) = (U R V)_tt / p^k = unit mod p.
//   6. Before writing, the program recomputes d w over the integers on all deg-cells, checks divisibility by
//      p^k and that f(z) is prime to p.  Correctness never depends on any of this: the verifier rechecks.
//
//  USAGE.   certq deg p k index outfile [K] < graph_with_n
//             index >= 0: certificate for the index-th diagonal entry of valuation k (0 = the first);
//             index = -1: GENERATOR MODE, see below.  K defaults to k + 1.
//  INPUT.   "V E n" followed by E lines "tail head".
//  OUTPUT.  certificate mode: the certificate file and one JSON line {certificate, order, residual, w_cells,
//           f_cells, z_support, pairing_mod_p}.
//           generator mode: for every diagonal entry with 1 <= dv_t < K a block
//               Z t p^dv_t            followed by "coef cell" lines: the integral cycle z_t
//               W coef cell ...        a chain w_t with d w_t = p^dv_t z_t (coefficients in [0, Q))
//               F coef cell ...        only when dv_t = k: a detecting cocycle mod p^k
//           These z_t generate the p-primary torsion of order < p^K (one generator per cyclic summand).
//  EXAMPLE. printf "4 11 5\n0 1\n0 2\n0 3\n1 2\n1 3\n2 3\n0 1\n0 2\n0 3\n1 2\n1 3\n" | certq 2 2 4 0 cert16.txt
//           (an element of order 16 in H_2(Conf_5) of K_4 with five doubled edges; about 4 minutes, 3 GB).
// =====================================================================================================
#include "swiat_cells.h"
#include "elim_modq.h"

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: certq deg p k index outfile [K] < graph\n"); return 1; }
    int deg = atoi(argv[1]);
    PP = strtoull(argv[2], 0, 10);
    int kk = atoi(argv[3]);
    int pick = atoi(argv[4]);
    string outfn = argv[5];
    int KK = argc > 6 ? atoi(argv[6]) : kk + 1;
    if (KK < kk + 1) { fprintf(stderr, "need K >= k+1\n"); return 1; }
    Q = ppow(KK);
    if (Q >= (1ULL << 31)) { fprintf(stderr, "modulus too large (entries are stored in 32 bits)\n"); return 1; }
    u64 PK = ppow(kk);
    if (!read_graph()) return 1;

    vector<vector<Key>> C(deg + 2);
    for (int k = 0; k <= deg + 1; k++) {
        enumerate(k, C[k]);
        sort(C[k].begin(), C[k].end());
        fprintf(stderr, "C_%d: %zu\n", k, C[k].size());
    }
    auto idx = [&](int k, Key key) {
        int j = lookup(C[k], key);
        if (j < 0) { fprintf(stderr, "missing face\n"); exit(2); }
        return j;
    };
    vector<pair<Key, int>> bd;

    // ---- step 1: d_1 .. d_deg mod Q, dropping pivot rows from the next map
    vector<char> drop;
    for (int k = 1; k <= deg; k++) {
        QElim M(C[k].size(), C[k - 1].size(), false);
        for (size_t i = 0; i < C[k].size(); i++) {
            boundary(decode(C[k][i]), bd);
            for (auto &x : bd) {
                int j = idx(k - 1, x.first);
                if (drop.empty() || !drop[j]) M.rows[i].push_back({j, modq_ll(x.second)});
            }
        }
        M.finalize();
        M.run();
        long long leftover = 0;
        for (size_t i = 0; i < C[k].size(); i++) if (!M.rowdead[i]) leftover += M.rows[i].size();
        drop.assign(C[k].size(), 0);
        for (size_t i = 0; i < C[k].size(); i++) if (M.rowdead[i]) drop[i] = 1;
        fprintf(stderr, "d_%d mod %llu: unit rank %lld, residual entries %lld\n", k, (unsigned long long)Q, M.pivots, leftover);
        if (leftover) fprintf(stderr, "WARNING: d_%d does not reduce completely; the certificate may fail verification\n", k);
    }

    // ---- step 2: D = d_{deg+1} on the kept deg-cells, eliminated with an operation log
    int K1 = deg + 1;
    vector<int> keepidx(C[deg].size(), -1), kept;
    for (size_t j = 0; j < C[deg].size(); j++) if (!drop[j]) { keepidx[j] = kept.size(); kept.push_back(j); }
    QElim D(C[K1].size(), kept.size(), true);
    for (size_t i = 0; i < C[K1].size(); i++) {
        boundary(decode(C[K1][i]), bd);
        for (auto &x : bd) { int j = idx(deg, x.first); if (keepidx[j] >= 0) D.rows[i].push_back({keepidx[j], modq_ll(x.second)}); }
    }
    D.finalize();
    fprintf(stderr, "eliminating d_%d: %d x %d\n", K1, D.nr, D.nc);
    D.run();
    fprintf(stderr, "unit rank %lld, %zu row operations\n", D.pivots, D.ops.size());

    // ---- step 3: residual block and its local Smith form with transformations
    vector<int> rr, rc, cpos;
    D.residual(rr, rc, cpos);
    int nr = rr.size(), nc = rc.size();
    fprintf(stderr, "residual block %d x %d\n", nr, nc);
    vector<vector<u64>> R(nr, vector<u64>(nc, 0)), U, Vm;
    for (int i = 0; i < nr; i++) for (auto &x : D.rows[rr[i]]) R[i][cpos[x.first]] = x.second;
    vector<int> dv = local_smith(R, nr, nc, U, Vm, true, KK);
    map<int, int> hist; for (int v : dv) hist[v]++;
    fprintf(stderr, "local Smith form of the residual (valuation: count):");
    for (auto &h : hist) fprintf(stderr, " %d:%d", h.first, h.second);
    fprintf(stderr, "  (+%d of valuation >= %d)\n", nc - (int)dv.size(), KK);

    // chain w_t (step 4): coefficients on the original (deg+1)-cells, mod Q
    auto chain_of = [&](int t) {
        vector<u64> c(D.nr, 0);
        for (int i = 0; i < nr; i++) c[rr[i]] = U[t][i];
        for (size_t o = D.ops.size(); o-- > 0;) { auto &op = D.ops[o]; if (c[op.t]) c[op.s] = (c[op.s] + Q - mulq(op.m, c[op.t])) % Q; }
        return c;
    };
    // cocycle f_t (step 5): values on the kept deg-cells, mod Q
    auto cocycle_of = [&](int t) {
        vector<u64> f(D.nc, 0);
        for (int j = 0; j < nc; j++) f[rc[j]] = Vm[j][t];
        for (size_t o = D.snapshots().size(); o-- > 0;) {
            int pc = D.snapshots()[o].first; u64 acc = 0, sc = 0;
            for (auto &x : D.snapshots()[o].second) { if (x.first == pc) sc = x.second; else acc = (acc + mulq(x.second, f[x.first])) % Q; }
            f[pc] = mulq((Q - acc) % Q, egcd_inv(sc));
        }
        return f;
    };
    // integral boundary of a chain, on all deg-cells
    auto bd_of = [&](const vector<u64> &c) {
        map<int, long long> dw;
        for (int i = 0; i < D.nr; i++)
            if (c[i]) { boundary(decode(C[K1][i]), bd); for (auto &x : bd) dw[idx(deg, x.first)] += (long long)c[i] * x.second; }
        return dw;
    };

    if (pick < 0) {
        // ---- generator mode
        FILE *fo = fopen(outfn.c_str(), "w");
        fprintf(fo, "GRAPH %d %d %d\n", V, E, N);
        for (int e = 0; e < E; e++) fprintf(fo, "EDGE %d %d\n", etail[e], ehead[e]);
        fprintf(fo, "DEGREE %d\n", deg);
        int ng = 0;
        for (int t = 0; t < (int)dv.size(); t++) {
            if (dv[t] < 1 || dv[t] >= KK) continue;
            vector<u64> c = chain_of(t);
            auto dw = bd_of(c);
            long long pv = ppow(dv[t]);
            fprintf(fo, "Z %d %llu\n", ng, (unsigned long long)pv);
            long long sz = 0;
            for (auto &kv : dw) if (kv.second) {
                if (kv.second % pv) { fprintf(stderr, "BUG: boundary not divisible\n"); return 3; }
                fprintf(fo, "%lld %s\n", kv.second / pv, describe(decode(C[deg][kv.first])).c_str());
                sz++;
            }
            for (int i = 0; i < D.nr; i++) if (c[i]) fprintf(fo, "W %llu %s\n", (unsigned long long)c[i], describe(decode(C[K1][i])).c_str());
            if (dv[t] == kk) {
                vector<u64> f = cocycle_of(t);
                for (int j = 0; j < D.nc; j++)
                    if (f[j] % (u64)pv) fprintf(fo, "F %llu %s\n", (unsigned long long)(f[j] % (u64)pv), describe(decode(C[deg][kept[j]])).c_str());
            }
            fprintf(stderr, "generator %d: order %lld, support %lld cells\n", ng, pv, sz);
            ng++;
        }
        fclose(fo);
        printf("{\"generators\":\"%s\",\"count\":%d}\n", outfn.c_str(), ng);
        return 0;
    }

    // ---- certificate mode
    int t = -1, seen = 0;
    for (int i = 0; i < (int)dv.size(); i++) if (dv[i] == kk && seen++ == pick) { t = i; break; }
    if (t < 0) { printf("{\"certificate\":null,\"reason\":\"no diagonal entry of valuation %d\"}\n", kk); return 1; }
    vector<u64> c = chain_of(t);
    vector<u64> f = cocycle_of(t);
    // ---- step 6: self-check over the integers
    auto dw = bd_of(c);
    long long nz = 0;
    __int128 pairing = 0;
    for (auto &kv : dw) {
        if (kv.second % (long long)PK) { fprintf(stderr, "BUG: d w not divisible by p^k\n"); return 3; }
        if (kv.second) nz++;
        if (keepidx[kv.first] >= 0) pairing += (__int128)(kv.second / (long long)PK) * (long long)f[keepidx[kv.first]];
    }
    long long pm = (long long)(((pairing % (__int128)PP) + PP) % PP);
    if (pm == 0) { fprintf(stderr, "BUG: pairing divisible by p\n"); return 3; }
    FILE *fo = fopen(outfn.c_str(), "w");
    fprintf(fo, "GRAPH %d %d %d\n", V, E, N);
    for (int e = 0; e < E; e++) fprintf(fo, "EDGE %d %d\n", etail[e], ehead[e]);
    fprintf(fo, "DEGREE %d\nPRIME %llu\n", deg, (unsigned long long)PK);
    long long nw = 0, nf = 0;
    for (int i = 0; i < D.nr; i++) if (c[i]) { fprintf(fo, "W %llu %s\n", (unsigned long long)c[i], describe(decode(C[K1][i])).c_str()); nw++; }
    for (int j = 0; j < D.nc; j++) if (f[j] % PK) { fprintf(fo, "F %llu %s\n", (unsigned long long)(f[j] % PK), describe(decode(C[deg][kept[j]])).c_str()); nf++; }
    fclose(fo);
    printf("{\"certificate\":\"%s\",\"order\":%llu,\"residual\":[%d,%d],\"w_cells\":%lld,\"f_cells\":%lld,\"z_support\":%lld,\"pairing_mod_p\":%lld}\n",
           outfn.c_str(), (unsigned long long)PK, nr, nc, nw, nf, nz, pm);
    return 0;
}
