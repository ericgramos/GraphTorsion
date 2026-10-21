// =====================================================================================================
//  reducemod.cpp  --  reduce many cycles modulo the boundaries of Gamma once, for fast span questions
// =====================================================================================================
//  Writes every input cycle as a short vector in the quotient  C_deg / B_deg  localised at p, modulo Q = p^K, so
//  that questions about spans of cycles modulo boundaries become small linear algebra (certificates/spanq).
//
//  METHOD.  As in orderq (no deleted edges): d_1 .. d_deg are reduced mod Q with unit pivots and their pivot rows
//  dropped (valid for cycles), and D = d_{deg+1} on the kept deg-cells is eliminated mod Q with unit pivots,
//  keeping a snapshot of every pivot row.  Modulo the span of the snapshots, every vector is congruent to a unique
//  vector supported on the LIVE (non-pivot) columns, obtained by reducing it against the snapshots in pivot order
//  (fact (3) of elim_modq.h).  The reduction here is sparse: a heap holds the pivot positions of the nonzero pivot
//  columns, and only snapshots that are actually needed are applied, so a cycle with few cells costs little.
//  The rows that survive the elimination (the residual rows R) are written too.  Then for any set S of cycles
//        C_deg / (B_deg + span S)   =   (Z/Q)^{live} / (span R + span of the reduced vectors of S).
//
//  USAGE.  reducemod deg p K blocks_file out_file < graph_with_n
//          blocks_file: blocks "Z header..." followed by "coef cell" lines (other tags are skipped); cells may be
//          written for a prefix of the edge list, in the default basis.
//  OUTPUT. out_file, text:
//            LIVE L Q                  number of live columns and the modulus
//            R                         one residual row per line: "col:val col:val ..." (live column indices)
//            ...
//            B <header of the block>   then one line "col:val ..." with the reduced vector
//  EXAMPLE. reducemod 2 2 8 tori.txt tori.red < graphs/K4dbl5_Hfirst_n5.txt
// =====================================================================================================
#include "swiat_cells.h"
#include "elim_modq.h"

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: reducemod deg p K blocks out < graph\n"); return 1; }
    int deg = atoi(argv[1]);
    PP = strtoull(argv[2], 0, 10);
    int KK = atoi(argv[3]);
    Q = ppow(KK);
    if (Q >= (1ULL << 31)) { fprintf(stderr, "modulus too large\n"); return 1; }
    string bfn = argv[4], ofn = argv[5];
    if (!read_graph()) return 1;
    vector<vector<Key>> C(deg + 2);
    for (int k = 0; k <= deg + 1; k++) { enumerate(k, C[k]); sort(C[k].begin(), C[k].end()); fprintf(stderr, "C_%d: %zu\n", k, C[k].size()); }
    vector<pair<Key, int>> bd;

    // ---- dropping (d_1 .. d_deg), as in orderq
    vector<char> drop;
    for (int k = 1; k <= deg; k++) {
        QElim M(C[k].size(), C[k - 1].size(), false);
        for (size_t i = 0; i < C[k].size(); i++) {
            boundary(decode(C[k][i]), bd);
            for (auto &x : bd) { int j = lookup(C[k - 1], x.first); if (drop.empty() || !drop[j]) M.rows[i].push_back({j, modq_ll(x.second)}); }
        }
        M.finalize(); M.run();
        long long left = 0;
        for (size_t i = 0; i < C[k].size(); i++) if (!M.rowdead[i]) left += M.rows[i].size();
        if (left) fprintf(stderr, "WARNING: d_%d does not reduce completely\n", k);
        drop.assign(C[k].size(), 0);
        for (size_t i = 0; i < C[k].size(); i++) if (M.rowdead[i]) drop[i] = 1;
    }
    vector<int> colof(C[deg].size(), -1);
    int ncol = 0;
    for (size_t j = 0; j < C[deg].size(); j++) if (!drop[j]) colof[j] = ncol++;

    // ---- eliminate D
    int K1 = deg + 1;
    QElim D(C[K1].size(), ncol, true);
    for (size_t i = 0; i < C[K1].size(); i++) {
        boundary(decode(C[K1][i]), bd);
        for (auto &x : bd) { int j = colof[lookup(C[deg], x.first)]; if (j >= 0) D.rows[i].push_back({j, modq_ll(x.second)}); }
    }
    D.finalize();
    D.run();
    vector<Op>().swap(D.ops);
    vector<int> pivorder(ncol, -1), liveidx(ncol, -1);
    for (size_t o = 0; o < D.snapshots().size(); o++) pivorder[D.snapshots()[o].first] = o;
    int L = 0;
    for (int c = 0; c < ncol; c++) if (!D.coldead[c]) liveidx[c] = L++;
    vector<u64> pinv(D.snapshots().size());
    for (size_t o = 0; o < D.snapshots().size(); o++) { int pc = D.snapshots()[o].first; for (auto &x : D.snapshots()[o].second) if (x.first == pc) pinv[o] = egcd_inv(x.second); }
    fprintf(stderr, "unit rank %lld, live columns %d\n", D.pivots, L);

    FILE *fo = fopen(ofn.c_str(), "w");
    fprintf(fo, "LIVE %d %llu\n", L, (unsigned long long)Q);
    long long nres = 0;
    for (int r = 0; r < D.nr; r++)
        if (!D.rowdead[r] && !D.rows[r].empty()) {
            fprintf(fo, "R");
            for (auto &x : D.rows[r]) fprintf(fo, " %d:%u", liveidx[x.first], x.second);
            fprintf(fo, "\n");
            nres++;
        }
    fprintf(stderr, "residual rows %lld\n", nres);

    // ---- reduce the blocks
    unordered_map<int, u64> z;                   // sparse vector over the kept columns
    string header;
    long long nb = 0;
    auto finish = [&]() {
        if (header.empty()) return;
        priority_queue<int, vector<int>, greater<int>> heap;   // pivot positions with a nonzero entry
        unordered_set<int> inheap;
        for (auto &kv : z) if (kv.second && pivorder[kv.first] >= 0) { heap.push(pivorder[kv.first]); inheap.insert(pivorder[kv.first]); }
        while (!heap.empty()) {
            int o = heap.top(); heap.pop(); inheap.erase(o);
            int pc = D.snapshots()[o].first;
            auto it = z.find(pc);
            if (it == z.end() || it->second == 0) continue;
            u64 m = mulq(it->second, pinv[o]);
            for (auto &x : D.snapshots()[o].second) {
                u64 &v = z[x.first];
                v = (v + Q - mulq(m, x.second)) % Q;
                int po = pivorder[x.first];
                if (v && po > o && !inheap.count(po)) { heap.push(po); inheap.insert(po); }
            }
            z[pc] = 0;
        }
        fprintf(fo, "B %s\n", header.c_str());
        vector<pair<int, u64>> out;
        for (auto &kv : z) if (kv.second) {
            if (liveidx[kv.first] < 0) { fprintf(stderr, "BUG: pivot column left\n"); exit(3); }
            out.push_back({liveidx[kv.first], kv.second});
        }
        sort(out.begin(), out.end());
        for (auto &x : out) fprintf(fo, "%d:%llu ", x.first, (unsigned long long)x.second);
        fprintf(fo, "\n");
        z.clear(); header.clear(); nb++;
        if (nb % 5000 == 0) fprintf(stderr, "%lld blocks reduced\n", nb);
    };
    ifstream in(bfn); string line; Cell c;
    while (getline(in, line)) {
        if (line.rfind("Z ", 0) == 0) { finish(); header = line.substr(2); continue; }
        if (header.empty() || !is_coef_line(line)) continue;
        size_t sp = line.find(' ');
        long long coef = atoll(line.substr(0, sp).c_str());
        if (!parse_cell(line.substr(sp + 1), c)) { fprintf(stderr, "bad cell\n"); return 2; }
        int jc = lookup(C[deg], encode(c));
        if (jc < 0) { fprintf(stderr, "cell not in the complex: %s\n", line.c_str()); return 2; }
        int j = colof[jc];
        if (j >= 0) z[j] = (z[j] + modq_ll(coef)) % Q;
    }
    finish();
    fclose(fo);
    fprintf(stderr, "%lld blocks written\n", nb);
    return 0;
}
