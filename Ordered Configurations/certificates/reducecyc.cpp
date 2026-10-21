// =====================================================================================================
//  reducecyc.cpp  --  coordinates of an integral cycle in a basis of H_deg, by exact reduction modulo boundaries
// =====================================================================================================
//  For a graph whose maps d_1 .. d_{deg+1} all reduce completely with pivots +-1 (so H_deg is free and has an
//  explicit basis), this reduces a given integral deg-cycle modulo boundaries:
//   1. d_1 .. d_deg are eliminated with pivots +-1 and their pivot rows dropped from the next map (valid for
//      cycles: a cycle is determined by its values on the kept cells, the dropped cells carrying a minor of
//      determinant +-1).  The kept deg-cells then give coordinates on Z_deg.
//   2. d_{deg+1} on the kept cells is eliminated with pivots +-1.  If nothing is left over, the non-pivot
//      columns index a basis of H_deg, and reducing the cycle against the pivot rows (in pivot order) leaves
//      exactly its coordinates in that basis.
//  The program prints the number of non-pivot columns (= rank of H_deg when the residual is empty), the
//  number of residual entries (0 means the statement above is exact) and the size of the result.
//  USAGE.  reducecyc deg cycle_file out_file < graph_with_n      (cycle file: "coef cell" lines)
//  OUTPUT. out_file: "coef cell" lines, one per nonzero coordinate (the cell names the basis element).
// =====================================================================================================
#include "swiat_cells.h"
#include "elim_z.h"

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: reducecyc deg cycle_file out_file < graph\n"); return 1; }
    int deg = atoi(argv[1]);
    string tfn = argv[2], outp = argv[3];
    if (!read_graph()) return 1;
    vector<vector<Key>> C(deg + 2);
    for (int k = 0; k <= deg + 1; k++) { enumerate(k, C[k]); sort(C[k].begin(), C[k].end()); fprintf(stderr, "C_%d: %zu\n", k, C[k].size()); }
    auto idx = [&](int k, Key key) { int j = lookup(C[k], key); if (j < 0) { fprintf(stderr, "missing cell\n"); exit(2); } return j; };
    vector<pair<Key, int>> bd;

    // ---- step 1: d_1 .. d_deg
    vector<char> drop;
    for (int k = 1; k <= deg; k++) {
        ZElim M(C[k].size(), C[k - 1].size());
        for (size_t i = 0; i < C[k].size(); i++) {
            boundary(decode(C[k][i]), bd);
            for (auto &x : bd) { int j = idx(k - 1, x.first); if (drop.empty() || !drop[j]) M.rows[i].push_back({j, (ent)x.second}); }
        }
        M.finalize(); M.run();
        long long left = 0;
        for (int r = 0; r < M.nr; r++) if (!M.rowdead[r]) left += M.rows[r].size();
        fprintf(stderr, "d_%d: %zu pivots, residual entries %lld\n", k, M.piv.size(), left);
        drop.assign(C[k].size(), 0);
        for (auto &pr : M.piv) drop[pr.second] = 1;
    }
    vector<int> colof(C[deg].size(), -1), colcell;
    for (size_t j = 0; j < C[deg].size(); j++) if (!drop[j]) { colof[j] = colcell.size(); colcell.push_back(j); }

    // ---- step 2: d_{deg+1} on the kept cells
    ZElim D(C[deg + 1].size(), colcell.size());
    for (size_t i = 0; i < C[deg + 1].size(); i++) {
        boundary(decode(C[deg + 1][i]), bd);
        for (auto &x : bd) { int j = colof[idx(deg, x.first)]; if (j >= 0) D.rows[i].push_back({j, (ent)x.second}); }
    }
    D.finalize(); D.run();
    long long left = 0, freecols = 0;
    for (int r = 0; r < D.nr; r++) if (!D.rowdead[r]) left += D.rows[r].size();
    for (int c = 0; c < D.nc; c++) if (!D.coldead[c]) freecols++;
    fprintf(stderr, "d_%d on kept cells: %zu pivots, residual entries %lld, non-pivot columns %lld\n", deg + 1, D.piv.size(), left, freecols);

    // ---- the cycle, reduced
    vector<I> x(D.nc, 0), g(D.nr, 0);
    {
        ifstream in(tfn); string line; Cell c;
        while (getline(in, line)) {
            if (!is_coef_line(line)) continue;
            size_t sp = line.find(' ');
            I coef = s2i(line.substr(0, sp));
            if (!parse_cell(line.substr(sp + 1), c)) { fprintf(stderr, "bad cell\n"); return 2; }
            int j = colof[idx(deg, encode(c))];
            if (j >= 0) x[j] += coef;
        }
    }
    D.reduce(x, g);
    FILE *fo = fopen(outp.c_str(), "w");
    long long nz = 0; I mx = 0;
    for (int c = 0; c < D.nc; c++)
        if (x[c]) {
            if (D.coldead[c]) { fprintf(stderr, "BUG: entry left on a pivot column\n"); return 3; }
            fprintf(fo, "%s %s\n", i2s(x[c]).c_str(), describe(decode(C[deg][colcell[c]])).c_str());
            nz++;
            I a = x[c] < 0 ? -x[c] : x[c];
            if (a > mx) mx = a;
        }
    fclose(fo);
    printf("{\"nonpivot_columns\":%lld,\"residual_entries\":%lld,\"coords_nonzero\":%lld,\"max_coef\":\"%s\"}\n", freecols, left, nz, i2s(mx).c_str());
    return 0;
}
