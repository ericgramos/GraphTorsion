// =====================================================================================================
//  exactsolve.cpp  --  write a given cycle x as y + d u over Z, with y supported on a subgraph H
// =====================================================================================================
//  SETTING.  The edges of H must be the FIRST h edges of the graph.  Then the base half-edge of every vertex
//  of H lies in H and the reduced complex of H is the subcomplex C(H) of cells that use only those edges
//  (no particle on a later edge, no generator half-edge on a later edge).  The input x is an integral
//  deg-cycle of Gamma whose class is known (from orderq) to lie in the image of H_deg(Conf_n H).  The output is
//  an integral (deg+1)-chain u and the cycle y = x - d u, which is checked to vanish on every cell outside H.
//  Then y is an explicit integral cycle of Conf_n(H) with i_*[y] = [x].
//
//  METHOD.
//   1. Dropping.  Eliminate d_deg restricted to (deg-cells outside H) x ((deg-1)-cells outside H) modulo a large
//      prime.  The pivot rows S carry a square minor that is nonzero mod p, hence has nonzero determinant over
//      Z.  Cells of H have no boundary outside H (C(H) is a subcomplex).  So if c is a deg-cycle vanishing on
//      every cell outside H except possibly on S, then restricting d c = 0 to the minor's columns gives
//      c_S * (nonsingular matrix) = 0, so c_S = 0.  Hence it is enough to make x - d u vanish outside H and
//      outside S.  This roughly halves the number of columns and keeps the residual block small.  (With h = 0
//      the subgraph is empty and the program decides whether x itself is a boundary, writing x = d u.)
//   2. Exact elimination over Z, pivots +-1 only (elim_z.h), of D = d_{deg+1} restricted to the remaining
//      columns, with an operation log.
//   3. Reduce x (restricted to these columns) against the pivot rows in pivot order; the multiples subtracted
//      are recorded as coefficients g on the final rows.  What is left of x lives on the residual columns
//      and must be an integral combination lam of the residual rows.  That small system is written to
//      out.res and solved by solve_residual.py (Hermite normal form); lam is added to g.
//   4. Replay the log backwards: u = g L is integral and x - d u vanishes on the columns of step 2.
//   5. Compute y = x - d u on all deg-cells exactly and check that it vanishes outside H (by step 1 it must).
//  All arithmetic is exact (32-bit matrix entries with overflow checks, 128-bit vectors).
//  The answer is NOT unique: y and u depend on the pivot choices (and on the dropped set S), so a different version or
//  ordering can return a different, equally valid y.  Check an answer by building a certificate, not by comparing.
//
//  USAGE.  exactsolve deg h target_file out_prefix < graph_with_n
//          target_file: "coef cell" lines, cells in the default basis of this edge order.
//          The environment variable SOLVER may give the path of solve_residual.py (default: current folder).
//  OUTPUT. out_prefix.y ("coef cell", the cycle on H), out_prefix.u ("coef cell"), out_prefix.res and
//          out_prefix.lam (the residual system and its solution), and one JSON line with sizes.
//  EXAMPLE (Section "the cycle on the triangle with a leaf" of the data README):
//          SOLVER=$PWD/solve_residual.py ./exactsolve 2 6 xA.txt solA < gp.txt   (about 10 minutes, 6 GB)
// =====================================================================================================
#include "swiat_cells.h"
#include "elim_z.h"
#include "elim_modq.h"

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: exactsolve deg h target out_prefix < graph\n"); return 1; }
    int deg = atoi(argv[1]), h = atoi(argv[2]);
    string tfn = argv[3], outp = argv[4];
    if (!read_graph()) return 1;
    for (int v = 0; v < V; v++) {
        bool hasH = false;
        for (int x : vhalf[v]) if (x / 2 < h) hasH = true;
        if (hasH && base[v] / 2 >= h) { fprintf(stderr, "vertex %d of H has its base outside H; list the edges of H first\n", v); return 1; }
    }
    auto inH = [&](const Cell &c) {
        for (int e = h; e < E; e++) if (!c.elist[e].empty()) return false;
        for (int v = 0; v < V; v++) if (c.vocc[v] >= 0 && (c.vh[v] / 2 >= h || base[v] / 2 >= h)) return false;
        return true;
    };
    vector<vector<Key>> C(deg + 2);
    for (int k = 0; k <= deg + 1; k++) { enumerate(k, C[k]); sort(C[k].begin(), C[k].end()); fprintf(stderr, "C_%d: %zu\n", k, C[k].size()); }
    auto idx = [&](int k, Key key) { int j = lookup(C[k], key); if (j < 0) { fprintf(stderr, "missing cell\n"); exit(2); } return j; };
    vector<pair<Key, int>> bd;
    vector<char> inHc(C[deg].size());
    for (size_t j = 0; j < C[deg].size(); j++) inHc[j] = inH(decode(C[deg][j]));

    // ---- step 1: dropping.  For k = 1 .. deg, eliminate d_k on (k-cells outside H) x ((k-1)-cells outside H and not
    // dropped at the previous step) modulo the prime 1000003, and drop its pivot rows.  At k = deg the pivot rows S
    // carry a square minor (on some (deg-1)-cells outside H) that is nonzero mod p, hence nonzero over Z, which is all
    // the dropping argument needs (see the header).  The earlier steps only make that elimination cleaner.
    vector<char> dropS(C[deg].size(), 0);
    {
        Q = PP = 1000003;
        vector<char> prevdrop;
        for (int k = 1; k <= deg; k++) {
            vector<int> lo(C[k - 1].size(), -1);
            int nlo = 0;
            for (size_t j = 0; j < C[k - 1].size(); j++) if (!inH(decode(C[k - 1][j])) && (prevdrop.empty() || !prevdrop[j])) lo[j] = nlo++;
            QElim M(C[k].size(), nlo, false);
            for (size_t i = 0; i < C[k].size(); i++) {
                Cell ck = decode(C[k][i]);
                if (inH(ck)) continue;
                boundary(ck, bd);
                for (auto &x : bd) { int j = lo[idx(k - 1, x.first)]; if (j >= 0) M.rows[i].push_back({j, modq_ll(x.second)}); }
            }
            M.finalize();
            M.run();
            prevdrop.assign(C[k].size(), 0);
            for (size_t i = 0; i < C[k].size(); i++) if (M.rowdead[i]) prevdrop[i] = 1;
            fprintf(stderr, "d_%d outside H: %lld pivots\n", k, M.pivots);
        }
        dropS = prevdrop;
        fprintf(stderr, "dropped %lld %d-cells outside H\n", (long long)count(dropS.begin(), dropS.end(), 1), deg);
    }
    for (int k = 0; k < deg - 1; k++) vector<Key>().swap(C[k]);
    vector<Key>().swap(C[deg - 1]);                // the (deg-1)-cells are no longer needed
    vector<int> colof(C[deg].size(), -1), colcell;
    for (size_t j = 0; j < C[deg].size(); j++) if (!inHc[j] && !dropS[j]) { colof[j] = colcell.size(); colcell.push_back(j); }

    // ---- step 2: exact elimination of D
    int K1 = deg + 1;
    ZElim D(C[K1].size(), colcell.size());
    for (size_t i = 0; i < C[K1].size(); i++) {
        boundary(decode(C[K1][i]), bd);
        for (auto &x : bd) { int j = colof[idx(deg, x.first)]; if (j >= 0) D.rows[i].push_back({j, (ent)x.second}); }
    }
    D.finalize();
    fprintf(stderr, "matrix %d x %d\n", D.nr, D.nc);
    D.run();
    fprintf(stderr, "+-1 pivots %zu, ops %zu\n", D.piv.size(), D.ops.size());

    // ---- step 3: reduce the target
    vector<I> xfull(C[deg].size(), 0);
    {
        ifstream in(tfn); string line; Cell c;
        while (getline(in, line)) {
            if (!is_coef_line(line)) continue;
            size_t sp = line.find(' ');
            I coef = s2i(line.substr(0, sp));
            if (!parse_cell(line.substr(sp + 1), c)) { fprintf(stderr, "bad cell\n"); return 2; }
            xfull[idx(deg, encode(c))] += coef;
        }
    }
    vector<I> x(D.nc, 0), g(D.nr, 0);
    for (size_t j = 0; j < C[deg].size(); j++) if (colof[j] >= 0) x[colof[j]] = xfull[j];
    D.reduce(x, g);
    vector<int> rr, rc, cpos(D.nc, -1);
    for (int r = 0; r < D.nr; r++)
        if (!D.rowdead[r] && !D.rows[r].empty()) { rr.push_back(r); for (auto &e : D.rows[r]) if (cpos[e.first] < 0) { cpos[e.first] = 0; rc.push_back(e.first); } }
    long long stray = 0;                            // target entries on live columns that no residual row touches
    for (int c = 0; c < D.nc; c++) if (x[c] && cpos[c] < 0 && !D.coldead[c]) stray++;
    sort(rc.begin(), rc.end());
    for (size_t j = 0; j < rc.size(); j++) cpos[rc[j]] = j;
    fprintf(stderr, "residual %zu x %zu, target entries outside the residual columns %lld\n", rr.size(), rc.size(), stray);
    if (stray) { fprintf(stderr, "target is not in the span\n"); return 4; }
    {
        FILE *fo = fopen((outp + ".res").c_str(), "w");
        fprintf(fo, "%zu %zu\n", rr.size(), rc.size());
        for (size_t i = 0; i < rr.size(); i++) for (auto &e : D.rows[rr[i]]) fprintf(fo, "%zu %d %d\n", i, cpos[e.first], (int)e.second);
        fprintf(fo, "X\n");
        for (size_t j = 0; j < rc.size(); j++) if (x[rc[j]]) fprintf(fo, "%zu %s\n", j, i2s(x[rc[j]]).c_str());
        fclose(fo);
    }
    for (auto &r : D.rows) vector<pair<int, ent>>().swap(r);   // free the matrix before the solver runs
    for (auto &c : D.colrows) vector<int>().swap(c);
    string cmd = "python3 " + string(getenv("SOLVER") ? getenv("SOLVER") : "solve_residual.py") + " " + outp + ".res " + outp + ".lam";
    if (system(cmd.c_str()) != 0) { fprintf(stderr, "residual solve failed\n"); return 5; }
    { ifstream in(outp + ".lam"); string a, b; while (in >> a >> b) g[rr[stoi(a)]] += s2i(b); }

    // ---- step 4: u = g L
    D.replay_back(g);

    // ---- step 5: y = x - d u, exactly, and the check that it lives on H
    vector<I> y = xfull;
    for (int i = 0; i < D.nr; i++)
        if (g[i]) { boundary(decode(C[K1][i]), bd); for (auto &b2 : bd) { int j = idx(deg, b2.first); y[j] -= g[i] * b2.second; chk(y[j]); } }
    long long off = 0, ny = 0, nu = 0;
    I mx = 0, mu = 0;
    for (size_t j = 0; j < C[deg].size(); j++) if (y[j]) { if (!inHc[j]) off++; ny++; }
    fprintf(stderr, "y: %lld cells, %lld of them outside H\n", ny, off);
    if (off) return 6;
    FILE *fy = fopen((outp + ".y").c_str(), "w");
    for (size_t j = 0; j < C[deg].size(); j++)
        if (y[j]) { fprintf(fy, "%s %s\n", i2s(y[j]).c_str(), describe(decode(C[deg][j])).c_str()); mx = max(mx, y[j] < 0 ? -y[j] : y[j]); }
    fclose(fy);
    FILE *fu = fopen((outp + ".u").c_str(), "w");
    for (int i = 0; i < D.nr; i++)
        if (g[i]) { fprintf(fu, "%s %s\n", i2s(g[i]).c_str(), describe(decode(C[K1][i])).c_str()); nu++; mu = max(mu, g[i] < 0 ? -g[i] : g[i]); }
    fclose(fu);
    printf("{\"y_cells\":%lld,\"y_max_coef\":\"%s\",\"u_cells\":%lld,\"u_max_coef\":\"%s\"}\n", ny, i2s(mx).c_str(), nu, i2s(mu).c_str());
    return 0;
}
