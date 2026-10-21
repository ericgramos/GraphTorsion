// =====================================================================================================
//  orderq.cpp  --  where does a cycle live?  Coordinates of cycles in H_deg localised at p, optionally modulo
//                  the image of a subgraph and modulo further relations
// =====================================================================================================
//  THE MODULE.  Let H be the subgraph Gamma minus a set of deleted edges, C(H) the subcomplex of cells that
//  use only edges and half-edges of H (with base half-edges inside H this is exactly the reduced complex of
//  H), and Rel an optional list of extra cycles.  The program studies
//        M = C_deg / ( B_deg + C_deg(H) + span(Rel) ).
//  For a cycle z:  [z] lies in im( H_deg(Conf_n H) -> H_deg(Conf_n Gamma) ) + span[Rel]  iff  z = 0 in M.
//  (If z = y + d u + r with y in C(H), then d y = d z - d d u - d r = 0, so y is a cycle of H.)
//  With nothing deleted the pivot rows of d_deg are dropped from the columns, as in the engine; this is valid
//  because every element of B + span(Rel) and every input cycle is a cycle.  With deletions all columns of
//  cells outside H are used.
//
//  WHAT IS COMPUTED.  Rows = boundaries of the (deg+1)-cells plus the relations, columns = the deg-cells in
//  play.  Unit-pivot elimination mod Q = p^K (elim_modq.h) is followed by the local Smith form of the residual
//  block with its column transformation V.  An input cycle is reduced against the pivot rows (fact (3) of
//  elim_modq.h: in pivot order this clears every pivot column), and its residual part x' gives coordinates
//  y = x' V.  For j below the number of nonzero diagonal entries, y_j mod p^a_j is the coordinate on the
//  summand Z/p^a_j; the order of the image of the cycle in the torsion of M is max p^(a_j - v_p(y_j)).
//  "free_nonzero" counts coordinates on zero columns.  For an input that is TORSION in M these are zero over
//  Z_(p); mod Q they can be nonzero multiples of p^(K - order) because the local Smith form mod Q does not
//  clear entries that are 0 mod Q but not 0 over Z_(p).  The torsion coordinates are not affected.
//
//  EXACTNESS.  Everything is the reduction of a Z_(p) computation.  To conclude that a torsion element is
//  zero in M from its vanishing mod p^K one needs the torsion of M to have exponent < p^K; check this by
//  comparing (unit pivots + number of diagonal entries) with the rank over a large prime (run with that prime
//  and K = 1).  This was done for every conclusion reported in the data folder.
//
//  INPUT CYCLES.  Blocks "Z id order" followed by "coef cell" lines (the generator files of certq; other tags
//  such as W and F are skipped).  Cells may be written for a prefix of the edge list.  They are taken to be
//  written in the default basis (lowest half-edge as base) and are rewritten in the basis used here with
//  h - b = (h - b') - (b - b') at every occupied vertex whose base was overridden.
//
//  USAGE.  orderq deg p K cycles_file deleted_edges base_overrides [relations_file] < graph_with_n
//          deleted_edges: comma separated edge indices, or "-"; base_overrides: "v:h,v:h", or "-" (needed when
//          the lowest half-edge of a vertex of H lies on a deleted edge).
//  OUTPUT. one JSON line {"diag": [a_j], "cycles": [{id, order_given, cells, ignored_cells, order_in_M,
//          coords, free_nonzero, pivot_nonzero}]}.  ignored_cells counts cells of the input that lie in C(H)
//          or among the dropped cells.  pivot_nonzero must be 0 (it is a self-check).
// =====================================================================================================
#include "swiat_cells.h"
#include "elim_modq.h"

int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage: orderq deg p K cycles deleted base [relations] < graph\n"); return 1; }
    int deg = atoi(argv[1]);
    PP = strtoull(argv[2], 0, 10);
    int KK = atoi(argv[3]);
    Q = ppow(KK);
    if (Q >= (1ULL << 31)) { fprintf(stderr, "modulus too large\n"); return 1; }
    string cycfn = argv[4], dels = argv[5], bases = argv[6];
    if (!read_graph()) return 1;
    vector<int> base0 = base;                       // basis in which the input cycles are written
    vector<char> deleted(E, 0);
    bool anydel = false;
    if (dels != "-") { stringstream ss(dels); string t; while (getline(ss, t, ',')) { deleted[atoi(t.c_str())] = 1; anydel = true; } }
    if (bases != "-") {
        stringstream ss(bases); string t;
        while (getline(ss, t, ',')) {
            int v, h; sscanf(t.c_str(), "%d:%d", &v, &h);
            if (hvert[h] != v) { fprintf(stderr, "bad base override\n"); return 1; }
            base[v] = h;
        }
    }
    for (int v = 0; v < V; v++) {
        bool hasH = false;
        for (int h : vhalf[v]) if (!deleted[h / 2]) hasH = true;
        if (hasH && deleted[base[v] / 2]) { fprintf(stderr, "base of vertex %d lies on a deleted edge; override it\n", v); return 1; }
    }
    auto inH = [&](const Cell &c) {
        for (int e = 0; e < E; e++) if (deleted[e] && !c.elist[e].empty()) return false;
        for (int v = 0; v < V; v++) if (c.vocc[v] >= 0 && (deleted[c.vh[v] / 2] || deleted[base[v] / 2])) return false;
        return true;
    };

    int lo = anydel ? deg : 0;                      // lower degrees are only needed for dropping
    vector<vector<Key>> C(deg + 2);
    for (int k = lo; k <= deg + 1; k++) { enumerate(k, C[k]); sort(C[k].begin(), C[k].end()); fprintf(stderr, "C_%d: %zu\n", k, C[k].size()); }
    vector<pair<Key, int>> bd;

    // ---- columns: kept deg-cells (no deletion) or deg-cells outside H (deletion)
    vector<int> colof(C[deg].size(), -1), colcell;
    if (!anydel) {
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
        for (size_t j = 0; j < C[deg].size(); j++) if (!drop[j]) { colof[j] = colcell.size(); colcell.push_back(j); }
    } else {
        for (size_t j = 0; j < C[deg].size(); j++) if (!inH(decode(C[deg][j]))) { colof[j] = colcell.size(); colcell.push_back(j); }
    }

    // ---- reading chains in the basis used here
    // Adds coef * (cell rewritten in the current basis) to the column vector `add(col, value)`.
    auto add_chain_line = [&](const string &line, const function<void(int, uint32_t)> &add, long long &ignored) {
        size_t sp = line.find(' ');
        long long coef = atoll(line.substr(0, sp).c_str());
        Cell c;
        if (!parse_cell(line.substr(sp + 1), c)) { fprintf(stderr, "bad cell %s\n", line.c_str()); exit(2); }
        vector<int> occ;
        for (int v = 0; v < V; v++) if (c.vocc[v] >= 0) occ.push_back(v);
        // product over occupied vertices of (h - b') - (b0 - b'); a term with h = b' vanishes
        function<void(int, long long, Cell &)> rec = [&](int i, long long sg, Cell &cc) {
            if (i == (int)occ.size()) {
                int j = lookup(C[deg], encode(cc));
                if (j < 0) { fprintf(stderr, "cell not in complex\n"); exit(2); }
                if (colof[j] >= 0) add(colof[j], modq_ll(sg * coef)); else ignored++;
                return;
            }
            int v = occ[i], h = c.vh[v];
            if (h != base[v]) { cc.vh[v] = h; rec(i + 1, sg, cc); }
            if (base0[v] != base[v]) { cc.vh[v] = base0[v]; rec(i + 1, -sg, cc); }
            cc.vh[v] = h;
        };
        Cell cc = c;
        rec(0, 1, cc);
    };

    // ---- rows: boundaries of (deg+1)-cells, then relations
    int K1 = deg + 1;
    QElim D(C[K1].size(), colcell.size(), true);
    for (size_t i = 0; i < C[K1].size(); i++) {
        boundary(decode(C[K1][i]), bd);
        for (auto &x : bd) { int j = colof[lookup(C[deg], x.first)]; if (j >= 0) D.rows[i].push_back({j, modq_ll(x.second)}); }
    }
    if (argc > 7) {
        ifstream in(argv[7]); string line; vector<pair<int, uint32_t>> row; long long nrel = 0, ign = 0;
        auto flushrow = [&]() { if (!row.empty()) { D.add_row(row); nrel++; row.clear(); } };
        while (getline(in, line)) {
            if (line.rfind("Z ", 0) == 0) { flushrow(); continue; }
            if (!is_coef_line(line)) continue;
            add_chain_line(line, [&](int col, uint32_t v) { row.push_back({col, v}); }, ign);
        }
        flushrow();
        fprintf(stderr, "added %lld relation rows\n", nrel);
    }
    D.finalize();
    fprintf(stderr, "M: %d x %d\n", D.nr, D.nc);
    D.run();
    vector<Op>().swap(D.ops);                       // the log is not needed here, only the snapshots

    // ---- residual and its local Smith form (column transformation only)
    vector<int> rr, rc, cpos;
    D.residual(rr, rc, cpos);
    int nr = rr.size(), nc = rc.size();
    long long livezero = 0;
    for (int c = 0; c < D.nc; c++) if (!D.coldead[c] && cpos[c] < 0) livezero++;
    fprintf(stderr, "unit rank %lld, residual %d x %d, zero columns %lld\n", D.pivots, nr, nc, livezero);
    vector<vector<u64>> R(nr, vector<u64>(nc, 0)), U, Vm;
    for (int i = 0; i < nr; i++) for (auto &x : D.rows[rr[i]]) R[i][cpos[x.first]] = x.second;
    vector<int> dv = local_smith(R, nr, nc, U, Vm, false, KK);
    map<int, int> hist; for (int v : dv) hist[v]++;
    fprintf(stderr, "torsion of M below p^K (valuation:count):");
    for (auto &h : hist) fprintf(stderr, " %d:%d", h.first, h.second);
    fprintf(stderr, "\n");
    vector<u64> pinv(D.snapshots().size());                 // inverse of each pivot entry, computed once
    for (size_t o = 0; o < D.snapshots().size(); o++) { int pc = D.snapshots()[o].first; for (auto &x : D.snapshots()[o].second) if (x.first == pc) pinv[o] = egcd_inv(x.second); }

    // ---- the input cycles, one block at a time
    printf("{\"diag\":[");
    for (size_t i = 0; i < dv.size(); i++) printf("%s%d", i ? "," : "", dv[i]);
    printf("],\"cycles\":[");
    ifstream in(cycfn); string line; bool first = true;
    vector<u64> z(D.nc, 0);
    string id = ""; long long ord = 0, ncells = 0, ignored = 0;
    auto finish = [&]() {
        if (id.empty()) return;
        for (size_t o = 0; o < D.snapshots().size(); o++) {            // reduce against the pivot rows in order
            int pc = D.snapshots()[o].first; u64 a = z[pc];
            if (!a) continue;
            u64 m = mulq(a, pinv[o]);
            for (auto &x : D.snapshots()[o].second) z[x.first] = (z[x.first] + Q - mulq(m, x.second)) % Q;
        }
        long long badpiv = 0, freez = 0;
        for (int c = 0; c < D.nc; c++) if (z[c]) { if (D.coldead[c]) badpiv++; else if (cpos[c] < 0) freez++; }
        vector<u64> y(nc, 0);                                  // coordinates y = x' V
        for (int j = 0; j < nc; j++) if (z[rc[j]]) { u64 a = z[rc[j]]; for (int i = 0; i < nc; i++) if (Vm[j][i]) y[i] = (y[i] + mulq(a, Vm[j][i])) % Q; }
        int ordexp = 0; string coords = "";
        for (int i = 0; i < nc; i++) {
            if (i < (int)dv.size()) {
                u64 m = ppow(dv[i]), t = y[i] % m;
                coords += (i ? "," : "") + to_string(t);
                if (t) ordexp = max(ordexp, dv[i] - val(t));
            } else if (y[i] % Q) freez++;
        }
        printf("%s{\"id\":\"%s\",\"order_given\":%lld,\"cells\":%lld,\"ignored_cells\":%lld,\"order_in_M\":%llu,\"coords\":[%s],\"free_nonzero\":%lld,\"pivot_nonzero\":%lld}",
               first ? "" : ",", id.c_str(), ord, ncells, ignored, (unsigned long long)ppow(ordexp), coords.c_str(), freez, badpiv);
        fflush(stdout);
        first = false;
        fill(z.begin(), z.end(), 0);
        id = ""; ncells = ignored = 0;
    };
    while (getline(in, line)) {
        if (line.rfind("Z ", 0) == 0) { finish(); stringstream ss(line.substr(2)); ss >> id >> ord; continue; }
        if (id.empty() || !is_coef_line(line)) continue;
        ncells++;
        add_chain_line(line, [&](int col, uint32_t v) { z[col] = (z[col] + v) % Q; }, ignored);
    }
    finish();
    printf("]}\n");
    return 0;
}
