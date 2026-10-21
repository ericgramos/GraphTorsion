// =====================================================================================================
//  spanq.cpp  --  fast span queries on the output of reducemod
// =====================================================================================================
//  reducemod has written, for every cycle of interest, its reduced vector in (Z/Q)^{live}, together with the
//  residual rows R; then C_deg / (B_deg + span S) = (Z/Q)^{live} / (span R + span of the reduced S).
//  Blocks whose header starts with "gen" are the generators z_t of the torsion T (header "gen id order");
//  the others ("tor ...") are candidate relations, numbered 0, 1, 2, ... in file order.
//
//  The program reads the file once and then answers queries from stdin, one per line.
//    span i j k ...      Rows = R + the listed relation blocks.  Prints the orderq-style JSON for the generators
//                        (diag of the local Smith form, coordinates of each z_t), so that kernel.py gives
//                        T cap (B + span S): the torsion elements that are combinations of the listed blocks.
//  CAUTION.  Everything is mod Q.  A torsion element that vanishes mod Q in C/(B + span S) is zero there only if
//  the torsion of that module has exponent below Q.  When the images of the listed cycles are divisible by high
//  powers of p this fails, and the answer can be a false "yes" (with Q = 2^8 the search of minimal_tori.py was
//  fooled this way).  Use a large K (up to 30 for p = 2) and confirm every final answer exactly, as
//  data/survey/k4_family/support/tori/exact_check.py does (rational kernel of the images and the cocycle).
//  USAGE.  spanq reduced_file < queries
// =====================================================================================================
#include "elim_modq.h"
#include <sstream>

typedef vector<pair<int, uint32_t>> SV;
// Parse a line of 'col:val' pairs into a sparse vector.
SV parse_vec(const string &line, size_t from) {
    SV v; stringstream ss(line.substr(from)); string tok;
    while (ss >> tok) { size_t c = tok.find(':'); v.push_back({stoi(tok.substr(0, c)), (uint32_t)stoull(tok.substr(c + 1))}); }
    return v;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: spanq reduced_file < queries\n"); return 1; }
    ifstream in(argv[1]);
    string line;
    int L = 0;
    vector<SV> R, rel, gen;
    vector<long long> genord;
    vector<string> relhdr;
    getline(in, line);
    { stringstream ss(line); string t; unsigned long long q; ss >> t >> L >> q; Q = q; }
    // the prime: smallest divisor of Q
    for (u64 d = 2; d <= Q; d++) if (Q % d == 0) { PP = d; break; }
    int KK = 0; for (u64 t = Q; t > 1; t /= PP) KK++;
    string hdr;
    while (getline(in, line)) {
        if (line.rfind("R", 0) == 0 && (line.size() == 1 || line[1] == ' ')) { R.push_back(parse_vec(line, 1)); continue; }
        if (line.rfind("B ", 0) == 0) { hdr = line.substr(2); continue; }
        SV v = parse_vec(line, 0);
        if (hdr.rfind("gen", 0) == 0) { stringstream ss(hdr); string g; long long id, od; ss >> g >> id >> od; gen.push_back(v); genord.push_back(od); }
        else { rel.push_back(v); relhdr.push_back(hdr); }
    }
    fprintf(stderr, "live %d, residual rows %zu, generators %zu, relation blocks %zu\n", L, R.size(), gen.size(), rel.size());
    fprintf(stdout, "READY\n"); fflush(stdout);

    while (getline(cin, line)) {
        stringstream ss(line); string cmd; ss >> cmd;
        vector<int> ids;
        if (cmd != "span") { printf("null\n"); fflush(stdout); continue; }
        { int i; while (ss >> i) ids.push_back(i); }
        // rows: R first, then the selected relations.  The pivot snapshots are needed to reduce the generators,
        // so the eliminator must keep them (keep = true); the operation log is discarded afterwards.
        QElim M(R.size() + ids.size(), L, true);
        for (size_t i = 0; i < R.size(); i++) M.rows[i] = R[i];
        for (size_t k = 0; k < ids.size(); k++) M.rows[R.size() + k] = rel[ids[k]];
        M.finalize(); M.run();
        vector<Op>().swap(M.ops);
        vector<int> rr, rc, cpos;
        M.residual(rr, rc, cpos);
        int nr = rr.size(), nc = rc.size();
        vector<vector<u64>> A(nr, vector<u64>(nc, 0)), U, Vm;
        for (int i = 0; i < nr; i++) for (auto &x : M.rows[rr[i]]) A[i][cpos[x.first]] = x.second;
        vector<int> dv = local_smith(A, nr, nc, U, Vm, false, KK);
        vector<u64> pinv(M.snapshots().size());
        for (size_t o = 0; o < M.snapshots().size(); o++) { int pc = M.snapshots()[o].first; for (auto &x : M.snapshots()[o].second) if (x.first == pc) pinv[o] = egcd_inv(x.second); }
        // reduce a dense vector against the snapshots in pivot order (fact (3) of elim_modq.h)
        auto reduce = [&](vector<u64> &z) {
            for (size_t o = 0; o < M.snapshots().size(); o++) {
                int pc = M.snapshots()[o].first; u64 a = z[pc];
                if (!a) continue;
                u64 m = mulq(a, pinv[o]);
                for (auto &x : M.snapshots()[o].second) z[x.first] = (z[x.first] + Q - mulq(m, x.second)) % Q;
            }
        };
        {
            printf("{\"diag\":[");
            for (size_t i = 0; i < dv.size(); i++) printf("%s%d", i ? "," : "", dv[i]);
            printf("],\"cycles\":[");
            for (size_t t = 0; t < gen.size(); t++) {
                vector<u64> z(L, 0);
                for (auto &x : gen[t]) z[x.first] = x.second;
                reduce(z);
                long long freez = 0;
                for (int c = 0; c < L; c++) if (z[c] && !M.coldead[c] && cpos[c] < 0) freez++;
                vector<u64> y(nc, 0);
                for (int j = 0; j < nc; j++) if (z[rc[j]]) for (int i = 0; i < nc; i++) if (Vm[j][i]) y[i] = (y[i] + mulq(z[rc[j]], Vm[j][i])) % Q;
                int ordexp = 0; string coords;
                for (int i = 0; i < nc; i++) {
                    if (i < (int)dv.size()) { u64 m = ppow(dv[i]), tt = y[i] % m; coords += (i ? "," : "") + to_string(tt); if (tt) ordexp = max(ordexp, dv[i] - val(tt)); }
                    else if (y[i] % Q) freez++;
                }
                printf("%s{\"id\":\"%zu\",\"order_given\":%lld,\"order_in_M\":%llu,\"coords\":[%s],\"free_nonzero\":%lld}", t ? "," : "", t, genord[t], (unsigned long long)ppow(ordexp), coords.c_str(), freez);
            }
            printf("]}\n");
            fflush(stdout);
        }
    }
    return 0;
}
