// =====================================================================================================
//  cert.cpp  --  finder for mod-2 torsion certificates (see certp.cpp for any prime; same algorithm)
// =====================================================================================================
//  Builds the reduced ordered Swiatkowski complex directly (its own copy of the engine's cell code), then
//  searches for a pair (w, f) proving that a class of order exactly 2 exists in H_deg; see certp.cpp.
//  USAGE.  cert deg seed outfile < graph_with_n   (input "V E n" + edges, as for the engine)
// =====================================================================================================

#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

int V, E, N;
bool checksq = false;
vector<int> pcolor;
bool colored = false;
vector<int> etail, ehead;
vector<int> hvert;         // half-edge -> vertex; half-edge 2e = tail end, 2e+1 = head end
vector<vector<int>> vhalf; // vertex -> half-edges
vector<int> base;          // vertex -> base half-edge (-1 if no half-edges)

// cell encoding: 2N bytes; for particle p: loc (half-edge id if at vertex, 2E+e if on edge), rank
struct Cell {
    vector<int> vocc;          // vertex -> particle or -1
    vector<int> vh;            // vertex -> half-edge used
    vector<vector<int>> elist; // edge -> particles tail..head
};

typedef unsigned __int128 Key;
Key encode(const Cell &c) {
    Key s = 0;
    for (int v = 0; v < V; v++)
        if (c.vocc[v] >= 0) {
            int p = c.vocc[v];
            s |= (Key)(unsigned)(c.vh[v]) << (12 * p);
        }
    for (int e = 0; e < E; e++)
        for (int r = 0; r < (int)c.elist[e].size(); r++) {
            int p = c.elist[e][r];
            s |= (Key)(unsigned)((2 * E + e) | (r << 8)) << (12 * p);
        }
    return s;
}

void canonicalize(Cell &c) {
    if (!colored)
        return;
    // pools of labels per colour, in increasing order
    map<int, vector<int>> pool;
    for (int p = 0; p < N; p++)
        pool[pcolor[p]].push_back(p);
    map<int, int> used;
    auto relabel = [&](int p) {
        int col = pcolor[p];
        return pool[col][used[col]++];
    };
    for (int v = 0; v < V; v++)
        if (c.vocc[v] >= 0)
            c.vocc[v] = relabel(c.vocc[v]);
    for (int e = 0; e < E; e++)
        for (auto &p : c.elist[e])
            p = relabel(p);
}
Cell decode(Key s) {
    Cell c;
    c.vocc.assign(V, -1);
    c.vh.assign(V, -1);
    c.elist.assign(E, {});
    vector<vector<pair<int, int>>> tmp(E);
    for (int p = 0; p < N; p++) {
        unsigned w = (unsigned)((s >> (12 * p)) & 0xFFF);
        int loc = w & 0xFF, r = w >> 8;
        if (loc < 2 * E) {
            int v = hvert[loc];
            c.vocc[v] = p;
            c.vh[v] = loc;
        } else
            tmp[loc - 2 * E].push_back({r, p});
    }
    for (int e = 0; e < E; e++) {
        sort(tmp[e].begin(), tmp[e].end());
        for (auto &x : tmp[e])
            c.elist[e].push_back(x.second);
    }
    return c;
}

// enumerate all cells of degree k
void enumerate(int k, vector<Key> &out) {
    Cell c;
    c.vocc.assign(V, -1);
    c.vh.assign(V, -1);
    c.elist.assign(E, {});
    int moving = 0;
    function<void(int)> rec = [&](int p) {
        if (p == N) {
            if (moving == k) {
                if (colored) {
                    Cell d = c;
                    canonicalize(d);
                    out.push_back(encode(d));
                } else
                    out.push_back(encode(c));
            }
            return;
        }
        int remaining = N - p;
        // place at a vertex
        if (moving < k) {
            for (int v = 0; v < V; v++)
                if (c.vocc[v] < 0) {
                    for (int h : vhalf[v])
                        if (h != base[v]) {
                            c.vocc[v] = p;
                            c.vh[v] = h;
                            moving++;
                            rec(p + 1);
                            moving--;
                            c.vocc[v] = -1;
                            c.vh[v] = -1;
                        }
                }
        }
        // place on an edge (only if we can still reach k moving)
        if (k - moving <= remaining - 1) {
            for (int e = 0; e < E; e++) {
                int L = c.elist[e].size();
                for (int pos = 0; pos <= L; pos++) {
                    c.elist[e].insert(c.elist[e].begin() + pos, p);
                    rec(p + 1);
                    c.elist[e].erase(c.elist[e].begin() + pos);
                }
            }
        }
    };
    rec(0);
}

// ================= F_2 elimination with stored pivots (for kernels) =================
struct F2Elim {
    int nr, nc;
    vector<vector<int>> rows;
    vector<vector<int>> colrows;
    vector<int> colcnt;
    vector<char> rowdead, coldead;
    vector<pair<int, vector<int>>> piv; // (pivot column, pivot row at pivot time)
    bool keep;
    F2Elim(int r, int c, bool keep_)
        : nr(r), nc(c), rows(r), colrows(c), colcnt(c, 0), rowdead(r, 0), coldead(c, 0), keep(keep_) {}
    void finalize() {
        for (int r = 0; r < nr; r++) {
            auto &R = rows[r];
            sort(R.begin(), R.end());
            vector<int> M;
            for (size_t i = 0; i < R.size();) {
                size_t j = i;
                while (j < R.size() && R[j] == R[i])
                    j++;
                if ((j - i) & 1)
                    M.push_back(R[i]);
                i = j;
            }
            R.swap(M);
            for (int c : R) {
                colrows[c].push_back(r);
                colcnt[c]++;
            }
        }
    }
    bool has(int r, int c) const {
        auto &R = rows[r];
        return binary_search(R.begin(), R.end(), c);
    }
    void pivot(int r, int c) {
        vector<int> &R = rows[r];
        vector<int> targets;
        for (int r2 : colrows[c])
            if (r2 != r && !rowdead[r2])
                targets.push_back(r2);
        sort(targets.begin(), targets.end());
        targets.erase(unique(targets.begin(), targets.end()), targets.end());
        for (int r2 : targets) {
            if (!has(r2, c))
                continue;
            auto &S = rows[r2];
            vector<int> out;
            out.reserve(S.size() + R.size());
            size_t i = 0, j = 0;
            while (i < S.size() || j < R.size()) {
                if (j == R.size() || (i < S.size() && S[i] < R[j])) {
                    out.push_back(S[i]);
                    i++;
                } else if (i == S.size() || R[j] < S[i]) {
                    out.push_back(R[j]);
                    colrows[R[j]].push_back(r2);
                    colcnt[R[j]]++;
                    j++;
                } else {
                    colcnt[S[i]]--;
                    i++;
                    j++;
                }
            }
            S.swap(out);
        }
        for (int x : R)
            colcnt[x]--;
        rowdead[r] = 1;
        coldead[c] = 1;
        colcnt[c] = 0;
        if (keep)
            piv.push_back({c, std::move(rows[r])});
        vector<int>().swap(rows[r]);
        vector<int>().swap(colrows[c]);
    }
    void run() {
        bool any = true;
        while (any) {
            any = false;
            typedef pair<int, int> PI;
            priority_queue<PI, vector<PI>, greater<PI>> pq;
            for (int c = 0; c < nc; c++)
                if (!coldead[c] && colcnt[c] > 0)
                    pq.push({colcnt[c], c});
            while (!pq.empty()) {
                auto [cnt, c] = pq.top();
                pq.pop();
                if (coldead[c] || colcnt[c] != cnt || cnt == 0)
                    continue;
                if (colrows[c].size() > 2 * (size_t)colcnt[c] + 8) {
                    vector<int> M;
                    for (int r : colrows[c])
                        if (!rowdead[r] && has(r, c))
                            M.push_back(r);
                    sort(M.begin(), M.end());
                    M.erase(unique(M.begin(), M.end()), M.end());
                    colrows[c].swap(M);
                }
                int best = -1;
                size_t bl = SIZE_MAX;
                for (int r : colrows[c])
                    if (!rowdead[r] && rows[r].size() < bl && has(r, c)) {
                        bl = rows[r].size();
                        best = r;
                    }
                if (best < 0)
                    continue;
                vector<int> touched(rows[best].begin(), rows[best].end());
                pivot(best, c);
                any = true;
                for (int c2 : touched)
                    if (!coldead[c2] && colcnt[c2] > 0)
                        pq.push({colcnt[c2], c2});
            }
        }
    }
    // random vector in the right kernel {x : rows . x = 0}
    vector<char> kernel_vector(mt19937_64 &rng) {
        vector<char> x(nc, 0);
        for (int c = 0; c < nc; c++)
            if (!coldead[c])
                x[c] = rng() & 1;
        for (int i = (int)piv.size() - 1; i >= 0; i--) {
            int c = piv[i].first;
            char s = 0;
            for (int j : piv[i].second)
                if (j != c)
                    s ^= x[j];
            x[c] = s;
        }
        return x;
    }
};

string describe(const Cell &c) {
    string s = "V[";
    for (int v = 0; v < V; v++) {
        if (v)
            s += ",";
        if (c.vocc[v] < 0)
            s += "-";
        else
            s += to_string(c.vocc[v]) + "/" + to_string(c.vh[v]);
    }
    s += "] E[";
    for (int e = 0; e < E; e++) {
        if (e)
            s += "|";
        for (size_t r = 0; r < c.elist[e].size(); r++) {
            if (r)
                s += " ";
            s += to_string(c.elist[e][r]);
        }
    }
    return s + "]";
}

// signed boundary of cell (degree k) as (face key, coefficient)
void boundary(const Cell &c, vector<pair<Key, int>> &out) {
    out.clear();
    int idx = 0;
    for (int v = 0; v < V; v++)
        if (c.vocc[v] >= 0) {
            int p = c.vocc[v], h = c.vh[v], h0 = base[v];
            int sgn = (idx % 2 == 0) ? 1 : -1;
            idx++;
            for (int t = 0; t < 2; t++) {
                int hh = t == 0 ? h : h0;
                Cell f = c;
                f.vocc[v] = -1;
                f.vh[v] = -1;
                int e = hh / 2;
                if (hh % 2 == 0)
                    f.elist[e].insert(f.elist[e].begin(), p);
                else
                    f.elist[e].push_back(p);
                out.push_back({encode(f), t == 0 ? sgn : -sgn});
            }
        }
}

int main(int argc, char **argv) {
    int deg = atoi(argv[1]); // torsion degree i: certificate for 2-torsion in H_i
    unsigned long long seed = argc > 2 ? strtoull(argv[2], 0, 10) : 1;
    string outfn = argc > 3 ? argv[3] : "cert.txt";
    if (scanf("%d %d %d", &V, &E, &N) != 3)
        return 1;
    etail.resize(E);
    ehead.resize(E);
    hvert.resize(2 * E);
    vhalf.assign(V, {});
    for (int e = 0; e < E; e++) {
        if (scanf("%d %d", &etail[e], &ehead[e]) != 2)
            return 1;
        hvert[2 * e] = etail[e];
        hvert[2 * e + 1] = ehead[e];
        vhalf[etail[e]].push_back(2 * e);
        vhalf[ehead[e]].push_back(2 * e + 1);
    }
    base.assign(V, -1);
    for (int v = 0; v < V; v++)
        if (!vhalf[v].empty())
            base[v] = vhalf[v][0];
    mt19937_64 rng(seed);
    // cells of degree 0..deg+1
    vector<vector<Key>> C(deg + 2);
    for (int k = 0; k <= deg + 1; k++) {
        enumerate(k, C[k]);
        sort(C[k].begin(), C[k].end());
        fprintf(stderr, "C_%d: %zu\n", k, C[k].size());
    }
    auto idx = [&](int k, Key key) {
        auto it = lower_bound(C[k].begin(), C[k].end(), key);
        if (it == C[k].end() || *it != key) {
            fprintf(stderr, "missing face\n");
            exit(2);
        }
        return (int)(it - C[k].begin());
    };
    vector<pair<Key, int>> bd;
    // F_2 reductions of d_1..d_deg with column dropping, to find dropped deg-cells R_deg
    vector<char> drop; // dropped (k-1)-cells
    for (int k = 1; k <= deg; k++) {
        F2Elim M(C[k].size(), C[k - 1].size(), false);
        for (size_t i = 0; i < C[k].size(); i++) {
            boundary(decode(C[k][i]), bd);
            for (auto &x : bd) {
                int j = idx(k - 1, x.first);
                if (drop.empty() || !drop[j])
                    M.rows[i].push_back(j);
            }
        }
        M.finalize();
        M.run();
        drop.assign(C[k].size(), 0);
        long long nd = 0;
        for (size_t i = 0; i < C[k].size(); i++)
            if (M.rowdead[i]) {
                drop[i] = 1;
                nd++;
            }
        fprintf(stderr, "d_%d over F2: rank %lld\n", k, nd);
    }
    // M' = d_{deg+1}^T with dropped columns removed (rows: (deg+1)-cells, cols: kept deg-cells)
    int K1 = deg + 1;
    vector<int> keepidx(C[deg].size(), -1);
    vector<int> kept;
    for (size_t j = 0; j < C[deg].size(); j++)
        if (!drop[j]) {
            keepidx[j] = kept.size();
            kept.push_back(j);
        }
    vector<vector<pair<int, int>>> full(C[K1].size()); // full integer boundary of each (deg+1)-cell
    F2Elim Mf(C[K1].size(), kept.size(), true);
    vector<vector<int>> colsT(kept.size());
    for (size_t i = 0; i < C[K1].size(); i++) {
        boundary(decode(C[K1][i]), bd);
        for (auto &x : bd) {
            int j = idx(deg, x.first);
            full[i].push_back({j, x.second});
            if (keepidx[j] >= 0) {
                Mf.rows[i].push_back(keepidx[j]);
                colsT[keepidx[j]].push_back(i);
            }
        }
    }
    Mf.finalize();
    F2Elim Mw(kept.size(), C[K1].size(), true);
    for (size_t j = 0; j < kept.size(); j++)
        Mw.rows[j] = colsT[j];
    vector<vector<int>>().swap(colsT);
    Mw.finalize();
    fprintf(stderr, "eliminating cocycle system\n");
    Mf.run();
    fprintf(stderr, "eliminating cycle system\n");
    Mw.run();
    long long kerw = 0;
    for (int c = 0; c < Mw.nc; c++)
        if (!Mw.coldead[c])
            kerw++;
    fprintf(stderr, "dim ker d_%d mod 2 = %lld\n", K1, kerw);
    for (int tw = 0; tw < 40; tw++) {
        vector<char> w = Mw.kernel_vector(rng);
        // integral boundary of w (as 0/1 chain)
        unordered_map<int, long long> dw;
        for (size_t i = 0; i < w.size(); i++)
            if (w[i])
                for (auto &x : full[i])
                    dw[x.first] += x.second;
        bool even = true;
        long long nz = 0;
        for (auto &kv : dw) {
            if (kv.second % 2)
                even = false;
            if (kv.second)
                nz++;
        }
        if (!even) {
            fprintf(stderr, "BUG: boundary of mod-2 cycle not even\n");
            return 3;
        }
        if (nz == 0)
            continue; // w lifts to an integral cycle with zero boundary: no information
        for (int tf = 0; tf < 20; tf++) {
            vector<char> fk = Mf.kernel_vector(rng);
            long long pairing = 0;
            for (auto &kv : dw)
                if (keepidx[kv.first] >= 0 && fk[keepidx[kv.first]])
                    pairing += kv.second / 2;
            if (((pairing % 2) + 2) % 2 == 1) {
                FILE *fo = fopen(outfn.c_str(), "w");
                fprintf(fo, "GRAPH %d %d %d\n", V, E, N);
                for (int e = 0; e < E; e++)
                    fprintf(fo, "EDGE %d %d\n", etail[e], ehead[e]);
                fprintf(fo, "DEGREE %d\n", deg);
                long long nw = 0, nf = 0;
                for (size_t i = 0; i < w.size(); i++)
                    if (w[i]) {
                        fprintf(fo, "W %s\n", describe(decode(C[K1][i])).c_str());
                        nw++;
                    }
                for (size_t j = 0; j < kept.size(); j++)
                    if (fk[j]) {
                        fprintf(fo, "F %s\n", describe(decode(C[deg][kept[j]])).c_str());
                        nf++;
                    }
                fclose(fo);
                printf("{\"certificate\":\"%s\",\"w_cells\":%lld,\"f_cells\":%lld,\"z_support\":%lld,\"tries_"
                       "w\":%d,\"tries_f\":%d}\n",
                       outfn.c_str(), nw, nf, nz, tw + 1, tf + 1);
                return 0;
            }
        }
    }
    printf("{\"certificate\":null}\n");
    return 1;
}
