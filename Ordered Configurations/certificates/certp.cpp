// =====================================================================================================
//  certp.cpp  --  finder for mod-p torsion certificates
// =====================================================================================================
//  A CERTIFICATE that H_deg(Conf_n Gamma; Z) contains an element of order exactly p consists of
//    * an integral (deg+1)-chain w with d w = 0 (mod p); then z = d w / p is an integral deg-cycle and
//      p z = d w is a boundary;
//    * a mod-p deg-cochain f with f(d s) = 0 (mod p) for every (deg+1)-cell s (a cocycle), and f(z) != 0.
//  If z were a boundary d u then f(z) = f(d u) = 0; so [z] is nonzero and p [z] = 0.
//
//  ALGORITHM.  Reduce d_1, ..., d_deg mod p with unit pivots, dropping pivot rows from the next map's
//  columns (as in the engine).  Then form the matrix of d_{deg+1} restricted to the kept deg-cells, and
//  its transpose.  Random kernel vectors of the transpose give candidates w (columns = (deg+1)-cells
//  with d w = 0 mod p on the kept coordinates); random kernel vectors of the matrix itself give cocycles
//  f.  For each candidate w, the integral boundary d w is formed and checked to be divisible by p
//  (exit 3 otherwise), and the pairing f(d w / p) mod p is evaluated.  The first nonzero pairing is written
//  to the output file:
//        GRAPH V E n / EDGE u v (E lines) / DEGREE deg / PRIME p /
//        W coefficient cell-description   (one line per cell of w)
//        F coefficient cell-description   (one line per cell of f)
//  Cell descriptions list, for each vertex, the particle and half-edge, and for each edge its word.
//  The independent verifier is validation/verify_cert.py.
//  USAGE.  certp deg prime seed outfile < graph_with_n
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

// ================= F_p elimination with stored pivots (for kernels) =================
typedef unsigned long long u64;
u64 PR;
u64 mulm(u64 a, u64 b) {
    return (unsigned __int128)a * b % PR;
}
u64 powm(u64 a, u64 e) {
    u64 r = 1;
    while (e) {
        if (e & 1)
            r = mulm(r, a);
        a = mulm(a, a);
        e >>= 1;
    }
    return r;
}
u64 inv(u64 a) {
    return powm(a, PR - 2);
}
struct FpElim {
    int nr, nc;
    vector<vector<pair<int, uint32_t>>> rows;
    vector<vector<int>> colrows;
    vector<int> colcnt;
    vector<char> rowdead, coldead;
    vector<pair<int, vector<pair<int, uint32_t>>>> piv;
    bool keep;
    FpElim(int r, int c, bool keep_)
        : nr(r), nc(c), rows(r), colrows(c), colcnt(c, 0), rowdead(r, 0), coldead(c, 0), keep(keep_) {}
    void finalize() {
        for (int r = 0; r < nr; r++) {
            auto &R = rows[r];
            sort(R.begin(), R.end());
            vector<pair<int, uint32_t>> M;
            for (auto &x : R) {
                if (!M.empty() && M.back().first == x.first)
                    M.back().second = (uint32_t)((M.back().second + (u64)x.second) % PR);
                else
                    M.push_back({x.first, (uint32_t)(x.second % PR)});
            }
            vector<pair<int, uint32_t>> M2;
            for (auto &x : M)
                if (x.second)
                    M2.push_back(x);
            R.swap(M2);
            for (auto &x : R) {
                colrows[x.first].push_back(r);
                colcnt[x.first]++;
            }
        }
    }
    u64 get(int r, int c) const {
        auto &R = rows[r];
        auto it = lower_bound(R.begin(), R.end(), make_pair(c, (uint32_t)0));
        return (it != R.end() && it->first == c) ? it->second : 0;
    }
    void pivot(int r, int c) {
        auto &R = rows[r];
        u64 ui = inv(get(r, c));
        vector<int> targets;
        for (int r2 : colrows[c])
            if (r2 != r && !rowdead[r2])
                targets.push_back(r2);
        sort(targets.begin(), targets.end());
        targets.erase(unique(targets.begin(), targets.end()), targets.end());
        for (int r2 : targets) {
            u64 a = get(r2, c);
            if (!a)
                continue;
            u64 f = mulm(a, ui);
            auto &S = rows[r2];
            vector<pair<int, uint32_t>> out;
            out.reserve(S.size() + R.size());
            size_t i = 0, j = 0;
            while (i < S.size() || j < R.size()) {
                if (j == R.size() || (i < S.size() && S[i].first < R[j].first)) {
                    out.push_back(S[i]);
                    i++;
                } else if (i == S.size() || R[j].first < S[i].first) {
                    u64 v = (PR - mulm(f, R[j].second)) % PR;
                    if (v) {
                        out.push_back({R[j].first, (uint32_t)v});
                        colrows[R[j].first].push_back(r2);
                        colcnt[R[j].first]++;
                    }
                    j++;
                } else {
                    u64 v = (S[i].second + PR - mulm(f, R[j].second)) % PR;
                    if (v)
                        out.push_back({S[i].first, (uint32_t)v});
                    else
                        colcnt[S[i].first]--;
                    i++;
                    j++;
                }
            }
            S.swap(out);
        }
        for (auto &x : R)
            colcnt[x.first]--;
        rowdead[r] = 1;
        coldead[c] = 1;
        colcnt[c] = 0;
        if (keep)
            piv.push_back({c, std::move(rows[r])});
        vector<pair<int, uint32_t>>().swap(rows[r]);
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
                        if (!rowdead[r] && get(r, c))
                            M.push_back(r);
                    sort(M.begin(), M.end());
                    M.erase(unique(M.begin(), M.end()), M.end());
                    colrows[c].swap(M);
                }
                int best = -1;
                size_t bl = SIZE_MAX;
                for (int r : colrows[c])
                    if (!rowdead[r] && rows[r].size() < bl && get(r, c)) {
                        bl = rows[r].size();
                        best = r;
                    }
                if (best < 0)
                    continue;
                vector<int> touched;
                for (auto &x : rows[best])
                    touched.push_back(x.first);
                pivot(best, c);
                any = true;
                for (int c2 : touched)
                    if (!coldead[c2] && colcnt[c2] > 0)
                        pq.push({colcnt[c2], c2});
            }
        }
    }
    vector<u64> kernel_vector(mt19937_64 &rng) {
        vector<u64> x(nc, 0);
        for (int c = 0; c < nc; c++)
            if (!coldead[c])
                x[c] = rng() % PR;
        for (int i = (int)piv.size() - 1; i >= 0; i--) {
            int c = piv[i].first;
            u64 s = 0, ac = 0;
            for (auto &y : piv[i].second) {
                if (y.first == c)
                    ac = y.second;
                else
                    s = (s + mulm(y.second, x[y.first])) % PR;
            }
            x[c] = mulm((PR - s) % PR, inv(ac));
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
    int deg = atoi(argv[1]);
    PR = strtoull(argv[2], 0, 10);
    unsigned long long seed = argc > 3 ? strtoull(argv[3], 0, 10) : 1;
    string outfn = argc > 4 ? argv[4] : "cert.txt";
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
    vector<char> drop;
    for (int k = 1; k <= deg; k++) {
        FpElim M(C[k].size(), C[k - 1].size(), false);
        for (size_t i = 0; i < C[k].size(); i++) {
            boundary(decode(C[k][i]), bd);
            for (auto &x : bd) {
                int j = idx(k - 1, x.first);
                if (drop.empty() || !drop[j])
                    M.rows[i].push_back({j, (uint32_t)((x.second % (long long)PR + PR) % PR)});
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
        fprintf(stderr, "d_%d over F_%llu: rank %lld\n", k, PR, nd);
    }
    int K1 = deg + 1;
    vector<int> keepidx(C[deg].size(), -1);
    vector<int> kept;
    for (size_t j = 0; j < C[deg].size(); j++)
        if (!drop[j]) {
            keepidx[j] = kept.size();
            kept.push_back(j);
        }
    vector<vector<pair<int, int>>> full(C[K1].size());
    FpElim Mf(C[K1].size(), kept.size(), true);
    vector<vector<pair<int, uint32_t>>> colsT(kept.size());
    for (size_t i = 0; i < C[K1].size(); i++) {
        boundary(decode(C[K1][i]), bd);
        for (auto &x : bd) {
            int j = idx(deg, x.first);
            full[i].push_back({j, x.second});
            uint32_t v = (uint32_t)((x.second % (long long)PR + PR) % PR);
            if (keepidx[j] >= 0) {
                Mf.rows[i].push_back({keepidx[j], v});
                colsT[keepidx[j]].push_back({(int)i, v});
            }
        }
    }
    Mf.finalize();
    FpElim Mw(kept.size(), C[K1].size(), true);
    for (size_t j = 0; j < kept.size(); j++)
        Mw.rows[j] = colsT[j];
    vector<vector<pair<int, uint32_t>>>().swap(colsT);
    Mw.finalize();
    fprintf(stderr, "eliminating cocycle system\n");
    Mf.run();
    fprintf(stderr, "eliminating cycle system\n");
    Mw.run();
    long long kerw = 0;
    for (int c = 0; c < Mw.nc; c++)
        if (!Mw.coldead[c])
            kerw++;
    fprintf(stderr, "dim ker d_%d mod %llu = %lld\n", K1, PR, kerw);
    for (int tw = 0; tw < 40; tw++) {
        vector<u64> w = Mw.kernel_vector(rng);
        unordered_map<int, long long> dw;
        for (size_t i = 0; i < w.size(); i++)
            if (w[i])
                for (auto &x : full[i])
                    dw[x.first] += (long long)w[i] * x.second;
        long long nz = 0;
        for (auto &kv : dw) {
            if (((kv.second % (long long)PR) + PR) % PR) {
                fprintf(stderr, "BUG: boundary not divisible by p\n");
                return 3;
            }
            if (kv.second)
                nz++;
        }
        if (nz == 0)
            continue;
        for (int tf = 0; tf < 20; tf++) {
            vector<u64> fk = Mf.kernel_vector(rng);
            long long pairing = 0;
            for (auto &kv : dw)
                if (keepidx[kv.first] >= 0 && fk[keepidx[kv.first]])
                    pairing = (pairing + (kv.second / (long long)PR) % (long long)PR *
                                             (long long)fk[keepidx[kv.first]]) %
                              (long long)PR;
            pairing = (pairing % (long long)PR + PR) % PR;
            if (pairing != 0) {
                FILE *fo = fopen(outfn.c_str(), "w");
                fprintf(fo, "GRAPH %d %d %d\n", V, E, N);
                for (int e = 0; e < E; e++)
                    fprintf(fo, "EDGE %d %d\n", etail[e], ehead[e]);
                fprintf(fo, "DEGREE %d\nPRIME %llu\n", deg, PR);
                long long nw = 0, nf = 0;
                for (size_t i = 0; i < w.size(); i++)
                    if (w[i]) {
                        fprintf(fo, "W %llu %s\n", w[i], describe(decode(C[K1][i])).c_str());
                        nw++;
                    }
                for (size_t j = 0; j < kept.size(); j++)
                    if (fk[j]) {
                        fprintf(fo, "F %llu %s\n", fk[j], describe(decode(C[deg][kept[j]])).c_str());
                        nf++;
                    }
                fclose(fo);
                printf("{\"certificate\":\"%s\",\"prime\":%llu,\"w_cells\":%lld,\"f_cells\":%lld,\"z_"
                       "support\":%lld,\"tries_w\":%d,\"tries_f\":%d}\n",
                       outfn.c_str(), PR, nw, nf, nz, tw + 1, tf + 1);
                return 0;
            }
        }
    }
    printf("{\"certificate\":null}\n");
    return 1;
}
