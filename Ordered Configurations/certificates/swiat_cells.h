// =====================================================================================================
//  swiat_cells.h  --  cells of the reduced ordered Swiatkowski complex, shared by certq, orderq, exactsolve
//                     and reducecyc
// =====================================================================================================
//  CONVENTIONS (the same as engine/oconf.cpp and certificates/verify_cert.py).
//   * The graph has V vertices and E edges; edge e is oriented tail -> head.  Half-edge 2e is the tail end of
//     edge e and half-edge 2e+1 its head end.  The BASE half-edge h0(v) of a vertex is its lowest-numbered
//     half-edge.  (Reordering the edges therefore changes the bases; the tools that study a subgraph H list
//     the edges of H first so that every vertex of H has its base inside H.)
//   * A cell of the reduced complex of Conf_n: each vertex is empty or carries one particle p together with a
//     non-base half-edge h (this stands for the generator h - h0(v)); each edge carries an ordered word of
//     particles, listed from tail to head; every particle 0..n-1 occurs exactly once.  The degree of a cell
//     is the number of occupied vertices.
//   * Boundary: if v_0 < v_1 < ... are the occupied vertices and v_j carries (p, h), then
//        d = sum_j (-1)^j ( [p pushed into edge(h) at the end of h] - [p pushed into edge(h0) at the end of h0] ),
//     where pushing in at a tail end prepends p to the word and at a head end appends it.
//   * Text form of a cell (used in every certificate and cycle file):
//        V[a/b,-,...] E[w_0|w_1|...]   vertex v carries particle a with half-edge b, or is empty ("-");
//                                      w_e is the word on edge e, particles separated by spaces.
//
//  KEYS.  A cell is packed into a 128-bit integer: particle p owns bits [12p, 12p+12).  Its low 8 bits hold
//  its location (the half-edge h if it sits at a vertex, 2E + e if it lies on edge e) and the high 4 bits its
//  rank in the word of that edge.  So n <= 10, 2E + E < 256 and words have length < 16.  Sorting the keys of
//  all cells of one degree gives the column/row numbering used by the eliminators.
// =====================================================================================================
#pragma once
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef unsigned __int128 Key;

int V, E, N;                    // vertices, edges, particles
vector<int> etail, ehead;       // edge e = etail[e] -> ehead[e]
vector<int> hvert;              // half-edge -> its vertex
vector<vector<int>> vhalf;      // vertex -> its half-edges, in increasing order
vector<int> base;               // vertex -> base half-edge (-1 for an isolated vertex)

struct Cell {
    vector<int> vocc;           // vertex -> particle sitting there, or -1
    vector<int> vh;             // vertex -> the half-edge h of its generator h - h0 (when occupied)
    vector<vector<int>> elist;  // edge -> word of particles, tail to head
};

// Read "V E N" followed by E lines "tail head" from stdin and set up the half-edge tables.
bool read_graph() {
    if (scanf("%d %d %d", &V, &E, &N) != 3) return false;
    etail.resize(E); ehead.resize(E); hvert.resize(2 * E); vhalf.assign(V, {});
    for (int e = 0; e < E; e++) {
        if (scanf("%d %d", &etail[e], &ehead[e]) != 2) return false;
        hvert[2 * e] = etail[e];
        hvert[2 * e + 1] = ehead[e];
        vhalf[etail[e]].push_back(2 * e);
        vhalf[ehead[e]].push_back(2 * e + 1);
    }
    base.assign(V, -1);
    for (int v = 0; v < V; v++)
        if (!vhalf[v].empty()) base[v] = vhalf[v][0];   // lowest half-edge, as vhalf is filled in order
    if (N > 10 || 3 * E >= 256) { fprintf(stderr, "graph too large for the 128-bit keys\n"); return false; }
    return true;
}

// Pack a cell into a 128-bit key: 12 bits per particle, the low 8 its location (half-edge if at a vertex, 2E + e if on edge e) and the high 4 its rank in the word of that edge.
Key encode(const Cell &c) {
    Key s = 0;
    for (int v = 0; v < V; v++)
        if (c.vocc[v] >= 0) s |= (Key)(unsigned)(c.vh[v]) << (12 * c.vocc[v]);
    for (int e = 0; e < E; e++)
        for (int r = 0; r < (int)c.elist[e].size(); r++)
            s |= (Key)(unsigned)((2 * E + e) | (r << 8)) << (12 * c.elist[e][r]);
    return s;
}

// Inverse of encode: rebuild the occupants of the vertices, their half-edges and the edge words from a key.
Cell decode(Key s) {
    Cell c;
    c.vocc.assign(V, -1); c.vh.assign(V, -1); c.elist.assign(E, {});
    vector<vector<pair<int, int>>> tmp(E);           // (rank, particle) per edge
    for (int p = 0; p < N; p++) {
        unsigned w = (unsigned)((s >> (12 * p)) & 0xFFF);
        int loc = w & 0xFF, r = w >> 8;
        if (loc < 2 * E) { int v = hvert[loc]; c.vocc[v] = p; c.vh[v] = loc; }
        else tmp[loc - 2 * E].push_back({r, p});
    }
    for (int e = 0; e < E; e++) {
        sort(tmp[e].begin(), tmp[e].end());
        for (auto &x : tmp[e]) c.elist[e].push_back(x.second);
    }
    return c;
}

// All cells of degree k: particles are placed one at a time, either at an empty vertex with one of its non-base
// half-edges or at any position of any edge word.
void enumerate(int k, vector<Key> &out) {
    Cell c;
    c.vocc.assign(V, -1); c.vh.assign(V, -1); c.elist.assign(E, {});
    int moving = 0;                                   // particles placed at vertices so far
    function<void(int)> rec = [&](int p) {
        if (p == N) { if (moving == k) out.push_back(encode(c)); return; }
        if (moving < k)
            for (int v = 0; v < V; v++)
                if (c.vocc[v] < 0)
                    for (int h : vhalf[v])
                        if (h != base[v]) {
                            c.vocc[v] = p; c.vh[v] = h; moving++;
                            rec(p + 1);
                            moving--; c.vocc[v] = -1; c.vh[v] = -1;
                        }
        if (k - moving <= N - p - 1)                  // enough particles left to reach k occupied vertices
            for (int e = 0; e < E; e++)
                for (int pos = 0; pos <= (int)c.elist[e].size(); pos++) {
                    c.elist[e].insert(c.elist[e].begin() + pos, p);
                    rec(p + 1);
                    c.elist[e].erase(c.elist[e].begin() + pos);
                }
    };
    rec(0);
}

// Sorted cell list of one degree -> index lookup (-1 if absent).
int lookup(const vector<Key> &Ck, Key key) {
    auto it = lower_bound(Ck.begin(), Ck.end(), key);
    return (it == Ck.end() || *it != key) ? -1 : (int)(it - Ck.begin());
}

// Text form of a cell, V[p/h,...] E[w_0|w_1|...], the format of certificate and cycle files.
string describe(const Cell &c) {
    string s = "V[";
    for (int v = 0; v < V; v++) {
        if (v) s += ",";
        s += c.vocc[v] < 0 ? string("-") : to_string(c.vocc[v]) + "/" + to_string(c.vh[v]);
    }
    s += "] E[";
    for (int e = 0; e < E; e++) {
        if (e) s += "|";
        for (size_t r = 0; r < c.elist[e].size(); r++) { if (r) s += " "; s += to_string(c.elist[e][r]); }
    }
    return s + "]";
}

// Inverse of describe().  The E[...] part may list FEWER edges than the graph has (a cell written for a
// subgraph formed by the first edges), or more, provided the extra slots are empty.
bool parse_cell(const string &s, Cell &c) {
    c.vocc.assign(V, -1); c.vh.assign(V, -1); c.elist.assign(E, {});
    size_t a = s.find("V["), b = s.find("] E["), z = s.rfind("]");
    if (a == string::npos || b == string::npos || z == string::npos) return false;
    string vs = s.substr(a + 2, b - a - 2), es = s.substr(b + 4, z - b - 4);
    int v = 0;
    stringstream sv(vs); string tok;
    while (getline(sv, tok, ',')) {
        if (v >= V) return false;
        if (tok != "-") { int p, h; if (sscanf(tok.c_str(), "%d/%d", &p, &h) != 2) return false; c.vocc[v] = p; c.vh[v] = h; }
        v++;
    }
    int e = 0; size_t st = 0;
    while (true) {
        size_t en = es.find('|', st);
        string part = es.substr(st, en == string::npos ? string::npos : en - st);
        stringstream sp(part); int p;
        while (sp >> p) { if (e >= E) return false; c.elist[e].push_back(p); }
        e++;
        if (en == string::npos) break;
        st = en + 1;
    }
    return v == V;
}

// Signed boundary of a cell, as a list of (face key, +-1).
void boundary(const Cell &c, vector<pair<Key, int>> &out) {
    out.clear();
    int idx = 0;                                       // position of the vertex among the occupied ones
    for (int v = 0; v < V; v++)
        if (c.vocc[v] >= 0) {
            int p = c.vocc[v], h = c.vh[v], h0 = base[v];
            int sgn = (idx % 2 == 0) ? 1 : -1;
            idx++;
            for (int t = 0; t < 2; t++) {
                int hh = t == 0 ? h : h0;              // first the generator half-edge, then the base
                Cell f = c;
                f.vocc[v] = -1; f.vh[v] = -1;
                int e = hh / 2;
                if (hh % 2 == 0) f.elist[e].insert(f.elist[e].begin(), p);   // tail end: prepend
                else f.elist[e].push_back(p);                               // head end: append
                out.push_back({encode(f), t == 0 ? sgn : -sgn});
            }
        }
}

// A chain read from a text block: (cell index in the sorted list of its degree, integer coefficient).
// Lines are "coef cell"; lines starting with any other tag are skipped by the callers.
inline bool is_coef_line(const string &line) {
    return !line.empty() && (isdigit((unsigned char)line[0]) || line[0] == '-');
}
