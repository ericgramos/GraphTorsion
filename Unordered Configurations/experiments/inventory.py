#!/usr/bin/env python3
"""
A single accounting of every stored computation in data/*.json.

WHAT IT DOES
------------
Every sweep and every individual experiment writes its results to its own JSON file, in its own
format (the large search torsion_hunt.py writes JSON lines, data/torsion_hunt.jsonl; its
timed-out jobs, if any, are left out).  This script reads all of them, rebuilds the graph behind every record, identifies graphs
up to isomorphism (Weisfeiler--Lehman hash, then a full isomorphism test), and writes one table:

    data/inventory.csv     one row per stored computation, with columns
        source     the script (or sweep file) that produced the record
        graph      the label used by that script
        class      an integer naming the isomorphism class of the graph (shared across sources)
        N, M       number of vertices and edges
        k          number of particles
        computed   what was computed: "all degrees", "tors H_d" (the torsion of one degree only),
                   "H_n (top) rank only", or the order of one explicit class
        result     the torsion found, by degree, as a primary decomposition, or "torsion-free"

It then prints a summary: records, distinct graphs and distinct (graph, k) pairs per source, how
many carry torsion in some H_i with i >= 2, how many carry a summand other than Z/2, and the same
restricted to planar graphs.  The numbers printed when this file was last run were

    10,653 records, 8,140 distinct graphs, 10,062 distinct (graph, k) pairs, k from 2 to 12;
    3,991 graphs with torsion in some H_i, i >= 2;   207 graphs with a summand other than Z/2;
    5,714 planar graphs, 2,032 of them with torsion in some H_i, i >= 2, and none with a summand
    other than Z/2.

One record (the theta-graph generator of two_torsion_anatomy.py) has no single graph attached
and is given class -1.

THE ONE EXCEPTION TO "NO SCRIPT IMPORTS ANOTHER"
------------------------------------------------
Several experiments define their graphs inside the script (the planar families of
planar_odd_search.py, the controls of linkless_controls.py, ...).  Rather than duplicate those
definitions, this script loads those modules with importlib and calls their `cases()` /
`families()` / `controls()` functions.  Loading a module runs only its definitions: every
experiment guards its work with `if __name__ == "__main__"`.

WHAT IS PRINTED
---------------
A table with one row per source and columns

    source       the script or small-graph sweep file ("sweep_small_graphs:n7_k4_d1" etc.)
    records      number of stored records from that source
    graphs       number of distinct isomorphism classes among them
    k            the range of particle numbers
    tors H>=2    records showing torsion in some H_i with i >= 2 ("-" for sources that record
                 no torsion at all, such as the top-degree ranks)
    non-Z/2      records showing a summand Z/m with m >= 3

followed by the totals quoted above and "wrote data/inventory.csv".  Excerpt:

    source                               records graphs k         tors H>=2  non-Z/2
    sweep_complete                            46     10 2-12             10        0
    sweep_bipartite                           92     27 2-7              25        7
    large_cases                                4      4 3-4               4        1

Per-source columns count RECORDS, so a graph stored by two scripts is counted in both rows;
the totals at the bottom count distinct graphs.

USAGE
-----
    python3 experiments/inventory.py          (a few seconds, well under 1 GB; needs networkx)

No options.  It only reads data/*.json and overwrites data/inventory.csv; rerun it after any
new computation.  Exit code: 0, or a Python traceback if a stored record names a graph the
registry does not know (add the graph to `registry()` when a new experiment stores results).

PITFALLS
--------
* "Torsion in some H_i, i >= 2" and "summand other than Z/2" are read from the stored result
  strings by regular expressions (`torsion_above_one`, `not_two`).  If a new experiment
  stores results in a new string format, check that these still parse it.
* Graph identity is decided by isomorphism of the underlying simple graph; the labels that the
  different scripts use are ignored.
"""

import csv
import glob
import importlib.util
import json
import os
import re
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import warnings                                                                   # noqa: E402

import networkx as nx                                                             # noqa: E402

warnings.filterwarnings("ignore", message="The hashes produced")

from gbt import (Graph, antiprism, complete, complete_bip, from_networkx,         # noqa: E402
                 octahedron, petersen_family, platonic_solids, prism, prism_graph)
from gbt.chain import primary                                                     # noqa: E402


def load_experiment(name):
    """Import experiments/<name>.py as a module without running its main()."""
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, "experiments",
                                                                     name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def data(name):
    """The parsed contents of data/<name>.json."""
    with open(os.path.join(ROOT, "data", name + ".json")) as fh:
        return json.load(fh)


def to_nx(G):
    """The networkx graph on vertices 0..N-1 with the edges of the gbt Graph G."""
    g = nx.Graph()
    g.add_nodes_from(range(G.N))
    g.add_edges_from(G.edges)
    return g


def tors_str(t):
    """A list of invariant factors as a primary decomposition, or "none" if it is empty."""
    return primary(t) if t else "none"


def torsion_by_degree(H):
    """'H2:(Z/2)^4 + Z/3; H3:Z/2' for a stored homology dict, or 'torsion-free'."""
    parts = ["H%s:%s" % (d, tors_str(v[1]))
             for d, v in sorted(H.items(), key=lambda x: int(x[0])) if v[1]]
    return "; ".join(parts) if parts else "torsion-free"


def registry():
    """Every named graph that a stored record refers to, keyed as the records key them."""
    reg = {}
    for n in range(3, 14):
        reg["K%d" % n] = complete(n)
    for a in range(2, 8):
        for b in range(a, 12):
            reg["K%d,%d" % (a, b)] = complete_bip(a, b)
    for name, G in petersen_family():
        reg[("petersen_family", name)] = G
    for name, G, _ in load_experiment("planar_odd_search").families():
        reg[("planar_odd_search", name)] = G
    for name, G, _ in load_experiment("planar_torsion").cases():
        reg[("planar_torsion", name)] = G
    for name, G, _, _ in load_experiment("linkless_controls").controls():
        reg[("linkless_controls", name)] = G
    for name, E, _ in load_experiment("density_interpolation").cases():
        reg[("density_interpolation", name)] = from_networkx(nx.Graph(E))
    hd = load_experiment("higher_degree_torsion")
    reg.update({("higher_degree_torsion", "octahedron"): octahedron(),
                ("higher_degree_torsion", "octahedron + 1 x C_3"): hd.octahedron_plus_triangles(1),
                ("higher_degree_torsion", "octahedron + 2 x C_3"): hd.octahedron_plus_triangles(2),
                ("higher_degree_torsion", "oct-C_3 joined"): hd.octahedron_triangle_joined()})
    solids = dict(platonic_solids())
    for n in range(3, 7):
        solids["prism_%d" % n] = prism_graph(n)
        solids["antiprism_%d" % n] = antiprism(n)
    for name, G in solids.items():
        reg[("platonic_solids", name)] = G
    reg.update({("two_torsion_anatomy", "K_4"): complete(4),
                ("two_torsion_anatomy", "K_5"): complete(5),
                ("two_torsion_anatomy", "K_{3,3}"): complete_bip(3, 3),
                ("two_torsion_anatomy", "octahedron"): octahedron(),
                ("two_torsion_anatomy", "prism"): prism()})
    return reg


def collect():
    """[(source, label, Graph or None, k, what was computed, result)] over all of data/."""
    reg = registry()
    rows = []

    def add(src, label, G, k, scope, res):
        """Append one record (source, label, Graph, k, what was computed, result string)."""
        rows.append(dict(source=src, label=label, G=G, k=int(k), scope=scope, result=res))

    for src in ("sweep_complete", "sweep_bipartite"):
        for key, v in data(src).items():
            g, k = key.split("|")
            add(src, g, reg[g], k, "all degrees", torsion_by_degree(v["H"]))
    large = {"K12-3": ("K12", 3), "K8-4": ("K8", 4), "K55": ("K5,5", 4), "K56": ("K5,6", 4)}
    for key, v in data("large_cases").items():
        g, k = large[key]
        add("large_cases", g, reg[g], k, "all degrees", torsion_by_degree(v["H"]))
    for key, v in data("top_degree_vanishing").items():
        g, n = key.split("|")
        add("top_degree_vanishing", g, reg[g], n, "H_%s (top) rank only" % n,
            "rank %d" % v["H_n"][0])
    for key, v in data("prime_census").items():
        g, k, d = key.split("|")
        add("prime_census", g, reg[g], k, "tors H_%s" % d, tors_str(v["torsion"]))
    for src in ("petersen_family", "linkless_controls", "planar_torsion", "planar_odd_search",
                "density_interpolation"):
        for key, v in data(src).items():
            g, k = key.split("|")
            add(src, g, reg[(src, g)], k, "all degrees", torsion_by_degree(v["H"]))
    for key, v in data("platonic_solids").items():
        g, k = key.split("|")
        G = reg[("platonic_solids", g)]
        if "H" in v:
            add("platonic_solids", g, G, k, "all degrees", torsion_by_degree(v["H"]))
        else:
            add("platonic_solids", g, G, k, "tors H_2", tors_str(v["H2_torsion"]))
    for key, v in data("higher_degree_torsion").items():
        g, k, d = key.split("|")
        add("higher_degree_torsion", g, reg[("higher_degree_torsion", g)], k, "tors H_%s" % d,
            tors_str(v["torsion"]))
    for key, v in data("genuine_higher_torsion").items():
        if key.startswith("deep"):
            _, g, k, d = key.split("|")
            G = {"octahedron": octahedron(), "antiprism_4": antiprism(4)}[g]
            add("genuine_higher_torsion", g, G, k, "tors H_%s" % d, tors_str(v["torsion"]))
        else:
            G = Graph(v["N"], [tuple(e) for e in v["edgelist"]])
            add("genuine_higher_torsion", "planar N=%d M=%d" % (v["N"], v["M"]), G, v["k"],
                "tors H_%d" % v["degree"], tors_str(v["torsion"]))
    for path in sorted(glob.glob(os.path.join(ROOT, "data", "small_graphs_*.json"))):
        name = os.path.basename(path)[:-5]
        for v in data(name):
            G = Graph(len(v["degrees"]), [tuple(e) for e in v["edgelist"]])
            add(name.replace("small_graphs_", "sweep_small_graphs:"), "N=%d M=%d" % (G.N, G.M),
                G, v["k"], "all degrees",
                v["torsion"] if v["torsion"] != "none" else "torsion-free")
    # random supergraph walks from known carriers (supergraph_walks.py): torsion of H_2 only
    if os.path.exists(os.path.join(ROOT, "data", "supergraph_walks.json")):
        for key, v in data("supergraph_walks").items():
            G = Graph(v["N"], [tuple(e) for e in v["edges"]])
            add("supergraph_walks", key, G, v["k"], "tors H_2", tors_str(v["torsion"]["2"]))
    # the large search (torsion_hunt.py): torsion of H_2 only; one JSON object per line
    hunt = os.path.join(ROOT, "data", "torsion_hunt.jsonl")
    if os.path.exists(hunt):
        with open(hunt) as fh:
            for line in fh:
                v = json.loads(line)
                if v.get("summary") is None:          # a job that timed out
                    continue
                g = nx.convert_node_labels_to_integers(nx.Graph([tuple(e) for e in v["edges"]]))
                G = Graph(g.number_of_nodes(), sorted((min(a, b), max(a, b)) for a, b in g.edges()))
                add("torsion_hunt:" + v["family"], v["name"], G, v["k"], "tors H_2",
                    "none" if v["summary"] == "0" else v["summary"])
    for key, v in data("two_torsion_anatomy").items():
        if "|" in key:
            part, g = key.split("|")
            add("two_torsion_anatomy", "%s:%s" % (part, g), reg[("two_torsion_anatomy", g)],
                2 if part == "A" else 4, "order of an explicit class", v.get("order", str(v)))
        else:
            add("two_torsion_anatomy", "B:theta(3,3)", None, 4, "H_2 generator", str(v)[:40])
    return rows


def classify(rows):
    """Attach an isomorphism-class number and a planarity flag to every row."""
    by_hash = {}                       # WL hash -> [(class number, networkx graph)]
    count = 0
    for r in rows:
        if r["G"] is None:
            r["cls"], r["planar"] = -1, None
            continue
        g = to_nx(r["G"])
        h = nx.weisfeiler_lehman_graph_hash(g)
        for i, gg in by_hash.get(h, []):     # only graphs with the same hash can be isomorphic
            if nx.is_isomorphic(g, gg):
                r["cls"] = i
                break
        else:
            by_hash.setdefault(h, []).append((count, g))
            r["cls"] = count
            count += 1
        r["planar"] = nx.check_planarity(g)[0]


def torsion_above_one(r):
    """True/False: does the record show torsion in some H_i with i >= 2?  None: not applicable."""
    s, scope = r["result"], r["scope"]
    if scope.startswith("tors H_"):
        return s not in ("none", "torsion-free") and int(scope.split("_")[1]) >= 2
    if "rank" in s or scope.startswith("order") or scope.startswith("H_2 gen"):
        return None
    return bool(re.search(r"H_?([2-9])", s))


def not_two(r):
    """Does the record show a torsion summand other than Z/2?"""
    return bool(re.search(r"Z/([3-9]|\d\d)", r["result"]))


def main():
    """Collect and classify every record, write data/inventory.csv, print the summary."""
    rows = collect()
    classify(rows)

    out = os.path.join(ROOT, "data", "inventory.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["source", "graph", "class", "N", "M", "k", "computed", "result"])
        for r in rows:
            G = r["G"]
            w.writerow([r["source"], r["label"], r["cls"], G.N if G else "", G.M if G else "",
                        r["k"], r["scope"], r["result"]])

    by_src = defaultdict(list)
    for r in rows:
        by_src[r["source"]].append(r)
    print("%-36s %7s %6s %-8s %10s %8s" % ("source", "records", "graphs", "k", "tors H>=2",
                                           "non-Z/2"))
    for src, rs in by_src.items():
        ks = sorted({r["k"] for r in rs})
        flags = [torsion_above_one(r) for r in rs]
        above = sum(1 for f in flags if f) if any(f is not None for f in flags) else "-"
        print("%-36s %7d %6d %-8s %10s %8d"
              % (src, len(rs), len({r["cls"] for r in rs}),
                 "%d-%d" % (ks[0], ks[-1]) if len(ks) > 1 else str(ks[0]), above,
                 sum(1 for r in rs if not_two(r))))

    real = [r for r in rows if r["cls"] >= 0]
    print("\ntotal: %d records, %d distinct graphs, %d distinct (graph, k) pairs, k from %d to %d"
          % (len(rows), len({r["cls"] for r in real}), len({(r["cls"], r["k"]) for r in real}),
             min(r["k"] for r in rows), max(r["k"] for r in rows)))
    print("graphs with torsion in some H_i, i >= 2 : %d"
          % len({r["cls"] for r in real if torsion_above_one(r)}))
    print("graphs with a summand other than Z/2   : %d"
          % len({r["cls"] for r in real if not_two(r)}))
    planar = [r for r in real if r["planar"]]
    print("planar graphs: %d; with torsion in some H_i, i >= 2: %d; with a summand other "
          "than Z/2: %d"
          % (len({r["cls"] for r in planar}),
             len({r["cls"] for r in planar if torsion_above_one(r)}),
             len({r["cls"] for r in planar if not_two(r)})))
    print("\nwrote %s" % os.path.relpath(out, ROOT))


if __name__ == "__main__":
    main()
