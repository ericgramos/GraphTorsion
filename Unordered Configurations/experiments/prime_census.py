#!/usr/bin/env python3
"""
Which cyclic groups occur, and in which graphs?  The evidence for the prime dichotomy.

THE QUESTION
------------
The flip-chain class puts Z/n in H_{n-1}(B_n K_{n+1,2n-1}), and the whole torsion of that
degree is the cokernel of an inclusion matrix.  Is that family the source of ALL the torsion we
have seen, or only some of it?

The answer splits by prime, not by degree.

    q >= 3   every single Z/q summand in the entire data set occurs in a graph that contains
             K_{4,5} as a subgraph.  The Z/3's live in K_{4,b} and K_{5,b}; the Z/4's in
             K_{5,7} and K_{5,8}.  Nothing else, anywhere, has a summand other than Z/2.

    q = 2    not bipartite at all.  Ko--Park's Z/2 is the q = 2 case of the same inclusion
             matrix (coker W_{1,2}(3) = Z/2 for K_{3,3}), but the planar mechanism is
             independent of it: the octahedron, cube, dodecahedron and icosahedron are planar,
             so they have no K_{3,3} minor whatsoever, and they all carry Z/2.  K_5 is another
             witness in degree one -- it has only five vertices, so no K_{3,3} minor, and it
             carries the Ko--Park class.

So Z/2 is generic and everything else is structural.  The conjecture this supports: for
q >= 3, a Z/q summand forces a K_{q+1,2q-1} minor.

WHY K_{q+1,2q-1} IS THE RIGHT GRAPH FOR Z/q
-------------------------------------------
By Wilson's diagonal form, a prime power q is the order of a primary summand of
H_{n-1}(B_n K_{n+1,b}) exactly when q <= n and b >= n+q-1.  Minimising over (n, b) gives n = q, b = 2q-1: the graph K_{q+1,2q-1} at
k = q particles, in degree q-1.  That is the flip-chain family, so it is a carrier; the conjecture
is that for q >= 3 it is the only one, up to minors.

WHAT THIS SCRIPT DOES
---------------------
    (1) Computes tors H_3(B_4 K_{5,7}) and tors H_{3}(B_4 K_{5,8}) directly and compares them
        with coker W_{3,4}(b).  These are the Z/4 cases -- the first non-squarefree torsion in
        a graph braid group -- and they are reachable only through tors H_d = tors coker d_{d+1} plus the sweeping
        peel, since the full complexes have 2.6 and 4.6 million cells.
    (2) Takes a census of the torsion summands stored in data/*.json, by primary type and by graph
        family, and checks that no non-Z/2 summand occurs outside the complete bipartite
        families.
    (3) Lists the Z/2 witnesses with no K_{3,3} minor, which are what make q = 2 exceptional.

    python3 experiments/prime_census.py [--skip-compute]

The two direct computations take under a minute together; the census is instant.  Exits
non-zero if the census finds a non-Z/2 summand in a graph with no K_{4,5} subgraph, so it is a
regression test on the conjecture.
"""

import argparse
import glob
import itertools
import json
import os
import sys
import time
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import complete_bip, lean, paths, primary, reduced_cell_count   # noqa: E402
from gbt.linalg import sparse_smith                                      # noqa: E402

DATA_FILE = "prime_census.json"

#: the (a, b, k, degree) cases where coker W predicts a non-squarefree summand
Z4_CASES = [(5, 7, 4, 3), (5, 8, 4, 3)]
#: smaller bipartite cases, for the Z/3 end of the same family
Z3_CASES = [(4, 5, 3, 2), (4, 9, 3, 2), (4, 11, 3, 2)]


def primary_types(invariant_factors):
    """The multiset of prime powers in a primary decomposition."""
    out = Counter()
    for t in invariant_factors:
        x, p = t, 2
        while x > 1:
            while p * p <= x and x % p:
                p += 1
            if p * p > x:
                p = x
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            out[p ** e] += 1
            p = 2
    return out


def wilson_cokernel_torsion(n, b):
    """Non-unit invariant factors of W_{n-1,n}(b): the predicted torsion of H_{n-1}."""
    rows = list(itertools.combinations(range(b), n - 1))
    ri = {s: i for i, s in enumerate(rows)}
    cols = [{ri[S]: 1 for S in itertools.combinations(T, n - 1)}
            for T in itertools.combinations(range(b), n)]
    _, divs = sparse_smith(cols, verbose=False)
    return [x for x in divs if x not in (0, 1)]


def direct_computations(res):
    """Compute the Z/4 and Z/3 cases and check them against the inclusion matrix."""
    print("== (1) the bipartite family, computed directly and checked against coker W ==\n")
    hdr = ("%-12s %3s %4s %13s %-30s %-30s %s"
           % ("graph", "k", "deg", "cells", "computed", "coker W predicts", "match"))
    print(hdr)
    print("-" * len(hdr))
    ok = True
    for (a, b, k, d) in Z4_CASES + Z3_CASES:
        G = complete_bip(a, b)
        key = "K%d,%d|%d|%d" % (a, b, k, d)
        if key in res:
            t = res[key]["torsion"]
        else:
            t = lean.torsion_in_degree(G, k, d)
            res[key] = {"torsion": t, "a": a, "b": b, "k": k, "degree": d}
            paths.save_json(DATA_FILE, res)
        pred = wilson_cokernel_torsion(a - 1, b)
        same = (primary(t) == primary(pred))
        ok = ok and same
        print("%-12s %3d %4d %13s %-30s %-30s %s"
              % ("K_{%d,%d}" % (a, b), k, d, format(reduced_cell_count(G, k), ","),
                 primary(t) if t else "none", primary(pred), "yes" if same else "NO"),
              flush=True)
    print("-" * len(hdr))
    print("\nThe two Z/4 rows are the first non-squarefree torsion in a graph braid group, and")
    print("they are reachable only through tors H_d = tors coker d_{d+1}: the full complexes have 2.6 and 4.6")
    print("million cells, but only two degree slices are needed.")
    return ok


def census():
    """Every torsion summand in data/*.json, by primary type and family."""
    print("\n== (2) census of the torsion summands stored in data/ ==\n")

    def family_of(fname, key):
        spec = key.split("|")[0]
        if fname.startswith("sweep_bipartite") or fname.startswith("prime_census"):
            return "complete bipartite"
        if fname.startswith("sweep_complete"):
            return "complete K_n"
        if fname.startswith("large_cases"):
            return "large cases (bipartite and complete)"
        if fname.startswith("petersen"):
            return "Petersen family"
        if fname.startswith("linkless"):
            return "linkless controls"
        if fname.startswith("density"):
            return "density interpolation (K_{4,5} to K_9)"
        if fname.startswith("small_graphs"):
            return "exhaustive small graphs"
        if (fname.startswith("planar") or fname.startswith("platonic")
                or fname.startswith("higher_degree") or fname.startswith("genuine")):
            return "planar families"
        del spec
        return fname

    tally = {}
    nonstandard = []
    for f in sorted(glob.glob(os.path.join(paths.DATA, "*.json"))):
        fname = os.path.basename(f)
        with open(f) as fh:
            d = json.load(fh)
        items = d.items() if isinstance(d, dict) else [(str(i), r) for i, r in enumerate(d)]
        for key, rec in items:
            if not isinstance(rec, dict):
                continue
            blocks = []
            if rec.get("H"):
                blocks = [(int(deg), v[1]) for deg, v in rec["H"].items()]
            elif isinstance(rec.get("torsion"), list):
                blocks = [(rec.get("degree", -1), rec["torsion"])]
            elif isinstance(rec.get("H2_torsion"), list):
                blocks = [(2, rec["H2_torsion"])]
            fam = family_of(fname, key)
            for deg, t in blocks:
                if not t:
                    continue
                types = primary_types(t)
                tally.setdefault(fam, Counter()).update(types)
                for q in types:
                    if q != 2:
                        nonstandard.append((fam, fname, key, deg, q, primary(t)))

    print("%-38s %s" % ("family", "primary summands found (Z/q with multiplicity)"))
    print("-" * 96)
    for fam in sorted(tally):
        c = tally[fam]
        print("%-38s %s" % (fam, "   ".join("Z/%d x%d" % (q, m)
                                            for q, m in sorted(c.items()))))
    print("-" * 96)
    total = sum(sum(c.values()) for c in tally.values())
    odd = sum(m for c in tally.values() for q, m in c.items() if q != 2)
    print("total: %d primary summands counted, %d of them not Z/2" % (total, odd))
    print("(one count per stored record: a graph stored in two files is counted twice; the")
    print(" small-graph sweep files store only a summary string and are not included)")

    print("\nevery non-Z/2 summand, with the graph it came from "
          "(checked: the graph contains K_{4,5} as a subgraph):")
    bad = []
    for fam, fname, key, deg, q, shown in sorted(set(nonstandard)):
        G = graph_of(fname, key)
        ok = G is not None and contains_K45(G)
        print("   %-34s H_%-2d Z/%d   %-26s %s"
              % (key, deg, q, shown,
                 "contains K_{4,5}" if ok else "  <== NO K_{4,5} SUBGRAPH FOUND"))
        if not ok:
            bad.append((key, deg, q))
    if not nonstandard:
        print("   (none)")
    return bad


def graph_of(fname, key):
    """Rebuild the graph behind a data record, for the records that can carry odd torsion.

    Returns None if the key is not recognised, which the census then reports as a failure
    rather than silently passing.
    """
    import re
    spec = key.split("|")[0]
    m = re.fullmatch(r"K(\d+),(\d+)", spec)
    if m:
        return complete_bip(int(m.group(1)), int(m.group(2)))
    legacy = {"K55": (5, 5), "K56": (5, 6), "K57": (5, 7)}
    if spec in legacy:
        return complete_bip(*legacy[spec])
    if fname.startswith("density"):
        # the six graphs of density_interpolation.py, rebuilt here so that this script does
        # not import another one: K_{4,5} on A = {0..3}, B = {4..8}, plus edges inside a side
        A, B = range(4), range(4, 9)
        base = [(a, b) for a in A for b in B]
        AA = list(itertools.combinations(A, 2))
        BB = list(itertools.combinations(B, 2))
        cases = {"K_{4,5}": base, "K_{4,5} + 3 in B": base + BB[:3],
                 "K_{4,5} + 6 in B": base + BB[:6], "K_{4,5} + all A": base + AA,
                 "K_{4,5} + all B": base + BB, "K_9": base + AA + BB}
        if spec in cases:
            from gbt import Graph
            return Graph(9, cases[spec])
    return None


def contains_K45(G):
    """True if G has a subgraph isomorphic to K_{4,5} (networkx subgraph monomorphism)."""
    import networkx as nx
    from networkx.algorithms import isomorphism
    big = nx.Graph(G.edges)
    return isomorphism.GraphMatcher(big, nx.complete_bipartite_graph(4, 5)) \
        .subgraph_is_monomorphic()


def two_is_different():
    """The Z/2 witnesses that have no K_{3,3} minor at all."""
    print("\n== (3) why q = 2 is exceptional ==\n")
    import networkx as nx
    from gbt import cube, dodecahedron, icosahedron, octahedron, complete
    rows = [("K_5", complete(5), 2, 1, "Ko--Park class; only 5 vertices, so no K_{3,3} minor"),
            ("octahedron", octahedron(), 4, 2, "planar"),
            ("cube Q_3", cube(), 4, 2, "planar"),
            ("dodecahedron", dodecahedron(), 4, 2, "planar"),
            ("icosahedron", icosahedron(), 4, 2, "planar")]
    print("%-16s %3s %4s %-9s %-9s %s"
          % ("graph", "k", "deg", "planar", "K_{3,3}?", "note"))
    for name, G, k, d, note in rows:
        Gx = nx.Graph(G.edges)
        pl = nx.check_planarity(Gx)[0]
        # a planar graph has no K_{3,3} minor; and a graph on < 6 vertices cannot have one
        no_k33 = pl or G.N < 6
        print("%-16s %3d %4d %-9s %-9s %s"
              % (name, k, d, pl, "no" if no_k33 else "maybe", note))
    print("\nAll five carry Z/2 and none can have a K_{3,3} minor, so the q = 2 case of")
    print("the odd-torsion conjecture is false -- which is exactly the point of the dichotomy.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-compute", action="store_true",
                    help="census only, no new homology computations")
    args = ap.parse_args()

    res = paths.load_json(DATA_FILE)
    ok = True
    if not args.skip_compute:
        t0 = time.time()
        ok = direct_computations(res)
        print("   [%.0fs]" % (time.time() - t0))
    bad = census()
    two_is_different()

    print()
    if not ok:
        print("A direct computation DISAGREED with coker W.")
        return 1
    if bad:
        print("CONJECTURE VIOLATED: non-Z/2 torsion with no bipartite structure: %s" % bad)
        return 1
    print("Every non-Z/2 summand in the stored data lives in a graph containing K_{4,5}.")
    print("Every other stored torsion summand is Z/2.")
    if not args.skip_compute:
        print("\nwrote data/%s" % DATA_FILE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
