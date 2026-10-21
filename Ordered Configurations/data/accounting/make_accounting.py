"""make_accounting.py -- one table of every homology computation recorded in data/ (and the history logs).

Each row is one computation: a graph, a number of particles n, and a method (exact over Z, ranks over
fields, local Smith form, a quotient complex, the vertex explosion shortcut, ...).  The script reads the
recorded result files, normalises them, decides planarity of each graph with networkx, identifies
isomorphic graphs across sources (so that "distinct graphs" is counted correctly), and writes
    accounting.csv          one row per computation
    accounting_graphs.csv   one row per graph up to isomorphism (with all the names it was run under)
    accounting_summary.txt  counts by category, by n and by planarity
Run from this folder:  python3 make_accounting.py
"""
import ast, csv, itertools, json, os, re, collections
import networkx as nx
from networkx.algorithms.isomorphism import MultiGraphMatcher

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..")
EARLY = os.path.join(D, "nonplanar", "earlier")
HLOG = os.path.join(D, "..", "history", "logs")

K = lambda m: list(itertools.combinations(range(m), 2))
Kb = lambda a, b: [(i, a + j) for i in range(a) for j in range(b)]
T = lambda a, b, c: [(0, 1)] * a + [(0, 2)] * b + [(1, 2)] * c
C = lambda m: [(i, (i + 1) % m) for i in range(m)]
NAMED = {"K4": K(4), "K5": K(5), "K6": K(6), "K33": Kb(3, 3)}

# ---------------------------------------------------------------- graph identification
_reps = []          # list of (nx.MultiGraph, label)
def gid(edges, label):
    """Return a canonical identifier for the isomorphism class of the multigraph (loops allowed)."""
    G = nx.MultiGraph()
    vs = sorted({x for e in edges for x in e})
    G.add_nodes_from(vs)
    G.add_edges_from([tuple(e) for e in edges])
    inv = (G.number_of_nodes(), G.number_of_edges(), tuple(sorted(d for _, d in G.degree())))
    for i, (H, inv2, lab) in enumerate(_reps):
        if inv2 == inv and MultiGraphMatcher(G, H).is_isomorphic():
            return i, lab
    _reps.append((G, inv, label))
    return len(_reps) - 1, label

def planar(edges):
    """Is the underlying simple graph planar? (networkx)."""
    G = nx.Graph()
    G.add_edges_from([tuple(e) for e in edges if e[0] != e[1]])
    return nx.check_planarity(G)[0]

rows = []
def add(cat, label, edges, n, method, degrees, result, status, source):
    """Append one computation to the accounting table, with the isomorphism class of its graph."""
    i, canon = gid(edges, label)
    rows.append(dict(category=cat, graph=label, graph_class=f"G{i:03d} ({canon})", V=len({x for e in edges for x in e}),
                     E=len(edges), planar="yes" if planar(edges) else "no", n=n, method=method, degrees=degrees,
                     result=result, status=status, source=source))

def tors_str(t):
    """t: list per degree of elementary divisors (>1)."""
    parts = []
    for i, ds in enumerate(t):
        if ds:
            c = collections.Counter(ds)
            parts.append(f"H{i}: " + " + ".join(f"(Z/{d})^{m}" if m > 1 else f"Z/{d}" for d, m in sorted(c.items())))
    return "; ".join(parts) or "no torsion"

# ---------------------------------------------------------------- 1. the exact survey
for l in open(os.path.join(D, "survey", "survey_exact.jsonl")):
    d = json.loads(l)
    add("exact survey", d["name"], d["edges"], d["n"], "exact over Z", "all", tors_str(d["torsion"]), "complete",
        "survey/survey_exact.jsonl")

# ---------------------------------------------------------------- 2. triangle family
for l in open(os.path.join(D, "survey", "family.jsonl")):
    d = json.loads(l); a, b, c = d["abc"]
    add("triangle family", f"T({a},{b},{c})", T(a, b, c), d["n"], "exact over Z", "all", tors_str(d["torsion"]),
        "complete", "survey/family.jsonl")
for l in open(os.path.join(D, "survey", "triangle_local_snf.jsonl")):
    d = json.loads(l)
    tor = "; ".join(f"p={p}: " + tors_str(t) for p, t in d["torsion"].items())
    add("triangle family", "T(2,2,3)", d["edges"], d["n"], "local Smith form mod 2^4 and 3^3 + rational ranks", "all", tor,
        "complete", "survey/triangle_local_snf.jsonl")
for l in open(os.path.join(D, "survey", "modp.jsonl")):
    d = json.loads(l)
    add("triangle family", "T(2,2,3)", d["edges"], d["n"], "ranks over F_2, F_3, F_1000003 (universal coefficients)", "all",
        "torsion summands " + json.dumps(d["torsion_counts"]), "complete", "survey/modp.jsonl")
for l in open(os.path.join(D, "survey", "triangle_family_shortcut.jsonl")):
    d = json.loads(l); a, b, c = d["abc"]
    res = "; ".join(f"p={p}: {k} summands" for p, k in d["torsion"].items() if k) or "no torsion"
    add("triangle family", f"T({a},{b},{c})", T(a, b, c), d["n"], "vertex explosion shortcut (prod2 + g2tors), p = 2,3,5,7",
        "H2", "H2: " + res, "complete", "survey/triangle_family_shortcut.jsonl")

# ---------------------------------------------------------------- 3. quotients by Young subgroups, unordered
for l in open(os.path.join(D, "survey", "colored.jsonl")):
    d = json.loads(l)
    add("Young quotients", d["name"], NAMED[d["name"]], len(d["colors"]), f"exact over Z, Young quotient {d['colors']}",
        "all", tors_str(d["torsion"]), "complete", "survey/colored.jsonl")
for l in open(os.path.join(D, "survey", "colored_T.jsonl")):
    d = json.loads(l); a, b, c = d["abc"]
    add("Young quotients", f"T({a},{b},{c})", T(a, b, c), len(d["colors"]), f"exact over Z, Young quotient {d['colors']}",
        "all", tors_str(d["torsion"]), "complete", "survey/colored_T.jsonl")
for l in open(os.path.join(D, "survey", "unordered_scan.jsonl")):
    d = json.loads(l)
    add("unordered", d["name"], NAMED[d["name"]], d["n"], "exact over Z, unordered complex (python/unordered.py)", "all",
        tors_str(d["torsion"]), "complete", "survey/unordered_scan.jsonl")

# ---------------------------------------------------------------- 4. fat graphs and their group quotients
FAT = {"C4 tripled": C(4) * 3, "K23 doubled": Kb(2, 3) * 2, "C4 doubled": C(4) * 2,
       "C4 doubled + diagonal": C(4) * 2 + [(0, 2)]}
for l in open(os.path.join(EARLY, "fat.jsonl")):
    d = json.loads(l)
    sig = {p: v for p, v in d["torsion_signal"].items() if v}
    add("fat graphs", d["name"], FAT[d["name"]], d["n"],
        f"ranks over F_p, p in {sorted(d['ranks'])}" + (f", quotient {d['quotient']}" if d["quotient"] != "none" else ""),
        "all", ("torsion: " + json.dumps(sig)) if sig else "no torsion", "complete", "nonplanar/earlier/fat.jsonl")
for line in open(os.path.join(EARLY, "grp.log")):
    m = re.match(r"(CONTROL )?(.*?) n=(\d+), (.*) p=(\d+): betti (.*)", line.strip())
    if not m: continue
    g = C(4) * 3 + ([(0, 4), (1, 4), (2, 4)] if "degree-3 vertex" in m.group(2) else [])
    add("fat graphs", m.group(2), g, int(m.group(3)), f"ranks over F_{m.group(5)}, quotient {m.group(4)}", "all",
        f"dims {m.group(6)}", "complete", "nonplanar/earlier/grp.log")
K4d = K(4) * 2
for line in open(os.path.join(EARLY, "vexp_bench.log")):
    d = json.loads(line)
    g = {"tripled C4": C(4) * 3, "doubled K4": K4d, "T(2,2,3)": T(2, 2, 3)}[d["case"]]
    if "direct" in d:
        add("fat graphs", d["case"], g, d["n"], f"local Smith form mod {d['p']}^k, full complex", "H2",
            "H2 torsion orders " + str(d["direct"]["H2_torsion"]), "complete", "nonplanar/earlier/vexp_bench.log")
    else:
        ok = not isinstance(d["vexp"]["H2_torsion"], str)
        add("fat graphs", d["case"], g, 5 if d["case"] != "T(2,2,3)" else 6, f"vertex explosion tool at v={d['v']}", "H2",
            str(d["vexp"]["H2_torsion"]), "complete" if ok else "failed (memory)", "nonplanar/earlier/vexp_bench.log")
for line, g, lab in [(l, K4d + [(0, 4), (1, 4), (2, 4)], "doubled K4 + vertex on 3") for l in open(os.path.join(EARLY, "run2.log")) if "vertex on 3" in l] + \
                    [(l, K(5) + K(4), "K5 with doubled K4") for l in open(os.path.join(EARLY, "run2.log")) if "K5 with" in l] + \
                    [(l, K(5) + K(4), "K5 with doubled K4") for l in open(os.path.join(EARLY, "apex.log"))] + \
                    [(l, K4d + [(0, 4)], "doubled K4 + one leaf") for l in open(os.path.join(EARLY, "run3.log"))]:
    p = re.search(r"p=(\d+)", line).group(1)
    q = ", quotient Z/5" if "Z/5" in line else ""
    add("fat graphs", lab, g, 5, f"ranks over F_{p}{q}", "H0-H2", "dims " + line.split("b0..b2")[1].strip(), "complete",
        "nonplanar/earlier/" + ("apex.log" if "mod Z/5" in line else "run2.log" if "run2" else "run3.log"))

# ---------------------------------------------------------------- 5. non-planar: n = 4, K5 doublings, queues, campaign
NP4 = {"K4 doubled (planar control)": K4d, "K33 doubled": Kb(3, 3) * 2, "K5 doubled": K(5) * 2,
       "K5 with a doubled K4": K(5) + K(4)}
for l in open(os.path.join(EARLY, "nonplanar4.jsonl")):
    d = json.loads(l)
    exc = any(any(v) for v in d["excess"].values())
    add("non-planar search", d["name"], NP4[d["name"]], d["n"], "ranks over F_2, F_3, F_5 vs large prime", "all",
        "torsion" if exc else "no torsion", "complete", "nonplanar/earlier/nonplanar4.jsonl")
runs = collections.defaultdict(dict)
for f in ("k5res.jsonl", "k5res2.jsonl"):
    for l in open(os.path.join(EARLY, f)):
        d = json.loads(l); runs[json.dumps(d["doubled"])][d["prime"]] = (d["ranks"], d["mode"])
for f in ("k5run.log", "k5run2.log", "k5run3.log"):
    for l in open(os.path.join(HLOG, f)):
        m = re.match(r"FAILED (\[.*\]) (\d+)( True| False)? -9", l)
        if m:
            runs[json.dumps(ast.literal_eval(m.group(1)))].setdefault("failed", []).append(int(m.group(2)))
for k, v in runs.items():
    dd = json.loads(k); g = K(5) + [tuple(e) for e in dd]
    ps = sorted(p for p in v if p != "failed")
    if not ps:
        add("non-planar search", f"K5 + doubled {dd}", g, 5, "ranks of d1-d3 over F_p, quotient Z/5", "H0-H2",
            "out of memory", "failed (memory)", "history/logs/k5run*.log")
        continue
    agree = len({json.dumps(v[p][0]) for p in ps}) == 1
    mode = v[ps[0]][1]
    status = "complete" if set(ps) >= {2, 3, 1000003} else "partial (small primes ran out of memory)" if v.get("failed") else "partial"
    add("non-planar search", f"K5 + doubled {dd}", g, 5, f"ranks of d1-d3 over F_p, p in {ps}" + (", quotient Z/5" if mode == "Z5" else ""),
        "H0-H2", "no torsion" if agree else "ranks differ", status, "nonplanar/earlier/k5res*.jsonl")
for line in open(os.path.join(EARLY, "k33c.log")):
    p = line.split()[1]
    add("non-planar search", "K33 with a tripled C4 (k33c)", Kb(3, 3) + [(0, 3), (3, 1), (1, 4), (4, 0)] * 2, 5,
        f"ranks of d1-d3 over F_{p}, quotient Z/5", "H0-H2", "dims " + line.split("b0..b2")[1].split("peak")[0].strip(),
        "complete", "nonplanar/earlier/k33c.log")
QG = {"K33": Kb(3, 3), "K5": K(5), "K33+de": Kb(3, 3) + [(0, 3)], "K5+dtri": K(5) + [(0, 1), (0, 2), (1, 2)]}
for f in ("np_queue.log", "np_queue2.log", "np_queue3.log", "np_queue4.log"):
    for l in open(os.path.join(EARLY, f)):
        d = json.loads(l)
        m = re.match(r"(\S+) n=(\d+) (.*)", d["job"])
        res = d["result"]
        ok = not (isinstance(res, str) and res.startswith("FAILED"))
        rs = "no torsion" if ok and (not isinstance(res, dict) or not res.get("torsion")) else str(res)
        add("non-planar search", m.group(1), QG[m.group(1)], int(m.group(2)), f"{m.group(3)}, local mod p={d['p']}",
            "see method", rs if ok else "out of memory", "complete" if ok else "failed (memory)", f"nonplanar/earlier/{f}")
for r in json.load(open(os.path.join(D, "nonplanar", "campaign.json"))):
    o = " ".join(r["options"])
    deg = "H0-H%d" % (int(r["options"][r["options"].index("-d") + 1]) - 1) if "-d" in r["options"] else "all"
    if "dims" in r:
        res = ("torsion " + json.dumps(r["torsion"])) if r.get("torsion") else "no torsion"; st = "complete"
    else:
        res = r.get("result") or ""; st = "failed (memory)" if "memory" in res else "not completed"
    planar_ctrl = planar(r["edges"])
    add("non-planar campaign" if not planar_ctrl else "planar controls (campaign)", r["name"], r["edges"], r["n"],
        f"engine {o}", deg, res, st, "nonplanar/campaign.json")

# ---------------------------------------------------------------- 6. structure computations (analysis tools)
T222 = T(2, 2, 2)
for n, b2 in [(3, 39), (4, 697), (5, 7639)]:
    add("structure of H2", "T(2,2,2)", T222, n, "tori / star x star / Theta3 decomposition (decomp, phi, prod2, cutimg3)",
        "H2", f"b2 = {b2} decomposed", "complete", "nonplanar/earlier/cut5.log, prod_T222.log")
for n, b2 in [(4, 7), (5, 189), (6, 3016)]:
    add("structure of H2", "Theta3", [(0, 1)] * 3, n, "span of tori in H2 (decomp)", "H2", f"b2 = {b2}", "complete", "decomp runs (not logged separately)")
for n, b2 in [(3, 1), (4, 169), (5, 3204)]:
    add("structure of H2", "Theta4", [(0, 1)] * 4, n, "span of tori in H2 (decomp)", "H2", f"b2 = {b2}", "complete", "decomp runs (not logged separately)")
for line in open(os.path.join(EARLY, "split2.log")):
    m = re.match(r"(T\(\d,\d,\d\)) x=2 p=(\d)", line)
    a, b, c = map(int, re.findall(r"\d", m.group(1)))
    add("structure of H2", m.group(1), T(a, b, c), 6, f"split types of torsion (g2split), p={m.group(2)}", "H2",
        line.split(": ", 1)[1].strip()[:80], "complete", "nonplanar/earlier/split2.log")
for f in ("prod_T222.log", "prod_T223.log"):
    for line in open(os.path.join(EARLY, f)):
        if "{" not in line: continue
        g = T(2, 2, 2) if "222" in f else T(2, 2, 3)
        add("structure of H2", "T(2,2,2)" if "222" in f else "T(2,2,3)", g, 6,
            "star x star tori span H2 of the exploded graph (prod2), " + line.split(":")[0].split()[-1], "H2",
            "spans" if '"b2"' in line and json.loads(line[line.index("{"):])["b2"] == json.loads(line[line.index("{"):])["rankP"] else "?",
            "complete", "nonplanar/earlier/" + f)

# ---------------------------------------------------------------- write
fields = list(rows[0].keys())
with open(os.path.join(HERE, "accounting.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fields); w.writeheader(); w.writerows(rows)
with open(os.path.join(HERE, "accounting_graphs.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["graph_class", "first label", "other labels", "V", "E", "planar", "computations", "n values", "edges"])
    for i, (G, inv, lab) in enumerate(_reps):
        key = f"G{i:03d} ({lab})"
        sub = [r for r in rows if r["graph_class"] == key]
        others = sorted({r["graph"] for r in sub} - {lab})
        w.writerow([key, lab, "; ".join(others[:8]) + (" ..." if len(others) > 8 else ""), G.number_of_nodes(),
                    G.number_of_edges(), "yes" if planar(list(G.edges())) else "no", len(sub),
                    " ".join(map(str, sorted({int(r["n"]) for r in sub}))), json.dumps(sorted(list(G.edges())))])
cats = collections.Counter(r["category"] for r in rows)
with open(os.path.join(HERE, "accounting_summary.txt"), "w") as f:
    def P(*a):
        """Print to the screen and to the summary file."""
        print(*a); print(*a, file=f)
    P(f"computations recorded: {len(rows)}")
    P(f"distinct graphs (up to isomorphism): {len(_reps)}")
    P(f"  planar: {sum(1 for G, _, _ in _reps if planar(list(G.edges())))}, non-planar: {sum(1 for G, _, _ in _reps if not planar(list(G.edges())))}")
    P(f"distinct (graph, n) pairs: {len({(r['graph_class'], r['n']) for r in rows})}")
    P("by category:")
    for c, k in cats.most_common():
        sub = [r for r in rows if r["category"] == c]
        P(f"  {c:28s} {k:4d} computations, {len({r['graph_class'] for r in sub}):3d} graphs, "
          f"{sum(r['status'] == 'complete' for r in sub)} complete")
    P("by status:", dict(collections.Counter(r["status"] for r in rows)))
    P("by n:", dict(sorted(collections.Counter(r["n"] for r in rows).items())))
    tor = [r for r in rows if ("torsion" in r["result"] and "no torsion" not in r["result"]) or "Z/" in r["result"]]
    P(f"computations reporting torsion: {len(tor)}")
