#!/usr/bin/env python3
"""
Regenerate the summary tables from the stored JSON in data/.

The sweeps and the large individual cases write their results to data/*.json as they go.  This
script reads whatever is there and prints the tables, so the stored numbers can be
regenerated without recomputing anything -- and so that a reader who has run only part of the
sweep still gets a coherent, correctly labelled subset.

    python3 experiments/make_tables.py [--latex]

With --latex the tables come out as LaTeX tabular bodies, ready to paste.

The tables:

    1  complete graphs: torsion of H_i(B_k K_n) for i >= 2, by (n, k)
    2  complete graphs: the multiplicity pattern, (Z/2)^m in H_2 against rank H_1 = C(n-1,2)
    3  complete bipartite: every case with odd torsion
    4  complete bipartite: the a = n+1 family and the coker W prediction
    5  the "background" count: 2-torsion multiplicity in H_{k-1} against rank H_1
    6  top-degree vanishing: the trichotomy
    7  everything else stored (planar graphs, controls, Petersen family, large cases)
"""

import argparse
import os
import sys
from math import comb

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import paths                                 # noqa: E402
from gbt.chain import odd_part                        # noqa: E402
from gbt.report import torsion_summary                # noqa: E402


def load(fname):
    """{key: {degree: (free, invariant factors)}} from one data file."""
    raw = paths.load_json(fname)
    out = {}
    for key, rec in raw.items():
        H = rec.get("H") if isinstance(rec, dict) else None
        if H is None:
            continue
        out[key] = {int(d): (v[0], v[1]) for d, v in H.items()}
    return out


class Table(object):
    """A tiny table printer with a plain and a LaTeX mode."""

    def __init__(self, title, headers, latex=False):
        self.title = title
        self.headers = headers
        self.rows = []
        self.latex = latex

    def add(self, *cells):
        self.rows.append([str(c) for c in cells])

    def show(self):
        print()
        print("### %s" % self.title)
        if not self.rows:
            print("(no data -- run the corresponding sweep first)")
            return
        if self.latex:
            print(r"\begin{tabular}{%s}" % ("l" * len(self.headers)))
            print(r"\toprule")
            print(" & ".join(self.headers) + r" \\")
            print(r"\midrule")
            for r in self.rows:
                print(" & ".join(c.replace("_", r"\_") for c in r) + r" \\")
            print(r"\bottomrule")
            print(r"\end{tabular}")
            return
        widths = [max(len(h), max(len(r[i]) for r in self.rows))
                  for i, h in enumerate(self.headers)]
        line = "  ".join(h.ljust(w) for h, w in zip(self.headers, widths))
        print(line)
        print("-" * len(line))
        for r in self.rows:
            print("  ".join(c.ljust(w) for c, w in zip(r, widths)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--latex", action="store_true")
    args = ap.parse_args()
    L = args.latex

    K = load("sweep_complete.json")
    B = load("sweep_bipartite.json")
    big = load("large_cases.json")
    # large_cases.json uses its own keys ("K12-3", "K8-4", "K55", "K56"); rename them into the
    # "K12|3" / "K5,5|4" form of the two sweep files before merging.
    rename = {"K12-3": "K12|3", "K8-4": "K8|4", "K55": "K5,5|4", "K56": "K5,6|4"}
    for old_key, rec in big.items():
        new_key = rename.get(old_key)
        if new_key is None:
            continue
        (B if "," in new_key else K).setdefault(new_key, rec)

    # ------------------------------------------------------------------ 1
    ks = sorted({int(k.split("|")[1]) for k in K})
    ns = sorted({int(k.split("|")[0][1:].split(",")[0]) for k in K if "," not in k})
    t = Table("complete graphs: torsion of H_i(B_k K_n), i >= 2",
              ["K_n"] + ["k=%d" % k for k in ks], latex=L)
    for n in ns:
        row = ["K_%d" % n]
        for k in ks:
            H = K.get("K%d|%d" % (n, k))
            if H is None:
                row.append(".")
                continue
            s = {d: v for d, v in torsion_summary(H).items() if d >= 2}
            row.append("; ".join("H%d: %s" % (d, v) for d, v in s.items()) or "free")
        t.add(*row)
    t.show()

    # ------------------------------------------------------------------ 2
    t = Table("complete graphs at k=3: (Z/2)^m in H_2 vs rank H_1 = C(n-1,2)",
              ["K_n", "rank H_1", "m", "C(n-1,2)", "match"], latex=L)
    for n in ns:
        H = K.get("K%d|3" % n)
        if H is None:
            continue
        m = sum(1 for x in H.get(2, (0, []))[1] if x % 2 == 0)
        t.add("K_%d" % n, H[1][0], m, comb(n - 1, 2),
              "yes" if m == comb(n - 1, 2) else "")
    t.show()

    # ------------------------------------------------------------------ 3
    t = Table("complete bipartite: every case with ODD torsion",
              ["K_{a,b}", "k", "torsion"], latex=L)
    for key in sorted(B, key=lambda s: (len(s), s)):
        H = B[key]
        if not any(odd_part(x) > 1 for d in H for x in H[d][1]):
            continue
        spec, k = key.split("|")
        t.add(spec, k, "   ".join("H_%d: %s" % (d, v)
                                  for d, v in torsion_summary(H).items()))
    t.show()

    # ------------------------------------------------------------------ 4
    t = Table("the a = n+1 family: torsion of H_{k-1} and coker W_{k-1,k}(b)",
              ["K_{a,b}", "k", "b <= 2k-1", "H_{k-1} torsion", "dim H_k",
               "predicted dim H_k"], latex=L)
    for key in sorted(B, key=lambda s: (len(s), s)):
        spec, k = key.split("|")
        k = int(k)
        if "," not in spec:
            continue
        a, b = (int(x) for x in spec.lstrip("K").split(","))
        if a != k + 1:
            continue
        H = B[key]
        pred = max(0, comb(b, k) - comb(b, k - 1))
        t.add("K_{%d,%d}" % (a, b), k, "yes" if b <= 2 * k - 1 else "no",
              torsion_summary(H).get(k - 1, "none"),
              H.get(k, (0, []))[0], pred)
    t.show()

    # ------------------------------------------------------------------ 5
    t = Table("the background: 2-torsion multiplicity in H_{k-1} vs rank H_1",
              ["graph", "k", "rank H_1", "m(Z/2)", "odd", "generic"], latex=L)
    for src, label in ((K, "K"), (B, "K")):
        for key in sorted(src, key=lambda s: (len(s), s)):
            spec, k = key.split("|")
            k = int(k)
            H = src[key]
            if k - 1 not in H or 1 not in H:
                continue
            m = sum(1 for x in H[k - 1][1] if x % 2 == 0)
            odd = ",".join(str(odd_part(x)) for x in H[k - 1][1] if odd_part(x) > 1)
            t.add(spec, k, H[1][0], m, odd or "-",
                  "yes" if (m == H[1][0] and m) else "")
    t.show()

    # ------------------------------------------------------------------ 6
    vanish = paths.load_json("top_degree_vanishing.json")
    t = Table("top-degree vanishing: the trichotomy",
              ["K_{a,b}", "n", "dim H_n", "predicted", "case"], latex=L)
    for key in sorted(vanish, key=lambda s: (len(s), s)):
        rec = vanish[key]
        spec, n = key.split("|")
        n = int(n)
        a, b = (int(x) for x in spec.lstrip("K").split(","))
        case = ("a <= n" if a <= n else
                "a = n+1, b <= 2n-1" if (a == n + 1 and b <= 2 * n - 1) else
                "a = n+1, b >= 2n" if a == n + 1 else "a >= n+2 (conjecture)")
        t.add("K_{%d,%d}" % (a, b), n, rec["H_n"][0],
              "-" if rec["predicted"] is None else rec["predicted"], case)
    t.show()

    # ------------------------------------------------------------------ 7
    for fname, title in (("planar_torsion.json", "planar graphs"),
                         ("linkless_controls.json", "linklessly embeddable controls"),
                         ("petersen_family.json", "the Petersen family"),
                         ("density_interpolation.json", "K_{4,5} to K_9"),
                         ("large_cases.json", "the large individual cases")):
        D = load(fname)
        if not D:
            continue
        t = Table(title, ["case", "k", "homology torsion"], latex=L)
        for key in sorted(D, key=lambda s: (len(s), s)):
            H = D[key]
            parts = key.rsplit("|", 1)
            t.add(parts[0], parts[1] if len(parts) > 1 else "",
                  "   ".join("H_%d: %s" % (d, v)
                             for d, v in torsion_summary(H).items()) or "none")
        t.show()

    print()
    print("sources read from %s:" % paths.DATA)
    for f in (sorted(os.listdir(paths.DATA)) if os.path.isdir(paths.DATA) else []):
        print("   %s" % f)
    return 0


if __name__ == "__main__":
    sys.exit(main())
