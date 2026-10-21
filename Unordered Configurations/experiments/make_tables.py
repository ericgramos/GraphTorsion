#!/usr/bin/env python3
"""
Regenerate the summary tables from the stored JSON in data/.

WHAT IT IS FOR
--------------
The sweeps and the large individual cases write their results to data/*.json as they go.  This
script reads whatever is there and prints the tables, so the stored numbers can be regenerated
without recomputing anything -- and so that a reader who has run only part of a sweep still
gets a coherent, correctly labelled subset.  It computes nothing new: every number is read
from data/, except the predictions (binomial coefficients) in tables 2, 4 and 6.

THE TABLES
----------
Each table starts with a line "### title".

    1  complete graphs: torsion of H_i(B_k K_n), i >= 2, by (n, k)
         one row per n, one column per k; each cell lists the torsion in degrees >= 2
         ("H2: (Z/2)^21"), "free" if there is none, "." if (n, k) was not computed.
    2  complete graphs at k=3: (Z/2)^m in H_2 vs rank H_1 = C(n-1,2)
         columns K_n, rank H_1, m (the number of 2-primary cyclic summands of H_2, i.e. the
         number of even invariant factors), C(n-1,2), and "yes" when m = C(n-1,2).
    3  complete bipartite: every case with ODD torsion
         columns K_{a,b}, k, and all of its torsion by degree.
    4  the a = n+1 family: torsion of H_{k-1} and coker W_{k-1,k}(b)
         for every stored K_{k+1,b}: whether b <= 2k-1, the torsion of H_{k-1}, the computed
         rank of H_k, and the predicted rank max(0, C(b,k) - C(b,k-1)).
    5  the background: 2-torsion multiplicity in H_{k-1} vs rank H_1
         for every complete and complete bipartite case: rank H_1, the number m of
         2-primary summands of H_{k-1}, the odd parts of its invariant factors, and
         "yes" in the last column when m = rank H_1 > 0.
    6  top-degree vanishing: the trichotomy
         from data/top_degree_vanishing.json: dim H_n(B_n K_{a,b}), the prediction ("-"
         where only the conjecture applies), and which case of the trichotomy (a, b, n) is in.
    7  one table each for planar_torsion.json, linkless_controls.json, petersen_family.json,
         density_interpolation.json and large_cases.json: case, k, and all torsion.

The run ends with the list of files found in data/.  An empty source prints "(no data -- run
the corresponding sweep first)" under its title.  Excerpt of table 2:

    ### complete graphs at k=3: (Z/2)^m in H_2 vs rank H_1 = C(n-1,2)
    K_n   rank H_1  m   C(n-1,2)  match
    -----------------------------------
    K_7   15        1   15
    K_8   21        21  21        yes

The large cases are merged into tables 1-5 under the sweep keys ("K12-3" becomes "K12|3",
"K55" becomes "K5,5|4"); a sweep record with the same key takes precedence.

USAGE
-----
    python3 experiments/make_tables.py [--latex]

    --latex   print each table as a LaTeX tabular body (booktabs rules, one "l" column per
              header, underscores escaped) instead of aligned plain text.

Instant; reads data/ only and writes nothing.  Exit code: always 0.
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
    """{key: {degree: (free, invariant factors)}} from one data file.

    Records without a full homology dict "H" (e.g. single-degree results) are skipped, and a
    missing file gives {}.
    """
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
        """title: printed after "###"; headers: column names; latex: LaTeX output if True."""
        self.title = title
        self.headers = headers
        self.rows = []
        self.latex = latex

    def add(self, *cells):
        """Append one row; every cell is converted with str()."""
        self.rows.append([str(c) for c in cells])

    def show(self):
        """Print the table (or the no-data note), padding columns to their widest entry."""
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
    """Load the stored results and print tables 1-7."""
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
