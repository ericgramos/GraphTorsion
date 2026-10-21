#!/usr/bin/env python3
"""
sweep_template.py -- a template for your own experiment: run a family of graphs, store the
answers, resume after an interruption, and summarise.

WHAT IT IS FOR
--------------
Every script in experiments/ follows the same pattern, and this file is that pattern in its
simplest form, written to be copied and edited.  As shipped it computes H_*(B_k W_n) for the
wheels W_n (an n-cycle plus a hub joined to every rim vertex) for a range of n and k, and
reports which ones have torsion above H_1.  To study a different family, change the function
`jobs()`; everything else can stay as it is.

THE PATTERN
-----------
  1. Make a list of jobs (name, graph, k).
  2. Sort them by size, using the closed formula gbt.reduced_cell_count, so that the cheap
     answers arrive first and a cap on the size can be applied BEFORE anything is built.
  3. Load the results stored so far.  Skip every job that is already there.
  4. For each remaining job: compute, print one line, and save the whole results file at once,
     so that an interruption (Ctrl-C, a closed laptop, running out of memory) loses at most the
     job in progress.  Rerunning the script continues where it stopped.
  5. Print a summary from the stored results.

USAGE
-----
    python3 examples/sweep_template.py [--nmax 8] [--kmax 4] [--cap 300000]
                                       [--out examples/output/wheel_sweep.json]

    --nmax N    largest wheel W_N (default 8; the smallest is W_3 = K_4)
    --kmax K    largest number of particles (default 4; the smallest is 2)
    --cap C     skip any job whose fully reduced complex has more than C cells (default
                300,000, which keeps every job to a few seconds and well under 1 GB)
    --out PATH  the results file (default examples/output/wheel_sweep.json).  It is NOT in
                data/, so running the template never touches the stored results of the
                project.  Delete the file to start again from nothing.

OUTPUT
------
One line per job, then a summary, for example

    W_5    k=4      10,995 cells  H_1=Z^6  H_2=Z^40  H_3=0  H_4=0          torsion: none   [0.1s]
    W_6    k=4      23,473 cells  H_1=Z^7  H_2=Z^63  H_3=0  H_4=0          torsion: none   [0.2s]
    ...
    summary: 18 jobs stored; torsion above H_1 in: none

The whole default run takes about two seconds.  Run it a second time and it prints only the
summary, because every job is already stored.

(Wheels are planar, so H_1 is torsion-free by Ko--Park; whether a wheel ever has torsion in
H_2 is exactly the kind of question this template is for.)

The results file is JSON: {"W_5|4": {"N": 6, "M": 10, "k": 4, "cells": 10995,
"H": {"0": [1, []], "1": [6, []], ...}}, ...}, where each H entry is [free rank, [invariant
factors > 1]], the same format as the files in data/.
"""

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gbt import fast, reduced_cell_count, wheel                             # noqa: E402
from gbt.report import has_torsion, homology_line, torsion_line             # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def jobs(nmax, kmax):
    """The list of (name, graph, k) to compute.  EDIT THIS to study another family.

    The name is used as part of the key in the results file, so it must identify the graph
    uniquely; "W_5" for the wheel with a 5-cycle rim.
    """
    out = []
    for n in range(3, nmax + 1):
        for k in range(2, kmax + 1):
            out.append(("W_%d" % n, wheel(n), k))
    return out


def load(path):
    """The results stored so far: a dict {key: record}, or {} if the file does not exist."""
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return json.load(fh)


def save(path, results):
    """Write all results to `path`, creating the folder if needed.  Writing to a temporary
    file and renaming it means an interruption mid-write cannot corrupt the stored file."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(results, fh, indent=1, sort_keys=True)
    os.replace(tmp, path)


def main():
    """Run the sweep: sort by size, skip stored jobs, compute and save one job at a time."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--kmax", type=int, default=4)
    ap.add_argument("--cap", type=int, default=300_000)
    ap.add_argument("--out", default=os.path.join(HERE, "output", "wheel_sweep.json"))
    args = ap.parse_args()

    results = load(args.out)
    todo = []
    for name, G, k in jobs(args.nmax, args.kmax):
        size = reduced_cell_count(G, k)               # a closed formula: costs nothing
        if size <= args.cap:
            todo.append((size, name, G, k))
    todo.sort(key=lambda job: job[0])                 # cheapest first

    for size, name, G, k in todo:
        key = "%s|%d" % (name, k)
        if key in results:
            continue                                  # already computed in an earlier run
        t0 = time.time()
        H, _, _ = fast.compute(G, k)                  # {degree: (free rank, [factors > 1])}
        results[key] = {"N": G.N, "M": G.M, "k": k, "cells": size,
                        "H": {str(d): [H[d][0], list(H[d][1])] for d in H}}
        save(args.out, results)                       # save after EVERY job
        print("%-6s k=%-2d %10s cells  %-40s torsion: %s   [%.1fs]"
              % (name, k, format(size, ","), homology_line(H), torsion_line(H),
                 time.time() - t0), flush=True)

    # The summary is computed from the stored file, so it covers earlier runs too.
    carriers = []
    for key, rec in sorted(results.items()):
        H = {int(d): (v[0], v[1]) for d, v in rec["H"].items()}
        if has_torsion(H, min_degree=2):
            carriers.append(key)
    print("summary: %d jobs stored; torsion above H_1 in: %s"
          % (len(results), ", ".join(carriers) if carriers else "none"))
    print("results file: %s" % os.path.relpath(args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
