"""kernel.py -- which part of a finite abelian 2-group T maps to zero?

Input: the JSON written by certificates/orderq for a list of cycles z_t that generate the 2-primary torsion
T = (+)_t Z/o_t of H_2(Conf_n Gamma) (the generator file of "certq ... -1"), computed modulo the image of a
subgraph H (and possibly further relations).  orderq gives, for every z_t, its coordinates phi(z_t) in the
torsion of the quotient module M, which is (+)_j Z/2^a_j.  The map

        phi : T -> (+)_j Z/2^a_j,     x = sum x_t z_t  |->  sum x_t phi(z_t),

is well defined (each z_t has order o_t in T, so o_t phi(z_t) = 0).  Its kernel is T cap A, where A is the
image of H_2(Conf_n H) (plus the relations): exactly the torsion elements that come from H.

T is small here (|T| = 2^23 for the graph with Z/16), so the kernel is found by enumerating all of T in
chunks with numpy.  The structure of the kernel K is read off from the number of its elements of each order:
if K = (+) Z/2^e_i, the number of x in K with 2^i x = 0 is 2^(sum_e min(e, i)), so the number r_i of summands
with e >= i is log2 |K[2^i]| - log2 |K[2^(i-1)]|.

Usage:  python3 kernel.py orderq_output.json [...]
"""
import json
import sys
from collections import Counter

import numpy as np


def structure(order_counts):
    """order_counts: {2^e: number of elements of order exactly 2^e}.  Returns {2^e: number of Z/2^e summands}."""
    exps = sorted(int(round(np.log2(k))) for k in order_counts)
    top = max(exps) if exps else 0
    by_exp = {int(round(np.log2(k))): v for k, v in order_counts.items()}
    cum, logs = 0, []
    for i in range(top + 1):                        # logs[i] = log2 |K[2^i]|
        cum += by_exp.get(i, 0)
        logs.append(int(round(np.log2(cum))))
    r = [logs[i] - logs[i - 1] for i in range(1, top + 1)] + [0]   # r[i-1] = summands with exponent >= i
    return {2 ** i: r[i - 1] - r[i] for i in range(1, top + 1) if r[i - 1] - r[i]}


def analyse(fn, chunk=1 << 20):
    """Returns (orderq data, {order: count} for the kernel, |kernel|, |T|)."""
    d = json.load(open(fn))
    diag, cyc = d["diag"], d["cycles"]
    mods = np.array([2 ** a for a in diag], dtype=np.int64)
    o = np.array([c["order_given"] for c in cyc], dtype=np.int64)
    Phi = np.array([c["coords"][:len(diag)] for c in cyc], dtype=np.int64).reshape(len(cyc), len(diag))
    total = int(np.prod(o))
    counts = Counter()
    for start in range(0, total, chunk):
        idx = np.arange(start, min(total, start + chunk), dtype=np.int64)
        X = np.empty((len(idx), len(o)), dtype=np.int64)
        r = idx.copy()
        for t in range(len(o)):                     # mixed-radix digits: x_t in [0, o_t)
            X[:, t] = r % o[t]
            r //= o[t]
        if len(diag):
            ker = np.all((X @ Phi) % mods == 0, axis=1)
        else:
            ker = np.ones(len(idx), dtype=bool)
        K = X[ker]
        if len(K):
            orders = np.max(o // np.gcd(K, o), axis=1)   # order of sum x_t z_t in T
            for v, c in zip(*np.unique(orders, return_counts=True)):
                counts[int(v)] += int(c)
    return d, counts, sum(counts.values()), total


if __name__ == "__main__":
    for fn in sys.argv[1:]:
        d, oc, ks, tot = analyse(fn)
        print(fn.split("/")[-1],
              "| torsion of M below 2^K:", dict(Counter(2 ** a for a in d["diag"])),
              "| |T cap A| =", ks, "of |T| =", tot,
              "| T cap A =", structure(oc),
              "| orders of the generator images:", [c["order_in_M"] for c in d["cycles"]])
