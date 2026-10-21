"""minimal_tori.py -- find a small set of tori of the triangle with a leaf whose span contains the class of order 16.

Input: the output of "certificates/reducemod 2 2 K blocks red < graphs/K4dbl5_Hfirst_n5.txt", where blocks holds the
twelve generators z_t of the torsion T of H_2(Conf_5 Gamma_16) (headers "gen id order") followed by the tori of the
triangle with a leaf (headers "tor id 0 x=X A=... i1=.. i2=..", from DUMP=... analysis/decomp).

An oracle (certificates/spanq, kept running in the background) answers "is there an element of order 16 in
T cap (B + span S)?" for a set S of tori: it prints the coordinates of the z_t in C_2 / (B_2 + span S), and the
kernel of T -> that quotient is enumerated as in kernel.py.  An element of T has order 16 exactly when its
coefficient on the order-16 generator is odd, so only that is tested.

CAUTION: the oracle works mod Q = 2^K and can answer "yes" falsely when the images of the tori are divisible by
high powers of 2 (see certificates/spanq).  With K = 8 this happened; use K = 24 and confirm the final set with
exact_check.py.

Search: the property is monotone in S.  Starting from all tori, groups of tori are removed greedily while the
property survives: first whole types (star vertex x and particle set A at the star), then halves, quarters, ... of
what is left, and finally single tori.  The result is a set S that is minimal for removal of single tori.

Usage:  python3 minimal_tori.py spanq_binary reduced_file out.json
"""
import json
import subprocess
import sys
import numpy as np


class Oracle:
    def __init__(self, binary, red):
        """Start spanq on the reduced file and wait until it has read it."""
        self.p = subprocess.Popen([binary, red], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)
        assert self.p.stdout.readline().strip() == "READY"
        self.calls = 0

    def span(self, ids):
        """One 'span' query: the JSON answer of spanq for the listed relation blocks."""
        self.calls += 1
        self.p.stdin.write("span " + " ".join(map(str, ids)) + "\n")
        self.p.stdin.flush()
        return json.loads(self.p.stdout.readline())


def kernel_with_odd16(d, chunk=1 << 20):
    """An element of T cap (B + span S) whose coefficient on the generator of order 16 is odd, or None."""
    diag, cyc = d["diag"], d["cycles"]
    mods = np.array([2 ** a for a in diag], dtype=np.int64)
    o = np.array([c["order_given"] for c in cyc], dtype=np.int64)
    Phi = np.array([c["coords"][:len(diag)] for c in cyc], dtype=np.int64).reshape(len(cyc), len(diag))
    t16 = int(np.argmax(o))
    total = int(np.prod(o))
    for start in range(0, total, chunk):
        idx = np.arange(start, min(total, start + chunk), dtype=np.int64)
        X = np.empty((len(idx), len(o)), dtype=np.int64)
        r = idx.copy()
        for t in range(len(o)):
            X[:, t] = r % o[t]
            r //= o[t]
        X = X[X[:, t16] % 2 == 1]
        if not len(X):
            continue
        ker = np.all((X @ Phi) % mods == 0, axis=1) if len(diag) else np.ones(len(X), dtype=bool)
        if ker.any():
            return X[np.argmax(ker)].tolist()
    return None


def headers(red):
    """Headers of the torus blocks in a reducemod output, in order (their positions are the relation numbers of spanq)."""
    hs = []
    with open(red) as f:
        for line in f:
            if line.startswith("B tor"):
                hs.append(line[6:].strip())
    return hs


def main(binary, red, out):
    """The greedy search; writes the final set and the element found."""
    hs = headers(red)
    typ = []
    for h in hs:
        parts = dict(p.split("=") for p in h.split()[2:])
        typ.append((int(parts["x"]), parts["A"]))
    orc = Oracle(binary, red)
    has = lambda S: kernel_with_odd16(orc.span(S)) is not None
    S = list(range(len(hs)))
    assert has(S), "all tori: no element of order 16"
    log = {"all": len(S)}
    # 1. remove whole types
    types = sorted(set(typ))
    for ty in types:
        T2 = [i for i in S if typ[i] != ty]
        if len(T2) < len(S) and has(T2):
            S = T2
    log["after_types"] = len(S)
    log["types_left"] = sorted(set(typ[i] for i in S))
    print("types left:", log["types_left"], len(S), flush=True)
    # 2. remove chunks of decreasing size, then single tori
    size = max(1, len(S) // 2)
    while True:
        i = 0
        while i < len(S):
            T2 = S[:i] + S[i + size:]
            if has(T2):
                S = T2
            else:
                i += size
        print("chunk", size, "->", len(S), "tori,", orc.calls, "queries", flush=True)
        if size == 1:
            break
        size = max(1, size // 2)
    x = kernel_with_odd16(orc.span(S))
    res = {"tori": S, "headers": [hs[i] for i in S], "element": x, "log": log, "queries": orc.calls}
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res)[:2000])


if __name__ == "__main__":
    main(*sys.argv[1:4])
