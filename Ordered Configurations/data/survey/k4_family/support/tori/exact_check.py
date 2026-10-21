"""exact_check.py -- does a given set of tori contain, exactly, a combination of order 16?

For tori T_1..T_m (integral 2-cycles of the triangle with a leaf H, pushed into Gamma_16) we want an integer vector a
such that y = sum a_i T_i is TORSION in H_2(Conf_5 Gamma_16) and has order 16.

  1. Torsion.  y is torsion iff its image in H_2(Gamma_16; Q) vanishes.  certificates/reducemod with a large prime p
     and K = 1 writes every T_i as a vector modulo boundaries over F_p (for p large there is no p-torsion, so this is
     H_2 with F_p coefficients, and no residual rows appear).  The kernel of the m vectors is computed over F_p for
     two primes; its reduced row echelon basis is lifted by rational reconstruction and must agree for both primes.
     When every entry of that basis is an integer (as in all our cases), the basis vectors form a Z-basis of the
     integral relation lattice L = {a in Z^m : sum a_i T_i is torsion} (they are unit vectors in the free variables).
  2. Order 16.  f is the cocycle mod 16 written by certq for the generator of order 16.  For torsion classes, f is odd
     exactly on the elements of order 16 (f is even on the summands of order 2 and 8, see certq).  So L contains an
     element of order 16 iff f(sum a_i T_i) is odd for some basis vector a of L.
  The answer is then certified independently (make_certificate.py and certificates/verify_cert.py), so the modular
  steps above only have to find a candidate.

Usage:  python3 exact_check.py tori_blocks generators_with_F graph_file id [id ...]   (ids as in the tori file)
        tori_blocks: the DUMP file of analysis/decomp; generators_with_F: the certq generator file (index 11 = order 16).
Prints the dimension of L and one shortest basis vector with odd f (or none) as one JSON line; with the environment
variable EXACT_OUT=file the same JSON is also written to that file.
"""
import json
import math
import os
import subprocess
import sys
import tempfile
from fractions import Fraction

import flint

HERE = os.path.dirname(os.path.abspath(__file__))
REDUCEMOD = os.path.join(HERE, "..", "..", "..", "..", "..", "certificates", "reducemod")
PRIMES = (1000003, 998244353)


def extract(tfn, ids, out):
    """Copy the requested tori out of the DUMP file into a block file for reducemod."""
    want, keep = set(map(str, ids)), False
    with open(out, "w") as o:
        for line in open(tfn):
            if line.startswith("Z "):
                keep = line.split()[1] in want
                if keep:
                    o.write("Z tor " + line[2:])
                continue
            if keep:
                o.write(line)


def read_red(fn):
    """Read a reducemod output: block headers and reduced vectors (no residual rows are expected for a large prime)."""
    hdr, vec = [], []
    for line in open(fn):
        if line.startswith("LIVE"):
            continue
        if line.startswith("R"):
            raise SystemExit("residual rows present: the prime divides the torsion?")
        if line.startswith("B "):
            hdr.append(line.split()[2])
            continue
        vec.append({int(a): int(b) for a, b in (t.split(":") for t in line.split())})
    return hdr, vec


def ratrec(a, p):
    """Rational reconstruction of a residue mod p (numerator and denominator up to sqrt(p/2))."""
    B = math.isqrt(p // 2)
    r0, r1, s0, s1 = p, a % p, 0, 1
    while r1 > B:
        q = r0 // r1
        r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
    if s1 == 0 or abs(s1) > B:
        return None
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)


def kernel(vec, p):
    """Kernel of the m reduced vectors over F_p: the rank and the reduced row echelon kernel basis, lifted to Q."""
    m = len(vec)
    cols = sorted({c for v in vec for c in v})
    ci = {c: i for i, c in enumerate(cols)}
    M = flint.nmod_mat(max(1, len(cols)), m, p)
    for j, v in enumerate(vec):
        for c, x in v.items():
            M[ci[c], j] = x
    E, r = M.rref()
    piv, j = [], 0
    for i in range(r):
        while int(E[i, j]) == 0:
            j += 1
        piv.append(j)
        j += 1
    K = []
    for f in (j for j in range(m) if j not in piv):
        v = [Fraction(0)] * m
        v[f] = Fraction(1)
        for i, pj in enumerate(piv):
            e = int(E[i, f])
            if e:
                q = ratrec(e, p)
                if q is None:
                    raise SystemExit("rational reconstruction failed")
                v[pj] = -q
        K.append(v)
    return r, K


def fvalues(gfn, tfile):
    """Values of the mod-16 cocycle of the order-16 generator on each torus."""
    F, cur = {}, None
    for line in open(gfn):
        if line.startswith("Z "):
            cur = int(line.split()[1])
            continue
        if cur == 11 and line.startswith("F "):
            a, cell = line[2:].rstrip("\n").split(" ", 1)
            F[cell] = int(a)
    E = F and len(next(iter(F)).split("] E[")[1][:-1].split("|"))

    def pad(cell):
        """Write a cell for the full edge list (empty slots for the edges outside H)."""
        v, e = cell.split("] E[")
        slots = e[:-1].split("|")
        return v + "] E[" + "|".join(slots + [""] * (E - len(slots))) + "]"
    vals, cur = {}, None
    for line in open(tfile):
        if line.startswith("Z "):
            cur = line.split()[2]
            vals[cur] = 0
            continue
        a, cell = line.rstrip("\n").split(" ", 1)
        vals[cur] += int(a) * F.get(pad(cell), 0)
    return vals


def main(tfn, gfn, graph, ids):
    """The exact test; prints one JSON line."""
    tmp = tempfile.mkdtemp()
    tf = os.path.join(tmp, "tori.txt")
    extract(tfn, ids, tf)
    Ks = []
    for p in PRIMES:
        red = os.path.join(tmp, f"red{p}.txt")
        with open(graph) as g:
            subprocess.run([REDUCEMOD, "2", str(p), "1", tf, red], stdin=g, stderr=subprocess.DEVNULL, check=True)
        hdr, vec = read_red(red)
        r, K = kernel(vec, p)
        Ks.append((hdr, r, K))
    (hdr, r, K), (_, r2, K2) = Ks
    assert r == r2 and K == K2, "the two primes disagree"
    assert all(x.denominator == 1 for v in K for x in v), "non-integral kernel basis"
    fv = fvalues(gfn, tf)
    f = [fv[h] for h in hdr]
    odd = []
    for v in K:
        s = sum(int(x) * f[j] for j, x in enumerate(v) if x)
        if s % 2:
            odd.append((sum(1 for x in v if x), {hdr[j]: int(x) for j, x in enumerate(v) if x}, s % 16))
    odd.sort(key=lambda t: t[0])
    res = {"tori": len(hdr), "rank_of_images": r, "kernel_dim": len(K), "odd_vectors": len(odd),
           "best": odd[0] if odd else None}
    print(json.dumps(res))
    if os.environ.get("EXACT_OUT"):
        json.dump(res, open(os.environ["EXACT_OUT"], "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
