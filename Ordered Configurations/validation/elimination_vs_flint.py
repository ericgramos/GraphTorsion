"""elimination_vs_flint.py -- the engine's exact sparse elimination against FLINT's Smith normal form.

300 random sparse integer matrices (sizes up to 40 x 40, several densities and entry distributions) are
eliminated by the engine's generic matrix mode (oconf -m: unit pivots, then exact SNF of the residual by
FLINT inside python/gconf.sparse_snf) and by a dense FLINT Smith normal form.  Rank and the multiset of
nonunit elementary divisors must agree.  Prints "mismatches = 0".
"""
import os, random, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
import flint
from gconf import sparse_snf

random.seed(1)
bad = 0
for trial in range(300):
    r, c = random.randint(1, 40), random.randint(1, 40)
    dens = random.choice([0.05, 0.1, 0.3])
    vals = random.choice([[1, -1], [1, -1, 2], [1, -1, 2, -3, 6], [2, 4, -2, 3]])
    M = [[random.choice(vals) if random.random() < dens else 0 for _ in range(c)] for _ in range(r)]
    trip = [(i, j, M[i][j]) for i in range(r) for j in range(c) if M[i][j]]
    S = flint.fmpz_mat(M).snf()
    d = [abs(int(S[i, i])) for i in range(min(r, c)) if S[i, i] != 0]
    rk, dv = sparse_snf(r, c, trip)
    if rk != len(d) or sorted(dv) != sorted(x for x in d if x != 1):
        bad += 1
        print("MISMATCH", r, c, d, rk, dv)
print("elimination tests: mismatches =", bad)
