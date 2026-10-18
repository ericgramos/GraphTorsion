"""fisharp.py -- the binomial-inversion test (FI#-type decomposition).

If H_i(Conf_n G) were a sum of induced representations as for an FI#-module (Church-Ellenberg-Farb 2015),
H_i(Conf_n) = (+)_k Ind_{S_k x S_{n-k}}^{S_n} (W_k (x) trivial), then
    b_i(Conf_n) = sum_k C(n, k) dim W_k   and therefore   w_k = sum_j (-1)^(k-j) C(k, j) b_i(Conf_j)
are non-negative integers (w_k = dim W_k).  This script prints the sequences b_i(n) and w_k for a few
graphs.  In all cases tested the w_k were non-negative, e.g. b_1(Conf_n K_4) = 3n + C(n, 2).

Usage:  python3 fisharp.py      (uses explode.betti, i.e. the engine through gconf.homology)
"""
import itertools, math
from explode import betti

K = lambda m: list(itertools.combinations(range(m), 2))
graphs = {   # name: (edges, largest n to compute)
    "Y (tree)": ([(0, 1), (0, 2), (0, 3)], 6),
    "theta3+leaf": ([(0, 1), (0, 1), (0, 1), (0, 2)], 5),
    "K4+leaf": (K(4) + [(0, 4)], 5),
    "sun3": ([(0, 1), (1, 2), (2, 0), (0, 3), (1, 4), (2, 5)], 5),
    "K5+leaf": (K(5) + [(0, 5)], 4),
    "K4 (no leaf)": (K(4), 5),
    "theta3 (no leaf)": ([(0, 1)] * 3, 6),
}
for name, (g, N) in graphs.items():
    bs = [betti(g, n) for n in range(N + 1)]          # bs[n] = Betti numbers of Conf_n
    D = max(len(b) for b in bs)
    print(name)
    for i in range(D):
        b = [x[i] if i < len(x) else 0 for x in bs]
        w = [sum((-1) ** (n - k) * math.comb(n, k) * b[k] for k in range(n + 1)) for n in range(N + 1)]
        print(f"   H_{i}: b(n)={b}  w(k)={w}")
