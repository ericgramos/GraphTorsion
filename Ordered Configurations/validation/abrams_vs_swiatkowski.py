"""abrams_vs_swiatkowski.py -- compare the Swiatkowski engine with Abrams' discretized model.

Two independent chain complexes compute H_*(Conf_n G):
  * the reduced ordered Swiatkowski complex (engine/oconf, through python/gconf.homology), and
  * Abrams' cube complex D_n(G') of a sufficiently subdivided copy G' of G (python/gconf.abrams_homology,
    pure Python + FLINT, with no code shared with the engine).
Abrams (2000) shows D_n(G') is homotopy equivalent to Conf_n(G) once every edge path between essential
vertices has at least n - 1 edges and every cycle at least n + 1 (we subdivide each edge into m pieces,
chosen per test).  The script also checks the engine's Euler characteristic against Gal's formula.

Usage:  python3 abrams_vs_swiatkowski.py            (the standard list below, a few minutes)
        python3 abrams_vs_swiatkowski.py dt 4 3     (T(2,2,2) = "dt" or T(1,2,3) = "g26", n, m)
Every line must end with agree=True.
"""
import itertools, os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from gconf import homology, abrams_homology, subdivide

K = lambda m: list(itertools.combinations(range(m), 2))
Kb = lambda a, b: [(i, a + j) for i in range(a) for j in range(b)]
if len(sys.argv) == 4:
    tests = {sys.argv[1]: ({"dt": [(0, 1), (0, 1), (0, 2), (0, 2), (1, 2), (1, 2)],
                            "g26": [(0, 1), (0, 2), (0, 2), (1, 2), (1, 2), (1, 2)]}[sys.argv[1]],
                           int(sys.argv[2]), int(sys.argv[3]))}
else:
    tests = {   # name: (edges, n, subdivision m)
        "K5 n=2": (K(5), 2, 1), "K33 n=2": (Kb(3, 3), 2, 1), "K4 n=2": (K(4), 2, 1),
        "K4 n=3": (K(4), 3, 2), "theta3 n=3": ([(0, 1)] * 3, 3, 3),
        "fig8 n=3": ([(0, 0), (0, 0)], 3, 4), "lollipop n=3": ([(0, 1), (1, 1)], 3, 4),
        "H n=3": ([(0, 1), (0, 2), (0, 3), (3, 4), (3, 5)], 3, 2),
        "theta4 n=2": ([(0, 1)] * 4, 2, 3), "K4 n=4": (K(4), 4, 3),
    }
allok = True
for name, (g, n, m) in tests.items():
    t = time.time()
    a = homology(g, n)
    chi_ok = a["chi"] == a["gal"]
    b = abrams_homology(subdivide(g, m), n)
    agree = (b["betti"][:len(a["betti"])] == a["betti"] and all(x == 0 for x in b["betti"][len(a["betti"]):])
             and b["torsion"][:len(a["torsion"])] == a["torsion"] and chi_ok)
    allok &= agree
    print(f"{name}: Swiatkowski betti={a['betti']} torsion={a['torsion']} chi={'ok' if chi_ok else 'MISMATCH'} | "
          f"Abrams betti={b['betti']} torsion={b['torsion']} | agree={agree} ({time.time() - t:.1f}s)", flush=True)
print("ALL OK" if allok else "FAILURES")
