"""
gbt -- graph braid torsion: exact integral homology of unordered configuration
spaces of graphs.

For a finite simple graph Gamma and an integer k >= 1, let

    B_k(Gamma) = UConf_k(Gamma) = (Gamma^k - Delta) / Sigma_k

be the unordered configuration space of k points on Gamma.  This package computes

    H_*(B_k Gamma; Z)

exactly: no floating point, no modular reconstruction, no probabilistic rank
estimates.  Each homology group is returned as a free rank together with the list
of its invariant factors, so that torsion is reported rather than inferred.

WHAT IS IN HERE
---------------
    graphs        the Graph class, standard families, named small graphs
    swiatkowski   the Swiatkowski chain complex, with optional partial reduction
    chain         ChainComplex: d^2 = 0 test, discrete Morse reduction, homology
    smith         exact Smith normal form over Z
    linalg        integer kernel lattices, cokernel torsion, ranks over F_p
    report        one printing format for homology, shared by every experiment
    fast          integer-indexed pipeline, ~100x faster than `chain`
    lean          memory-lean pipeline for complexes of a few million cells
    discretised   Abrams' discretised cube complex -- an independent model
    sizes         exact cell counts, for ordering and budgeting jobs
    flipchains    the flip chains z(e), X_{B_0}, Y_T for K_{n+1,b}

TYPICAL USE
-----------
    >>> from gbt import complete, fast, fmt
    >>> H, raw, red = fast.compute(complete(7), 3)
    >>> fmt(*H[2])
    'Z^350 + Z/2'

See the repository README for which pipeline to use at which size, and
`experiments/` for one runnable script per computation.
"""

from .graphs import (                                            # noqa: F401
    Graph,
    complete,
    complete_bip,
    from_edges,
    from_networkx,
    subdivide,
    disjoint_union,
    tetrahedron,
    octahedron,
    cube,
    dodecahedron,
    icosahedron,
    platonic_solids,
    antiprism,
    prism_graph,
    prism,
    theta,
    wheel,
    petersen,
    petersen_family,
)
from .swiatkowski import swiatkowski                             # noqa: F401
from .chain import ChainComplex, fmt, primary                    # noqa: F401
from .smith import smith                                         # noqa: F401
from .sizes import cell_count, reduced_cell_count                # noqa: F401
from .discretised import abrams                                  # noqa: F401
from . import fast, lean, linalg, discretised, paths, report, flipchains  # noqa: F401

__all__ = [
    "Graph", "complete", "complete_bip", "from_edges", "from_networkx",
    "subdivide", "disjoint_union", "tetrahedron", "octahedron", "cube",
    "dodecahedron", "icosahedron", "platonic_solids", "antiprism", "prism_graph",
    "prism", "theta",
    "wheel", "petersen", "petersen_family",
    "swiatkowski", "ChainComplex", "fmt", "primary", "smith",
    "cell_count", "reduced_cell_count",
    "abrams", "fast", "lean", "linalg", "discretised", "paths", "report",
    "flipchains",
]

__version__ = "1.0"
