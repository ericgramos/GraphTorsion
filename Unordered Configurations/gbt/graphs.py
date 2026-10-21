"""
Finite simple graphs, the families used in the experiments, and a few named small graphs.

A graph here is a `Graph` object: vertices are 0 .. N-1, and the INDEX of an edge in
`G.edges` is that edge's name everywhere else in the package.  Both conventions matter,
because the Swiatkowski differential depends on the vertex order (Koszul signs) and on
the choice of a basepoint edge at each vertex (partial reduction).

Vertex order conventions used by the bipartite computations:

    complete_bip(a, b)   puts the a-side first, so  A = {0, ..., a-1}  and
                         B = {a, ..., a+b-1}.  Every vertex of A precedes every vertex
                         of B, which is what the sign bookkeeping in gbt.flipchains assumes.

Only two functions in this module need networkx -- `from_networkx` and
`petersen_family` -- and they import it lazily, so everything else (and the whole core of
the package) runs on a bare Python 3 installation.

HOW TO ENTER YOUR OWN GRAPH
---------------------------
    >>> from gbt import Graph, from_edges
    >>> G = Graph(4, [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)])   # vertices 0..3
    >>> G.N, G.M, G.edges[4], G.inc[0]
    (4, 5, (0, 2), [0, 3, 4])
    >>> H = from_edges([("a", "b"), ("b", "c"), ("c", "a")])     # any labels; relabelled
    >>> H.edges
    [(0, 1), (1, 2), (0, 2)]

Edges are stored with the smaller endpoint first, in the order given, and that order is
their index.  Reordering the edge list gives an isomorphic complex (different basepoint
edges, different generator names) and therefore the same homology.

Families and named graphs: complete, complete_bip, octahedron, tetrahedron, dodecahedron,
icosahedron, platonic_solids, antiprism, prism_graph, prism, theta, wheel, petersen,
petersen_family; operations: from_edges, from_networkx, disjoint_union, subdivide,
Graph.delete_edge.
"""

import itertools


class Graph:
    """A finite simple graph on the vertex set {0, 1, ..., N-1}.

    Attributes
    ----------
    N : int
        number of vertices.
    edges : list of 2-tuples
        the edges, each stored sorted; `edges[i]` is the edge with index i.
    M : int
        number of edges.
    inc : list of lists
        `inc[v]` is the list of indices of the edges incident to v, in increasing order.

    Multi-edges and loops are rejected: the Swiatkowski model as implemented here is for
    simple graphs.  (A loop or multi-edge can always be subdivided away; see `subdivide`.)
    """

    def __init__(self, N, edges):
        """N vertices 0..N-1 and an iterable of 2-element edges (pairs of vertices).

        Raises ValueError on a repeated edge, a loop, or an endpoint outside 0..N-1.
        Vertices with no edges are allowed (but see the pitfall in `gbt.fast`).
        """
        self.N = N
        self.edges = [tuple(sorted(e)) for e in edges]
        if len(set(self.edges)) != len(self.edges):
            raise ValueError("repeated edge: this package is for simple graphs")
        self.M = len(self.edges)
        self.inc = [[] for _ in range(N)]
        for i, (u, v) in enumerate(self.edges):
            if u == v:
                raise ValueError("loop at vertex %d" % u)
            if not (0 <= u < N and 0 <= v < N):
                raise ValueError("edge %r outside 0..%d" % ((u, v), N - 1))
            self.inc[u].append(i)
            self.inc[v].append(i)

    # -- basic queries ------------------------------------------------------
    def deg(self, v):
        """Degree of v."""
        return len(self.inc[v])

    def other(self, v, ei):
        """The endpoint of the edge with index `ei` other than `v`."""
        u, w = self.edges[ei]
        return w if u == v else u

    def nbrs(self, v):
        """Neighbours of v, in the order of the incident edge indices."""
        return [self.other(v, e) for e in self.inc[v]]

    def edge_index(self, u, v):
        """Index of the edge {u, v}; raises KeyError if there is no such edge."""
        key = tuple(sorted((u, v)))
        for i in self.inc[u]:
            if self.edges[i] == key:
                return i
        raise KeyError("no edge %r" % (key,))

    def is_connected(self):
        """True if the graph is connected (depth-first search from vertex 0); the empty
        graph counts as connected."""
        if self.N == 0:
            return True
        seen = {0}
        stack = [0]
        while stack:
            v = stack.pop()
            for w in self.nbrs(v):
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        return len(seen) == self.N

    def delete_edge(self, ei):
        """A copy of this graph with the edge of index `ei` removed (same vertex set)."""
        return Graph(self.N, [e for j, e in enumerate(self.edges) if j != ei])

    def __repr__(self):
        """Short description, e.g. 'Graph(N=6, M=12)'."""
        return "Graph(N=%d, M=%d)" % (self.N, self.M)


# =========================================================================
# Families
# =========================================================================

def complete(n):
    """The complete graph K_n on vertices 0..n-1, edges in lexicographic order."""
    return Graph(n, itertools.combinations(range(n), 2))


def complete_bip(a, b):
    """The complete bipartite graph K_{a,b}, with the a-side first.

    A = {0, ..., a-1} and B = {a, ..., a+b-1}; every vertex of A precedes every vertex
    of B.  The edge list is in lexicographic order of (A-endpoint, B-endpoint).
    """
    return Graph(a + b, [(i, a + j) for i in range(a) for j in range(b)])


def from_edges(edges):
    """Build a graph from an arbitrary edge list, relabelling the vertices to 0..N-1.

    The relabelling is order preserving: the sorted distinct endpoints are sent to
    0, 1, 2, ... in order, so a bipartite edge list whose parts are already separated
    stays separated.
    """
    vs = sorted({x for e in edges for x in e})
    idx = {v: i for i, v in enumerate(vs)}
    return Graph(len(vs), [(idx[u], idx[v]) for u, v in edges])


def from_networkx(Gx):
    """Convert a networkx graph (any hashable labels) to a `Graph`.

    The nodes are relabelled 0..N-1 in networkx's node order, and the edge order is
    networkx's edge order.  Multi-edges are collapsed (the input is converted to a simple
    nx.Graph first); a self-loop raises ValueError.  Requires networkx.
    """
    import networkx as nx
    Gx = nx.convert_node_labels_to_integers(nx.Graph(Gx))
    return Graph(Gx.number_of_nodes(), list(Gx.edges()))


def disjoint_union(G, H):
    """The disjoint union G + H, with H's vertices shifted up by G.N.

    The edges of G come first (same indices), then those of H.  By the Kunneth formula the
    homology of B_k(G + H) is a sum over k_1 + k_2 = k of tensor and Tor products of the
    homology of B_{k_1} G and B_{k_2} H.
    """
    return Graph(G.N + H.N,
                 list(G.edges) + [(u + G.N, v + G.N) for (u, v) in H.edges])


def subdivide(G, t):
    """Replace every edge of G by a path of t edges.

    Needed for Abrams' discretised model, which requires enough room on each edge;
    irrelevant to the Swiatkowski model, which never needs subdivision.  The original
    vertices keep their numbers 0..N-1; the new subdivision vertices are numbered N, N+1,
    ... in order along each edge, edge by edge.  t = 1 returns a copy of G.
    """
    if t < 1:
        raise ValueError("t must be at least 1")
    edges = []
    nxt = G.N
    for (u, v) in G.edges:
        if t == 1:
            edges.append((u, v))
            continue
        prev = u
        for _ in range(t - 1):
            edges.append((prev, nxt))
            prev = nxt
            nxt += 1
        edges.append((prev, v))
    return Graph(nxt, edges)


# =========================================================================
# Named small graphs that recur in the experiments
# =========================================================================

#: edge list of the octahedron K_{2,2,2} with the antipodal pairs {0,1}, {2,3}, {4,5}
OCTAHEDRON_EDGES = [(0, 2), (0, 3), (0, 4), (0, 5),
                    (1, 2), (1, 3), (1, 4), (1, 5),
                    (2, 4), (2, 5), (3, 4), (3, 5)]


def octahedron():
    """The octahedron K_{2,2,2}: planar, 6 vertices, 12 edges.

    Vertices 0..5 with the antipodal (non-adjacent) pairs {0,1}, {2,3}, {4,5}; the edge
    list is OCTAHEDRON_EDGES.  The carrier of 2-torsion in H_2(B_4) with the fewest
    vertices among all graphs on at most seven vertices (see
    experiments/sweep_small_graphs.py).
    """
    return Graph(6, OCTAHEDRON_EDGES)


def tetrahedron():
    """The tetrahedron, i.e. K_4 (this simply returns complete(4)): planar, 4 vertices,
    6 edges, 3-regular."""
    return complete(4)


#: edge list of the dodecahedron (3-regular, 20 vertices, 30 edges)
DODECAHEDRON_EDGES = [
    (0, 1), (0, 10), (0, 19), (1, 2), (1, 8), (2, 3), (2, 6), (3, 4), (3, 19),
    (4, 5), (4, 17), (5, 6), (5, 15), (6, 7), (7, 8), (7, 14), (8, 9), (9, 10),
    (9, 13), (10, 11), (11, 12), (11, 18), (12, 13), (12, 16), (13, 14),
    (14, 15), (15, 16), (16, 17), (17, 18), (18, 19)]

#: edge list of the icosahedron (5-regular, 12 vertices, 30 edges)
ICOSAHEDRON_EDGES = [
    (0, 1), (0, 5), (0, 6), (0, 7), (0, 10), (1, 2), (1, 5), (1, 7), (1, 11),
    (2, 3), (2, 7), (2, 8), (2, 11), (3, 4), (3, 8), (3, 9), (3, 11), (4, 5),
    (4, 9), (4, 10), (4, 11), (5, 10), (5, 11), (6, 7), (6, 8), (6, 9), (6, 10),
    (7, 8), (8, 9), (9, 10)]


def dodecahedron():
    """The dodecahedron: planar, 20 vertices, 30 edges, 3-regular."""
    return Graph(20, DODECAHEDRON_EDGES)


def icosahedron():
    """The icosahedron: planar, 12 vertices, 30 edges, 5-regular."""
    return Graph(12, ICOSAHEDRON_EDGES)


def platonic_solids():
    """The Platonic solids provided here, as (name, Graph) pairs, in order of edge count.

    The names are 'tetrahedron', 'octahedron', 'dodecahedron', 'icosahedron'; they are
    the keys used in data/platonic_solids.json.
    """
    return [("tetrahedron", tetrahedron()),
            ("octahedron", octahedron()),
            ("dodecahedron", dodecahedron()),
            ("icosahedron", icosahedron())]


def antiprism(n):
    """The n-antiprism: two n-cycles joined by a zigzag band; 4-regular, planar.

    antiprism(3) is the octahedron.  The antiprisms are the other infinite family of
    4-regular polyhedra, so they are the natural place to look for more octahedron-like
    behaviour.
    """
    if n < 3:
        raise ValueError("n must be at least 3")
    edges = []
    for i in range(n):
        edges.append((i, (i + 1) % n))                       # top cycle
        edges.append((n + i, n + (i + 1) % n))               # bottom cycle
        edges.append((i, n + i))                             # the band
        edges.append((i, n + (i - 1) % n))
    return Graph(2 * n, sorted(set(tuple(sorted(e)) for e in edges)))


def prism_graph(n):
    """The n-prism: two n-cycles joined by a perfect matching; 3-regular, planar.

    prism_graph(3) is the triangular prism.
    """
    if n < 3:
        raise ValueError("n must be at least 3")
    edges = []
    for i in range(n):
        edges.append((i, (i + 1) % n))
        edges.append((n + i, n + (i + 1) % n))
        edges.append((i, n + i))
    return Graph(2 * n, sorted(set(tuple(sorted(e)) for e in edges)))


def prism():
    """The triangular prism C_3 x K_2: planar, 6 vertices, 9 edges.

    A useful negative control: same vertex count as the octahedron, no H_2 torsion at
    k = 4.  Equal to prism_graph(3).
    """
    return Graph(6, [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5),
                     (0, 3), (1, 4), (2, 5)])


def theta(paths=3, length=2):
    """The theta graph: two hub vertices joined by `paths` internally disjoint paths.

    Each path has `length` edges, so `length` must be at least 2 for the result to be
    simple.  theta(3, 2) is one of An--Knudsen's three planar generators of H_2.
    """
    if length < 2:
        raise ValueError("length must be at least 2 (otherwise there are multi-edges)")
    u, w = 0, 1
    edges = []
    nxt = 2
    for _ in range(paths):
        prev = u
        for _ in range(length - 1):
            edges.append((prev, nxt))
            prev = nxt
            nxt += 1
        edges.append((prev, w))
    return Graph(nxt, edges)


def wheel(n):
    """The wheel W_n: a cycle of length n plus a hub joined to every rim vertex.

    Rim vertices are 0..n-1 (in cyclic order) and the hub is n; n >= 3.
    """
    edges = [(i, (i + 1) % n) for i in range(n)] + [(n, i) for i in range(n)]
    return Graph(n + 1, edges)


def petersen():
    """The Petersen graph, as the Kneser graph on 2-subsets of a 5-set.

    Vertices are the 2-subsets of {0,..,4} in lexicographic order, adjacent when disjoint:
    10 vertices, 15 edges, 3-regular, non-planar.
    """
    sets = list(itertools.combinations(range(5), 2))
    idx = {s: i for i, s in enumerate(sets)}
    edges = [(idx[s], idx[t]) for s, t in itertools.combinations(sets, 2)
             if not set(s) & set(t)]
    return Graph(10, edges)


# =========================================================================
# The Petersen family (needs networkx, for isomorphism testing)
# =========================================================================

def petersen_family():
    """The seven graphs of the Petersen family, as (name, Graph) pairs.

    By Robertson--Seymour--Thomas these are exactly the forbidden minors for linkless
    embeddability, obtained as the closure of K_6 under Delta-Y and Y-Delta exchanges.
    This function performs that closure rather than hard-coding the list, so the family
    is produced by the construction that characterises it.

    Returns the graphs sorted by (vertex count, degree sequence).  Requires networkx.
    """
    import networkx as nx

    def delta_Y(G):
        """Replace a triangle by a claw."""
        out = []
        for a, b, c in itertools.combinations(G.nodes(), 3):
            if G.has_edge(a, b) and G.has_edge(b, c) and G.has_edge(a, c):
                H = G.copy()
                H.remove_edges_from([(a, b), (b, c), (a, c)])
                x = max(H.nodes()) + 1
                H.add_edges_from([(x, a), (x, b), (x, c)])
                if min(dict(H.degree()).values()) >= 1:
                    out.append(nx.convert_node_labels_to_integers(H))
        return out

    def Y_delta(G):
        """Replace a degree-3 vertex by a triangle, skipping moves that create a
        multi-edge (the family is a family of simple graphs)."""
        out = []
        for x in list(G.nodes()):
            if G.degree(x) != 3:
                continue
            a, b, c = list(G.neighbors(x))
            if G.has_edge(a, b) or G.has_edge(b, c) or G.has_edge(a, c):
                continue
            H = G.copy()
            H.remove_node(x)
            H.add_edges_from([(a, b), (b, c), (a, c)])
            out.append(nx.convert_node_labels_to_integers(H))
        return out

    fam = [nx.complete_graph(6)]
    frontier = list(fam)
    while frontier:
        nxt = []
        for G in frontier:
            for H in delta_Y(G) + Y_delta(G):
                if not any(nx.is_isomorphic(H, F) for F in fam):
                    fam.append(H)
                    nxt.append(H)
        frontier = nxt
    fam.sort(key=lambda G: (G.number_of_nodes(),
                            sorted(d for _, d in G.degree())))

    def name_of(G):
        """A readable name for a member of the family: K_6, K_{3,3,1}, Petersen and
        K_{4,4}-e are recognised by isomorphism, the others are called G_{N,M}."""
        if nx.is_isomorphic(G, nx.complete_graph(6)):
            return "K_6"
        if nx.is_isomorphic(G, nx.complete_multipartite_graph(3, 3, 1)):
            return "K_{3,3,1}"
        if nx.is_isomorphic(G, nx.petersen_graph()):
            return "Petersen"
        K44 = nx.complete_bipartite_graph(4, 4)
        K44.remove_edge(0, 4)
        if nx.is_isomorphic(G, K44):
            return "K_{4,4}-e"
        return "G_%d,%d" % (G.number_of_nodes(), G.number_of_edges())

    out = []
    used = set()
    for G in fam:
        nm = name_of(G)
        while nm in used:
            nm += "'"
        used.add(nm)
        out.append((nm, from_networkx(G)))
    return out
