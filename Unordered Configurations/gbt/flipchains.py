"""
Flip chains: explicit chains in the unreduced Swiatkowski complex of Gamma = K_{n+1,b}.

THE ONE-PARTICLE FLIP CHAIN
---------------------------
For an edge e = {a, w} with a in A and w in B, put

        z(e)  =  h_{w,e}  -  h_{a,e}          (a particle crossing e, from either end).

In the unreduced Swiatkowski complex d h_{v,e} = e - v, so the two "park on e" terms cancel:

        d z(e)  =  a - w                       (sitting particles).

PRODUCTS
--------
Chains supported on disjoint vertex sets multiply (the complex is a tensor product over the
vertices), and the differential is a derivation for this product.  For a bijection

        g : B_0 u {*}  -->  A,      B_0 = {w_1 < ... < w_n},   g(*) = the "hole",

the edges e_i = {g(w_i), w_i} are disjoint, so Z_g = z(e_1) z(e_2) ... z(e_n) is defined.  Then

        X_{B_0} = sum_g sgn(g) Z_g ,                       (sgn computed on (w_1,...,w_n,*))
        Y_T     = sum_{g : T -> A injective} sgn(g; x, y) Z_g . (x - y),

where in Y_T the set A \\ g(T) = {x < y} and sgn(g; x, y) is the sign of the sequence
(g(t_1), ..., g(t_{n-1}), x, y) as a permutation of A.  Then

        d X_{B_0}  =  (-1)^{n-1}  sum_{T subset B_0, |T| = n-1}  Y_T ,

which is the column B_0 of the unsigned inclusion matrix W_{n-1,n}(b).

ENCODING
--------
Chains are returned as {cell: coefficient} with cells in the encoding of `gbt.swiatkowski`
(unreduced): (mono, st) with st[v] = 1 for a sitting particle and 2 + edge index for a moving
one.  Vertices are A = 0..n and B = n+1..n+b, exactly as `complete_bip(n+1, b)` builds them.
"""

import itertools


def perm_sign(seq):
    """Sign of the permutation presented by a sequence of distinct comparable items."""
    s = 1
    a = list(seq)
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if a[i] > a[j]:
                s = -s
    return s


def edge_index(n, b, a, w):
    """Index of the edge {a, w} in complete_bip(n+1, b): lexicographic in (a, w)."""
    return a * b + (w - (n + 1))


def word_to_cell(N, word):
    """Turn an ordered product of generators into (cell, sign).

    `word` is a list of (vertex, state) with state = 1 (sitting, degree 0) or 2 + e (moving
    along edge e, degree 1).  Sorting into increasing vertex order costs the sign of the
    permutation of the DEGREE-ONE entries only, since degree-zero entries commute freely.
    """
    st = [0] * N
    for v, s in word:
        assert st[v] == 0, "two generators at one vertex: not a product of disjoint chains"
        st[v] = s
    movers = [v for v, s in word if s >= 2]
    return ((), tuple(st)), perm_sign(movers)


def expand(N, factors):
    """Multiply out a list of factors, each a list of (coefficient, generator).

    Returns {cell: coefficient}.  The factors are multiplied in the order given.
    """
    out = {}
    for choice in itertools.product(*factors):
        coef = 1
        word = []
        for c, gen in choice:
            coef *= c
            word.append(gen)
        cell, sgn = word_to_cell(N, word)
        out[cell] = out.get(cell, 0) + coef * sgn
        if out[cell] == 0:
            del out[cell]
    return out


def z(n, b, a, w):
    """The flip chain z({a, w}) as a factor: h_{w,e} - h_{a,e}."""
    e = edge_index(n, b, a, w)
    return [(+1, (w, 2 + e)), (-1, (a, 2 + e))]


def add_into(acc, chain, c=1):
    for k, v in chain.items():
        acc[k] = acc.get(k, 0) + c * v
        if acc[k] == 0:
            del acc[k]
    return acc


def X_B0(n, b, B0):
    """X_{B_0} = sum over bijections g : B_0 u {*} -> A of sgn(g) z(e_1)...z(e_n)."""
    N = n + 1 + b
    B0 = sorted(B0)
    acc = {}
    for img in itertools.permutations(range(n + 1)):   # img = (g(w_1),...,g(w_n), g(*))
        factors = [z(n, b, img[i], B0[i]) for i in range(n)]
        add_into(acc, expand(N, factors), perm_sign(img))
    return acc


def Y_T(n, b, T):
    """Y_T = sum over injections g : T -> A of sgn(g; x, y) z(...)...z(...) (x - y)."""
    N = n + 1 + b
    T = sorted(T)
    acc = {}
    for img in itertools.permutations(range(n + 1), n - 1):
        x, y = sorted(set(range(n + 1)) - set(img))
        factors = [z(n, b, img[i], T[i]) for i in range(n - 1)]
        factors.append([(+1, (x, 1)), (-1, (y, 1))])
        add_into(acc, expand(N, factors), perm_sign(list(img) + [x, y]))
    return acc


def X(n, b):
    """X = sum over all n-subsets B_0 of X_{B_0}, with dX = +-(b-n+1) y for y primitive."""
    acc = {}
    for B0 in itertools.combinations(range(n + 1, n + 1 + b), n):
        add_into(acc, X_B0(n, b, B0))
    return acc


def closed_formula_X(n, b):
    """The coefficients of X from its closed formula, without expanding the products.

    For a matching configuration (S, t) with S_A = S n A, S_B = S n B, the coefficient is
    (-1)^{|S_A|} sgn sigma(S,t), where sigma lists A as
        [ S_A increasing | t(v) for v in S_B, in increasing order of v | the hole ].
    Returns {cell: coefficient} over all matching configurations, in the same encoding as X().
    """
    N = n + 1 + b
    A = list(range(n + 1))
    B = list(range(n + 1, N))
    out = {}
    for A0 in itertools.combinations(A, n):
        hole = [a for a in A if a not in A0][0]
        for B0 in itertools.combinations(B, n):
            for perm in itertools.permutations(range(n)):
                pairs = [(A0[i], B0[perm[i]]) for i in range(n)]
                for marks in itertools.product((0, 1), repeat=n):
                    t = {}
                    for (a, w), mk in zip(pairs, marks):
                        if mk == 0:
                            t[a] = w
                        else:
                            t[w] = a
                    SA = sorted(v for v in t if v <= n)
                    SB = sorted(v for v in t if v > n)
                    sigma = SA + [t[v] for v in SB] + [hole]
                    st = [0] * N
                    for v, u in t.items():
                        st[v] = 2 + edge_index(n, b, min(v, u), max(v, u))
                    out[((), tuple(st))] = (-1) ** len(SA) * perm_sign(sigma)
    return out
