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
Every chain here has no parked particles, so mono is always ().

SIGN CONVENTIONS
----------------
A product of generators written in some order is converted to a basis element by sorting
the generators into increasing vertex order; only the degree-1 (moving) generators
anticommute, so the sign is the sign of the permutation that sorts the moving vertices
(`word_to_cell`).  `perm_sign` computes signs by counting inversions, so a "permutation of
A" is any sequence of distinct elements of A, compared with its sorted order.

EXAMPLE
-------
    >>> from gbt import flipchains, swiatkowski, complete_bip
    >>> X = flipchains.X(2, 3)                    # the chain for K_{3,3}, 2 particles
    >>> X == flipchains.closed_formula_X(2, 3)
    True
    >>> C = swiatkowski(complete_bip(3, 3), 2)
    >>> dX = {}
    >>> for c, a in X.items():
    ...     for r, v in C.B[2][c].items():
    ...         dX[r] = dX.get(r, 0) + a * v
    >>> sorted({abs(v) for v in dX.values() if v})       # dX = +-(b - n + 1) y, b - n + 1 = 2
    [2]

`experiments/torsion_class_X.py` and `experiments/wilson_block.py` check these identities
systematically.
"""

import itertools


def perm_sign(seq):
    """Sign of the permutation presented by a sequence of distinct comparable items.

    That is, (-1)^(number of inversions): +1 if sorting the sequence takes an even number
    of transpositions.  Quadratic in the length, which is at most n + 1 here.
    """
    s = 1
    a = list(seq)
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if a[i] > a[j]:
                s = -s
    return s


def edge_index(n, b, a, w):
    """Index of the edge {a, w} in complete_bip(n+1, b): lexicographic in (a, w).

    a in A = {0..n}, w in B = {n+1..n+b}; the index is a*b + (w - (n+1)).
    """
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
    """acc += c * chain, in place, for chains as dicts {cell: coefficient} (zeros dropped).

    Returns acc so that calls can be nested.
    """
    for k, v in chain.items():
        acc[k] = acc.get(k, 0) + c * v
        if acc[k] == 0:
            del acc[k]
    return acc


def X_B0(n, b, B0):
    """X_{B_0} = sum over bijections g : B_0 u {*} -> A of sgn(g) z(e_1)...z(e_n).

    n is the number of particles (|A| = n + 1), b = |B|, and B0 an n-subset of
    B = {n+1, ..., n+b}.  A bijection g is enumerated as the tuple
    (g(w_1), ..., g(w_n), g(*)) of all of A, and sgn(g) is the sign of that tuple; e_i is
    the edge {g(w_i), w_i}.  Returns the chain as {cell: coefficient}; it has
    (n+1)! * 2^n terms before cancellation.
    """
    N = n + 1 + b
    B0 = sorted(B0)
    acc = {}
    for img in itertools.permutations(range(n + 1)):   # img = (g(w_1),...,g(w_n), g(*))
        factors = [z(n, b, img[i], B0[i]) for i in range(n)]
        add_into(acc, expand(N, factors), perm_sign(img))
    return acc


def Y_T(n, b, T):
    """Y_T = sum over injections g : T -> A of sgn(g; x, y) z(...)...z(...) (x - y).

    T is an (n-1)-subset of B; for each injection g, {x < y} = A minus g(T), the factor
    (x - y) is "sit at x" minus "sit at y" (degree 0), and sgn(g; x, y) is the sign of the
    tuple (g(t_1), ..., g(t_{n-1}), x, y).  Returns {cell: coefficient}.
    """
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
    """X = sum over all n-subsets B_0 of X_{B_0}, with dX = +-(b-n+1) y for y primitive.

    Primitive means the coefficients of y have gcd 1 (in fact they are all +-1), so when
    b <= 2n - 1 (where H_n(B_n K_{n+1,b}) = 0 and so nothing smaller kills y) the class [y]
    has order exactly b - n + 1 in H_{n-1}.  Returns {cell: coefficient}; the number of
    terms grows like C(b, n) (n+1)! 2^n, so keep n <= 4.
    """
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

    A matching configuration is a top-degree cell in which the n moving particles use n
    pairwise disjoint edges.  It is enumerated as: the n A-vertices used (A0, so the hole
    is the remaining one), the n B-vertices used (B0), a bijection pairing them, and for
    each pair which endpoint is the source (marks: 0 = the A end moves, 1 = the B end
    moves).  t maps each source to the other end of its edge.
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
                    SA = sorted(v for v in t if v <= n)       # sources in A
                    SB = sorted(v for v in t if v > n)        # sources in B
                    # sigma lists all of A once: the A-sources, then the A-targets of the
                    # B-sources in order of their B-source (NOT in order of their own
                    # value -- the other choice does not give a chain with dX ~ y), then
                    # the hole
                    sigma = SA + [t[v] for v in SB] + [hole]
                    st = [0] * N
                    for v, u in t.items():
                        st[v] = 2 + edge_index(n, b, min(v, u), max(v, u))
                    out[((), tuple(st))] = (-1) ** len(SA) * perm_sign(sigma)
    return out
