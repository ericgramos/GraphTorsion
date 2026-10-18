"""
Exact Smith normal form of an integer matrix.

Returns (rank, invariant factors).  The algorithm is the classical one:

  1. pick a pivot of least absolute value in the remaining block;
  2. Euclidean-reduce its row and its column, swapping the smaller remainder into the
     pivot position whenever one appears;
  3. if some entry of the remaining block is not divisible by the pivot, add that entry's
     row to the pivot row and go back to 2.

Step 3 is what makes the output the true Smith form (a divisibility chain) rather than
merely diagonal.  Everything is done in Python integers, so there is no overflow and no
rounding; the cost is that the matrices handed to this routine should already be small.
That is the job of the discrete Morse reduction in `chain`, `fast` and `lean`: by the
time a matrix reaches `smith` it is typically a few hundred rows across.

No pivoting heuristics beyond "least absolute value" are used.  On the residual matrices
that come out of the Morse reduction, entry growth has not been a problem in any
computation in this repository.
"""


def smith(M, nr, nc):
    """Smith normal form of the nr x nc integer matrix M (a list of row lists).

    Returns
    -------
    (rank, divs) : (int, list of int)
        `divs` is the list of diagonal entries in divisibility order, so the invariant
        factors are the entries that are not 0 or 1, and `rank` is len(divs).

    The argument is copied, not modified.
    """
    M = [row[:] for row in M]
    divs = []
    r0 = c0 = 0
    while r0 < nr and c0 < nc:
        # --- choose a pivot of least absolute value (short-circuit on a unit) ---
        piv = None
        best = None
        for i in range(r0, nr):
            Mi = M[i]
            for j in range(c0, nc):
                v = Mi[j]
                if v:
                    a = abs(v)
                    if best is None or a < best:
                        best, piv = a, (i, j)
                        if a == 1:
                            break
            if best == 1:
                break
        if piv is None:
            break                                  # the remaining block is zero
        i, j = piv
        M[r0], M[i] = M[i], M[r0]
        for row in M:
            row[c0], row[j] = row[j], row[c0]

        again = True
        while again:
            again = False
            # --- clear the pivot column ---
            p = M[r0][c0]
            for i in range(r0 + 1, nr):
                if M[i][c0]:
                    q = M[i][c0] // p
                    if q:
                        for j in range(c0, nc):
                            M[i][j] -= q * M[r0][j]
                    if M[i][c0]:                   # a remainder: swap it into the pivot
                        M[r0], M[i] = M[i], M[r0]
                        again = True
            # --- clear the pivot row ---
            p = M[r0][c0]
            for j in range(c0 + 1, nc):
                if M[r0][j]:
                    q = M[r0][j] // p
                    if q:
                        for i in range(r0, nr):
                            M[i][j] -= q * M[i][c0]
                    if M[r0][j]:
                        for i in range(r0, nr):
                            M[i][c0], M[i][j] = M[i][j], M[i][c0]
                        again = True
            # --- divisibility: the pivot must divide the rest of the block ---
            if not again:
                p = M[r0][c0]
                bad = None
                for i in range(r0 + 1, nr):
                    for j in range(c0 + 1, nc):
                        if M[i][j] % p:
                            bad = (i, j)
                            break
                    if bad:
                        break
                if bad:
                    i, _ = bad
                    for j in range(c0, nc):
                        M[r0][j] += M[i][j]
                    again = True

        divs.append(abs(M[r0][c0]))
        r0 += 1
        c0 += 1
    return r0, divs
