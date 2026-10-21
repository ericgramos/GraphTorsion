"""solve_residual.py -- exact integral solution of the small residual system left by exactsolve.

exactsolve writes a file
    nr nc                 size of the residual block R (rows = residual rows, columns = residual columns)
    i j value             the nonzero entries of R (one per line)
    X
    j value               the nonzero entries of the target vector x
and asks for an integer row vector lam with  lam * R = x.  This program writes the nonzero entries of lam
("i value" per line) to the output file.

Method.
  1. Choose columns J forming a basis of the column space of R: the pivot columns of the row-reduced form of
     R mod 1000003, checked to give the same rank mod 998244353.  Since x lies in the rational row space of R,
     every column of R outside J is a rational combination of the columns in J, and x satisfies the same
     relations; so any solution of lam * R_J = x_J solves the whole system.  (This is rechecked at the end
     over Z on all columns, so a wrong choice of J could only make the program fail, never give a wrong lam.)
  2. Hermite normal form of [R_J | I] (python-flint) gives T with T R_J = H, H in row echelon form.
  3. Solve mu * H = x_J by walking the pivots of H from left to right (each pivot column of H is touched
     only by its own row and earlier ones); a non-integral quotient means there is no integral solution.
  4. lam = mu * T.
Usage:  python3 solve_residual.py residual_file lam_file
"""
import sys
import flint


def read(fn):
    """Read the residual file: its size, the entries of R and the target vector x."""
    with open(fn) as f:
        nr, nc = map(int, f.readline().split())
        ent = []
        for line in f:
            if line.startswith("X"):
                break
            a, b, c = line.split()
            ent.append((int(a), int(b), int(c)))
        x = [0] * nc
        for line in f:
            if line.strip():
                a, b = line.split()
                x[int(a)] = int(b)
    return nr, nc, ent, x


def column_basis(nr, nc, ent, p):
    """Pivot columns of the reduced row echelon form of R mod p."""
    M = flint.nmod_mat(nr, nc, p)
    for a, b, c in ent:
        M[a, b] = c % p
    E = M.rref()[0]
    piv, j = [], 0
    for r in range(nr):
        while j < nc and int(E[r, j]) == 0:
            j += 1
        if j == nc:
            break
        piv.append(j)
        j += 1
    return piv


def main(fn, out):
    """Solve lam * R = x over Z and write lam (see the module docstring)."""
    nr, nc, ent, x = read(fn)
    if nr == 0:
        assert not any(x), "empty residual but nonzero target"
        open(out, "w").close()
        return
    J = column_basis(nr, nc, ent, 1000003)
    assert len(J) == len(column_basis(nr, nc, ent, 998244353)), "rank differs between the two primes"
    print(f"rank {len(J)} of {nr} x {nc}", file=sys.stderr, flush=True)
    pos = {j: k for k, j in enumerate(J)}
    r_ = len(J)
    A = flint.fmpz_mat(nr, r_ + nr)
    for a, b, c in ent:
        if b in pos:
            A[a, pos[b]] = c
    for r in range(nr):
        A[r, r_ + r] = 1
    Hf = A.hnf()
    del A
    # mu * H = x_J
    rem = [x[j] for j in J]
    mu = {}
    for r in range(nr):
        piv = next((j for j in range(r_) if Hf[r, j] != 0), None)
        if piv is None:
            break
        q, m = divmod(rem[piv], int(Hf[r, piv]))
        if m:
            sys.exit("no integral solution")
        if q:
            mu[r] = q
            for j in range(piv, r_):
                h = Hf[r, j]
                if h != 0:
                    rem[j] -= q * int(h)
    if any(rem):
        sys.exit("target is not in the row span")
    lam = [0] * nr
    for k, q in mu.items():
        for r in range(nr):
            t = Hf[k, r_ + r]
            if t != 0:
                lam[r] += q * int(t)
    # recheck on ALL columns, over Z
    chk = [0] * nc
    for a, b, c in ent:
        if lam[a]:
            chk[b] += lam[a] * c
    assert chk == x, "solution fails on the remaining columns"
    with open(out, "w") as f:
        for r, v in enumerate(lam):
            if v:
                f.write(f"{r} {v}\n")
    print(f"residual solved: {sum(1 for v in lam if v)} nonzero multipliers, max |lam| {max(abs(v) for v in lam)}",
          file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
