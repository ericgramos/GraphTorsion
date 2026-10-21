"""check_generator_tools.py -- end-to-end check of the generator, subgraph and torus tools (October 2026).

Reruns, from the files in the package, every computation behind the explicit generators of the class of order 16
(data/survey/k4_family/support) and compares the results with the saved ones.  Each step prints OK or FAIL.

  certs      verify every certificate in certificates/ with validation/verify_cert.py (about 10 minutes)
  gens       certq generator mode on the H-first edge order; the image test with orderq must reproduce
             image_triangle_leaf_Hfirst_K10.json exactly
  exact      exactsolve on x = z_2 + z_11 must write it as y + du with y on H, and the certificate built from that y
             must verify (y depends on pivot choices, so it is checked rather than compared); reducecyc must reproduce
             the saved homology coordinates of the saved y; exactsolve with h = 0 must write 16 y* as a boundary, and
             the certificate built from it must verify
  tori       orderq with the tori of H as relations must find y in their span (as saved)
  span       spanq against orderq: for random sets of tori, the torsion in their span modulo boundaries computed by the
             fast path (reducemod + spanq) and by the independent slow path (orderq with the tori as relations) must
             agree; plus the exact test (exact_check.py) on the two-tori set must give an odd cocycle value
Usage:  python3 validation/check_generator_tools.py [step ...]       (default: all steps; about 2 hours, 6 GB)
Run from the top folder of the package after `make`.  Work files go to a temporary folder.
"""
import gzip
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile

TOP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
C = os.path.join(TOP, "certificates")
SUP = os.path.join(TOP, "data", "survey", "k4_family", "support")
G = os.path.join(SUP, "graphs")
sys.path.insert(0, SUP)
sys.path.insert(0, os.path.join(SUP, "tori"))
from kernel import analyse, structure  # noqa: E402

WORK = tempfile.mkdtemp(prefix="gentools_")
results = []


def run(cmd, stdin=None, env=None, out=None):
    """Run a command with stdin from a file; optionally save its stdout."""
    e = dict(os.environ)
    if env:
        e.update(env)
    with open(stdin) if stdin else open(os.devnull) as fi:
        r = subprocess.run(cmd, stdin=fi, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=e)
    if out:
        open(out, "w").write(r.stdout)
    return r


def report(name, ok, detail=""):
    """Record and print the outcome of one check."""
    results.append((name, ok))
    print(f"{'OK  ' if ok else 'FAIL'} {name} {detail}", flush=True)


def gunzip(src, dst):
    """Decompress a .gz file."""
    with gzip.open(src, "rt") as f, open(dst, "w") as o:
        shutil.copyfileobj(f, o)


def torsion_of(json_text):
    """The structure of T cap (span) from an orderq/spanq JSON output (kernel.py)."""
    fn = os.path.join(WORK, "tmp_oq.json")
    open(fn, "w").write(json_text)
    d, oc, ks, tot = analyse(fn)
    return structure(oc)


def step_certs():
    """Verify every certificate of the package."""
    for f in sorted(os.listdir(C)):
        if f.startswith("cert_") and f.endswith(".txt.gz"):
            r = run(["python3", os.path.join(TOP, "validation", "verify_cert.py"), os.path.join(C, f)])
            report(f"verify {f}", r.returncode == 0 and "VERIFIED" in r.stdout)


def gens_file():
    """Generators of the torsion in the H-first edge order (certq generator mode), made once per run."""
    g = os.path.join(WORK, "gens_Hfirst.txt")
    if not os.path.exists(g):
        run([os.path.join(C, "certq"), "2", "2", "4", "-1", g], stdin=os.path.join(G, "K4dbl5_Hfirst_n5.txt"))
    return g


def step_gens():
    """The orderq image test for the triangle with a leaf must reproduce the saved output."""
    g = gens_file()
    out = os.path.join(WORK, "img.json")
    run([os.path.join(C, "orderq"), "2", "2", "10", g, "6,7,8,9,10", "-"], stdin=os.path.join(G, "K4dbl5_Hfirst_n5.txt"), out=out)
    saved = open(os.path.join(SUP, "image_triangle_leaf_Hfirst_K10.json")).read()
    report("orderq image of the triangle with a leaf", open(out).read() == saved)


def write_cert(path, W, gens, subgraph):
    """Write a certificate of order 16: the chain W (dict cell -> coef), the cocycle F of the order-16 generator."""
    nums = open(os.path.join(G, "K4dbl5_Hfirst_n5.txt")).read().split()
    with open(path, "w") as o:
        o.write(f"GRAPH {nums[0]} {nums[1]} {nums[2]}\n")
        for i in range(int(nums[1])):
            o.write(f"EDGE {nums[3 + 2 * i]} {nums[4 + 2 * i]}\n")
        o.write(f"DEGREE 2\nPRIME 16\nSUBGRAPH {subgraph}\n")
        for cell, a in W.items():
            if a:
                o.write(f"W {a} {cell}\n")
        cur = None
        for line in open(gens):
            if line.startswith("Z "):
                cur = int(line.split()[1])
            elif cur == 11 and line.startswith("F "):
                o.write(line)


def read_chain(path, prefix=None):
    """A chain file ("coef cell" lines, or lines "TAG coef cell" when prefix is given) as a dict."""
    out = {}
    for line in open(path):
        if prefix:
            if not line.startswith(prefix + " "):
                continue
            line = line[len(prefix) + 1:]
        a, cell = line.rstrip("\n").split(" ", 1)
        out[cell] = out.get(cell, 0) + int(a)
    return out


def step_exact():
    """exactsolve must write x = z_2 + z_11 as y + d u with y on H, and the certificate built from that y must verify
    (the solution is not unique: it depends on the pivot choices, so it is checked, not compared); reducecyc must
    reproduce the saved homology coordinates of the saved y; exactsolve with h = 0 must write 16 y* as a boundary,
    and the two-tori certificate built from it must verify."""
    g = gens_file()
    x = os.path.join(WORK, "x.txt")
    gunzip(os.path.join(SUP, "target_x_z2_plus_z11.txt.gz"), x)
    pre = os.path.join(WORK, "sol")
    env = {"SOLVER": os.path.join(C, "solve_residual.py")}
    r = run([os.path.join(C, "exactsolve"), "2", "6", x, pre], stdin=os.path.join(G, "K4dbl5_Hfirst_n5.txt"), env=env)
    report("exactsolve writes x = y + du with y on the triangle with a leaf", r.returncode == 0)
    if r.returncode == 0:
        # W = 8 w_2 + w_11 - 16 u, so that dW = 16 (z_2 + z_11) - 16 du = 16 y
        W, cur = {}, None
        for line in open(g):
            if line.startswith("Z "):
                cur = int(line.split()[1])
            elif line.startswith("W ") and cur in (2, 11):
                a, cell = line[2:].rstrip("\n").split(" ", 1)
                W[cell] = W.get(cell, 0) + (8 if cur == 2 else 1) * int(a)
        for cell, a in read_chain(pre + ".u").items():
            W[cell] = W.get(cell, 0) - 16 * a
        cert = os.path.join(WORK, "cert_triangle.txt")
        write_cert(cert, W, g, 6)
        r = run(["python3", os.path.join(TOP, "validation", "verify_cert.py"), cert])
        report("the certificate built from this y verifies", r.returncode == 0 and "VERIFIED" in r.stdout)
    y_saved = os.path.join(WORK, "y_saved.txt")
    gunzip(os.path.join(SUP, "triangle_leaf_cycle_y.txt.gz"), y_saved)
    co = os.path.join(WORK, "coords.txt")
    run([os.path.join(C, "reducecyc"), "2", y_saved, co], stdin=os.path.join(G, "triangle_leaf_n5.txt"))
    c_saved = os.path.join(WORK, "coords_saved.txt")
    gunzip(os.path.join(SUP, "triangle_leaf_cycle_y_homology_coords.txt.gz"), c_saved)
    report("reducecyc reproduces the homology coordinates of the saved y", open(co).read() == open(c_saved).read())
    ys = os.path.join(SUP, "tori", "two_tori_cycle_ystar.txt")
    y16 = os.path.join(WORK, "ystar16.txt")
    with open(y16, "w") as o:
        for cell, a in read_chain(ys).items():
            o.write(f"{16 * a} {cell}\n")
    pw = os.path.join(WORK, "solW")
    r = run([os.path.join(C, "exactsolve"), "2", "0", y16, pw], stdin=os.path.join(G, "K4dbl5_Hfirst_n5.txt"), env=env)
    ok = r.returncode == 0 and '"y_cells":0' in r.stdout
    report("exactsolve writes 16 y* as a boundary", ok)
    if ok:
        cert = os.path.join(WORK, "cert_two_tori.txt")
        write_cert(cert, read_chain(pw + ".u"), g, 6)
        r = run(["python3", os.path.join(TOP, "validation", "verify_cert.py"), cert])
        report("the rebuilt two-tori certificate verifies", r.returncode == 0 and "VERIFIED" in r.stdout)


def step_tori():
    """y must lie in the span of the tori of H."""
    tori = os.path.join(WORK, "tori.txt")
    gunzip(os.path.join(SUP, "tori_triangle_leaf_n5.txt.gz"), tori)
    yb = os.path.join(WORK, "yb.txt")
    gunzip(os.path.join(SUP, "tori", "y_block_six_edges.txt.gz"), yb)
    out = os.path.join(WORK, "ytori.json")
    run([os.path.join(C, "orderq"), "2", "2", "8", yb, "-", "-", tori], stdin=os.path.join(G, "triangle_leaf_n5.txt"), out=out)
    d = json.load(open(out))["cycles"][0]
    report("y lies in the span of the tori of H (mod 2^8, no free part)", d["order_in_M"] == 1 and d["free_nonzero"] == 0)


def step_span():
    """The fast span path (reducemod + spanq) must agree with orderq on several sets of tori; the exact check of the two-tori set."""
    g = gens_file()
    tori = os.path.join(WORK, "tori.txt")
    if not os.path.exists(tori):
        gunzip(os.path.join(SUP, "tori_triangle_leaf_n5.txt.gz"), tori)
    # blocks: generators then tori, as minimal_tori.py expects
    blocks = os.path.join(WORK, "blocks.txt")
    with open(blocks, "w") as o:
        for line in open(g):
            if line.startswith("Z "):
                _, i, od = line.split()
                o.write(f"Z gen {i} {od}\n")
            elif line[:1].isdigit() or line[:1] == "-":
                o.write(line)
        for line in open(tori):
            o.write(("Z tor " + line[2:]) if line.startswith("Z ") else line)
    red = os.path.join(WORK, "blocks.red")
    run([os.path.join(C, "reducemod"), "2", "2", "24", blocks, red], stdin=os.path.join(G, "K4dbl5_Hfirst_n5.txt"))
    from minimal_tori import Oracle
    orc = Oracle(os.path.join(C, "spanq"), red)
    ntor = sum(1 for line in open(tori) if line.startswith("Z "))
    random.seed(20261020)
    two = [18125, 18132, 18171, 18175, 18176, 18177, 19894, 19903, 19906, 19908, 19909, 19911, 19949, 19950]
    tests = [("the 14 tori of the two-tori class", two)] + \
            [(f"{k} random tori", sorted(random.sample(range(ntor), k))) for k in (10, 300, 3000)]
    for name, ids in tests:
        fast = torsion_of(json.dumps(orc.span(ids)))
        rel = os.path.join(WORK, "rel.txt")
        want = set(ids)
        keep = False
        with open(rel, "w") as o:
            for line in open(tori):
                if line.startswith("Z "):
                    keep = int(line.split()[1]) in want
                if keep:
                    o.write(line)
        out = os.path.join(WORK, "slow.json")
        run([os.path.join(C, "orderq"), "2", "2", "24", g, "-", "-", rel], stdin=os.path.join(G, "K4dbl5_Hfirst_n5.txt"), out=out)
        slow = torsion_of(open(out).read())
        report(f"spanq = orderq on {name}", fast == slow, f"(torsion from the span: {slow})")
    r = run(["python3", os.path.join(SUP, "tori", "exact_check.py"), tori, g, os.path.join(G, "K4dbl5_Hfirst_n5.txt")] + list(map(str, two)))
    res = json.loads(r.stdout.strip().splitlines()[-1]) if r.stdout.strip() else {}
    report("exact check of the two-tori set (one relation, odd cocycle value)", res.get("kernel_dim") == 1 and res.get("odd_vectors") == 1)


STEPS = {"certs": step_certs, "gens": step_gens, "exact": step_exact, "tori": step_tori, "span": step_span}

if __name__ == "__main__":
    for s in (sys.argv[1:] or list(STEPS)):
        STEPS[s]()
    print(f"\n{sum(ok for _, ok in results)} of {len(results)} checks passed; work files in {WORK}")
    sys.exit(0 if all(ok for _, ok in results) else 1)
