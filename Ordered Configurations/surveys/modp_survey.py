"""modp_survey.py -- torsion detection through field ranks (universal coefficients).

For each job: dim H_i over F_P for a large prime P = 1000003 (the rational Betti numbers b_i unless there
is P-torsion) and over F_p for small primes p.  By the universal coefficient theorem
     dim H_i(F_p) = b_i + t_i(p) + t_{i-1}(p),
t_i(p) = number of p-primary cyclic summands of H_i(Z); hence t_i(p) = sum_{j<=i} (-1)^(i-j) (dim H_j(F_p) - b_j).
Usage:  python3 modp_survey.py jobs.json out.jsonl [primes, default 2,3,5,7]
"""
import _common  # sets up the module path and ENGINE (see _common.py)
import json, sys, time, os
from gconf import homology

BIG = 1000003
jobs = json.load(open(sys.argv[1]))
outfn = sys.argv[2]
primes = [int(x) for x in (sys.argv[3] if len(sys.argv) > 3 else "2,3,5,7").split(",")]
done = set()
if os.path.exists(outfn):
    done = {(json.loads(l)["name"], json.loads(l)["n"]) for l in open(outfn)}
for job in jobs:
    if (job["name"], job["n"]) in done:
        continue
    t = time.time()
    rq = homology(job["edges"], job["n"], prime=BIG, tag="mp")
    b = rq["betti"]
    rec = dict(name=job["name"], edges=job["edges"], n=job["n"], betti=b, chi_ok=rq["chi"] == rq["gal"], torsion_counts={})
    for p in primes:
        rp = homology(job["edges"], job["n"], prime=p, tag="mp")["betti"]
        diff = [x - y for x, y in zip(rp, b)]
        tcount, prev = [], 0
        for d in diff:
            cur = d - prev
            tcount.append(cur)
            prev = cur
        rec["torsion_counts"][p] = tcount
    rec["secs"] = round(time.time() - t)
    with open(outfn, "a") as f:
        f.write(json.dumps(rec) + "\n")
    flag = {p: tc for p, tc in rec["torsion_counts"].items() if any(tc)}
    print(job["name"], job["n"], b, ("TORSION " + str(flag)) if flag else f"no p-torsion for p in {primes}", f"{rec['secs']}s", flush=True)
