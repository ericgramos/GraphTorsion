"""local_survey.py -- exact p-primary torsion (local Smith form, engine -q) for a list of jobs, plus the
rational Betti numbers from a large prime.  Usage:  python3 local_survey.py jobs.json out.jsonl [primes]
Result: data/survey/triangle_local_snf.jsonl, local.jsonl.
"""
import _common  # sets up the module path and ENGINE (see _common.py)
import json, sys, time, os
from gconf import homology, local_torsion
jobs = json.load(open(sys.argv[1])); outfn = sys.argv[2]
primes = [int(x) for x in (sys.argv[3] if len(sys.argv) > 3 else "2,3").split(",")]
done = set()
if os.path.exists(outfn):
    done = {(json.loads(l)["name"], json.loads(l)["n"]) for l in open(outfn)}
for job in jobs:
    if (job["name"], job["n"]) in done: continue
    t = time.time()
    rq = homology(job["edges"], job["n"], prime=1000003, tag="lsq")
    rec = dict(name=job["name"], edges=job["edges"], n=job["n"], betti=rq["betti"], chi_ok=rq["chi"] == rq["gal"], torsion={})
    for p in primes:
        r = local_torsion(job["edges"], job["n"], p, tag="lsl")
        rec["torsion"][p] = r["torsion"]
        # summands of order >= p^k are invisible mod p^k; they are absent exactly when the number of
        # elementary divisors of order < p^k (rank_local) equals the rational rank of every map
        rec.setdefault("orders_complete", {})[p] = r["rank_local"] == rq["rank"]
        if not rec["orders_complete"][p]:
            print(f"  WARNING: p={p}: some elementary divisors are divisible by p^k (ranks {r['rank_local']} vs {rq['rank']})")
    rec["secs"] = round(time.time() - t)
    with open(outfn, "a") as f: f.write(json.dumps(rec) + "\n")
    tors = {p: {i: t for i, t in enumerate(tt) if t} for p, tt in rec["torsion"].items()}
    tors = {p: v for p, v in tors.items() if v}
    print(job["name"], job["n"], rec["betti"], ("TORSION " + str(tors)) if tors else f"no p-torsion, p in {primes}", f"{rec['secs']}s", flush=True)
