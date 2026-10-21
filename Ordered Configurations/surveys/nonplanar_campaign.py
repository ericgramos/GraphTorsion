"""nonplanar_campaign.py -- every run of the non-planar torsion search, with its recorded result.

All runs of the campaign (October 12-13, 2026) are listed in data/nonplanar/campaign.json.  Each record has
    name        short identifier (also the name of its raw output files in data/nonplanar/raw/)
    V, n, edges the graph (edge list, oriented u -> v) and the number of particles
    options     the engine options used: coefficients (-q p K or -W), label group (-g), automorphisms
                (-a), character (-X, -Y), degree cap (-d), dense-phase settings (-D)
    cells, dims, torsion, modulus, group, rc, secs, peak_MB     the recorded outcome
dims[i] is the dimension of H_i of the (twisted) coinvariant complex over the residue field; torsion maps
"H_i" to the local Smith stage counts [number of Z/p summands, ...]; an empty dict means no torsion.

Usage:
    python3 nonplanar_campaign.py --list                 one line per run (name, graph size, result)
    python3 nonplanar_campaign.py --run NAME [NAME ...]  re-run with the current engine and compare
    python3 nonplanar_campaign.py --run-small            re-run every completed run with peak < 2 GB
    python3 nonplanar_campaign.py --rational NAME [...]  rank check: re-run over the large prime 1000003 with the
                                                         same group and character, and compare each map's rank
                                                         with the sum of its recorded local stages.  A local run
                                                         mod p^K cannot see summands of order >= p^K; they are
                                                         absent exactly when the two numbers agree for every map.
                                                         Agreeing ranks are saved in data/nonplanar/rational_ranks.json.
    python3 nonplanar_campaign.py --prime P NAME [...]   P-torsion test at another prime P (e.g. 7): ranks over F_P
                                                         compared with the saved rational ranks.

The runs were made with engine versions 17-26 (history/), which agree with engine/oconf on the regression
suite; '-F' (fast enumeration) is added automatically when a group is present, which does not change the
complex.  Large runs need up to 6 GB of memory and up to two hours.
"""
import json, os, sys, time
import _common

DATA = os.path.join(_common.HERE, "..", "data", "nonplanar", "campaign.json")
recs = json.load(open(DATA))


def show(r):
    """Print one line summarising a recorded campaign run."""
    res = r.get("result") or ("TORSION " + json.dumps(r["torsion"]) if r.get("torsion") else "no torsion")
    print(f"{r['name']:18s} V={r['V']:2d} E={len(r['edges']):2d} n={r['n']}  {' '.join(r['options'])[:70]:70s}  {res}")


def rerun(r):
    """Re-run a recorded campaign run with the current engine and compare the answer with the recorded one."""
    opts = list(r["options"])
    if ("-g" in opts or "-a" in opts) and "-F" not in opts:
        opts.append("-F")
    data, secs, peak, rc = _common.run_engine(opts, [tuple(e) for e in r["edges"]], r["n"], tag=r["name"])
    if rc != 0:
        print(f"{r['name']}: engine failed with code {rc}")
        return False
    c = data["cells"]
    rk = [0] * (len(c) + 1)
    tors = {}
    for m in data["maps"]:
        st = m.get("stages")
        rk[m["k"]] = st[0] if st else m["pivots"]
        if st and (len(st) > 1 or m.get("leftover_nnz")):
            tors[f"H{m['k'] - 1}"] = st[1:]
    dims = [c[i] - rk[i] - rk[i + 1] for i in range(len(data["maps"]))]
    ok = dims == r.get("dims") and tors == r.get("torsion")
    print(f"{r['name']}: dims {dims} torsion {tors}  recorded {r.get('dims')} {r.get('torsion')}  "
          f"{'AGREE' if ok else 'DIFFER'}  ({secs:.0f}s, {peak} MB)")
    return ok


RATIONAL = os.path.join(_common.HERE, "..", "data", "nonplanar", "rational_ranks.json")


def field_ranks(r, P, tag, local=False):
    """Re-run campaign record r over F_P (same group, character and degree cap) and return {k: rank of d_k}.
    local=True uses the local engine with K = 1 (-q P 1), which has the packed storage and the dense
    finishing phase and needs much less memory than the plain field routines on large complexes."""
    o = list(r["options"])
    i = o.index("-q")
    if local:
        o[i:i + 3] = ["-q", str(P), "1"]
    else:
        o[i:i + 3] = ["-p", str(P)]
        for flag, k in (("-D", 1), ("-M", 2), ("-C", 1)):     # dense-phase options only apply to local mode
            if flag in o:
                j = o.index(flag); del o[j:j + 1 + k]
    if ("-g" in o or "-a" in o) and "-F" not in o:
        o.append("-F")
    data, secs, peak, rc = _common.run_engine(o, [tuple(e) for e in r["edges"]], r["n"], tag=r["name"] + tag)
    if rc != 0:
        return None, secs, peak, rc
    from gconf import snf_divisors_modp
    ranks = {}
    if local:
        return {m["k"]: sum(m["stages"]) for m in data["maps"]}, secs, peak, rc
    for m in data["maps"]:                  # pivots found by the engine + rank of the residual block mod P
        extra = snf_divisors_modp(m["res_file"], P)[0] if m.get("res_rows") else 0
        ranks[m["k"]] = m["pivots"] + extra
        if os.path.exists(m["res_file"]):
            os.remove(m["res_file"])
    return ranks, secs, peak, rc


def rational(r):
    """Compare a recorded local run (mod p^K) with ranks over a large prime, to rule out summands of order >= p^K."""
    raw = os.path.join(_common.HERE, "..", "data", "nonplanar", "raw", r["name"] + ".json")
    if not os.path.exists(raw) or "-q" not in r["options"]:
        print(f"{r['name']}: no recorded local run to compare with")
        return None
    rec = json.load(open(raw))
    local = {m["k"]: sum(m["stages"]) for m in rec["maps"]}
    rat, secs, peak, rc = field_ranks(r, 1000003, "_rat")
    if rat is None:
        print(f"{r['name']}: engine failed with code {rc}")
        return None
    ok = all(rat[k] == local.get(k) for k in rat)
    print(f"{r['name']}: ranks over F_1000003 {[rat[k] for k in sorted(rat)]}, local stage sums "
          f"{[local.get(k) for k in sorted(rat)]}  {'AGREE: no summands of order >= p^K' if ok else 'DIFFER'}  ({secs:.0f}s, {peak} MB)")
    if ok:                                  # record the rational ranks for later --prime checks
        db = json.load(open(RATIONAL)) if os.path.exists(RATIONAL) else {}
        db[r["name"]] = {str(k): v for k, v in rat.items()}
        json.dump(db, open(RATIONAL, "w"), indent=1)
    return ok


def other_prime(r, P):
    """p-torsion test at another prime P: compare ranks over F_P with the recorded rational ranks.
    The number of cyclic summands of P-power order in H_{k-1} of the quotient complex equals
    rank_Q(d_k) - rank_P(d_k).  Requires P not dividing the group order (the engine refuses otherwise)."""
    db = json.load(open(RATIONAL)) if os.path.exists(RATIONAL) else {}
    if r["name"] not in db:
        print(f"{r['name']}: no recorded rational ranks; run --rational first")
        return None
    rat = {int(k): v for k, v in db[r["name"]].items()}
    rp, secs, peak, rc = field_ranks(r, P, f"_p{P}", local=P < 32)
    if rp is None:
        print(f"{r['name']}: engine failed with code {rc}" + (" (P divides the group order)" if rc == 10 else ""))
        return None
    tors = {f"H{k - 1}": rat[k] - rp[k] for k in sorted(rp) if rat.get(k) != rp[k]}
    print(f"{r['name']}: ranks over F_{P} {[rp[k] for k in sorted(rp)]}  "
          f"{'TORSION ' + json.dumps(tors) if tors else f'no {P}-torsion'}  ({secs:.0f}s, {peak} MB)", flush=True)
    return tors


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] == "--list":
        for r in recs:
            show(r)
    elif sys.argv[1] == "--run":
        by = {r["name"]: r for r in recs}
        for name in sys.argv[2:]:
            rerun(by[name])
    elif sys.argv[1] == "--rational":
        by = {r["name"]: r for r in recs}
        for name in sys.argv[2:]:
            rational(by[name])
    elif sys.argv[1] == "--prime":
        by = {r["name"]: r for r in recs}
        for name in sys.argv[3:]:
            other_prime(by[name], int(sys.argv[2]))
    elif sys.argv[1] == "--run-small":
        for r in recs:
            if "dims" in r and (r.get("peak_MB") or 0) < 2000:
                rerun(r)
