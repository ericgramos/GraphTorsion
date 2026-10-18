import json, subprocess, sys, itertools
S = "/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad"
def loc(edges, n, p, k, extra=()):
    V = 1 + max(max(e) for e in edges)
    inp = f"{V} {len(edges)} {n}\n" + "".join(f"{a} {b}\n" for a, b in edges)
    out = subprocess.run(["./oconf14", "-q", str(p), str(k), "-o", S + "/loc14", *extra], input=inp, capture_output=True, text=True)
    if out.returncode: return "FAILED %d" % out.returncode
    d = json.loads(out.stdout)
    tors = {}
    for m in d["maps"]:
        st = m["stages"]
        t = [p ** s for s, c in enumerate(st) for _ in range(c) if s > 0]
        if t or m["leftover_nnz"]: tors[f"H{m['k']-1}"] = {"orders": dict((x, t.count(x)) for x in set(t)), "leftover": m["leftover_nnz"]}
    return tors
if __name__ == "__main__":
    K = lambda m: list(itertools.combinations(range(m), 2))
    T = lambda a, b, c: [(0,1)]*a + [(0,2)]*b + [(1,2)]*c
    C4 = [(0,1),(1,2),(2,3),(3,0)]
    tests = [("UConf_2(K5)  expect H1: Z/2", K(5), 2, 2, 3, ("-k", "aa")),
             ("UConf_4(K6)  expect H1, H2: Z/2", K(6), 4, 2, 3, ("-k", "aaaa")),
             ("UConf_4(K6) mod 3^3  expect none", K(6), 4, 3, 3, ("-k", "aaaa")),
             ("doubled K4 n=5 /Z5  expect 9 two-primary summands in H2", K(4)*2, 5, 2, 4, ("-g", "1,2,3,4,0")),
             ("tripled C4 n=5 /Z5  expect (Z/3)^3 in H2", C4*3, 5, 3, 3, ("-g", "1,2,3,4,0")),
             ("T(1,2,3) n=6 /Z5  expect Z/2 in H2", T(1,2,3), 6, 2, 4, ("-g", "1,2,3,4,0,5")),
             ("T(2,2,2) n=6 full  expect (Z/2)^2", T(2,2,2), 6, 2, 4, ()),
             ("T(1,2,3) n=6 full  expect Z/2", T(1,2,3), 6, 2, 4, ()),
             ("T(1,3,3) n=6 full  expect (Z/3)^6", T(1,3,3), 6, 3, 3, ()),
             ("T(2,2,3) n=6 full mod 2^4  expect (Z/2)^12", T(2,2,3), 6, 2, 4, ()),
             ("T(2,2,3) n=6 full mod 3^3  expect (Z/3)^8", T(2,2,3), 6, 3, 3, ())]
    for name, g, n, p, k, extra in tests:
        import time; t = time.time()
        print(name, "->", loc(g, n, p, k, extra), f"[{time.time()-t:.0f}s]", flush=True)
