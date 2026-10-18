import itertools, math, json
from collections import Counter
K5 = list(itertools.combinations(range(5), 2))
def canon(S):
    return min(tuple(sorted(tuple(sorted((p[a], p[b]))) for a, b in S)) for p in itertools.permutations(range(5)))
classes = {}
for r in range(0, 11):
    for S in itertools.combinations(K5, r):
        c = canon(S); classes.setdefault(c, S)
def counts(edges, n):
    deg = Counter(x for e in edges for x in e); E = len(edges); es = [1]
    for v in range(5):
        a = deg[v]-1; es = [(es[i] if i < len(es) else 0)+(a*es[i-1] if i >= 1 else 0) for i in range(len(es)+1)]
    return [math.factorial(n)*es[k]*math.comb(n-k+E-1, E-1) for k in range(min(n,5)+1)]
rows = []
for c, S in classes.items():
    g = K5 + list(S)
    cc = counts(g, 5)
    sub = Counter(x for e in S for x in e)
    rows.append(dict(doubled=list(S), ndoubled=len(S), deg_le3=sum(cc[:4]), q=sum(cc[:4])//5, edges=g))
rows.sort(key=lambda r: r['q'])
json.dump(rows, open('k5sub.json', 'w'))
for r in rows: print(r['ndoubled'], f"{r['deg_le3']/1e6:6.1f}M full  {r['q']/1e6:5.1f}M /Z5", r['doubled'])
print(len(rows), 'classes')
