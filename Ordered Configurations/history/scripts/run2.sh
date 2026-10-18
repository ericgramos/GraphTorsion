S=/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad
mkdir -p $S
run() { ./oconf9 -g 1,2,3,4,0 -d 3 -p $2 -o $S/r2 < $1 2>/dev/null | python3 -c "
import json,sys; d=json.load(sys.stdin); c=d['cells']; r=[0]+[m['pivots'] for m in d['maps']]
print('$3 p=$2: cells',c,'ranks',r,'b0..b2',[c[k]-r[k]-r[k+1] for k in range(3)], flush=True)"; }
run apex5.in 3 "K5 with doubled K4 (non-planar)"
run ctrl5.in 2 "doubled K4 + vertex on 3 (planar)"
run ctrl5.in 1000003 "doubled K4 + vertex on 3 (planar)"
run ctrl5.in 3 "doubled K4 + vertex on 3 (planar)"
