S=/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad
while pgrep -f run2.sh > /dev/null; do sleep 5; done
for p in 2 1000003; do ./oconf9 -g 1,2,3,4,0 -d 3 -p $p -o $S/r3 < leaf5.in 2>/dev/null | python3 -c "
import json,sys; d=json.load(sys.stdin); c=d['cells']; r=[0]+[m['pivots'] for m in d['maps']]
print('doubled K4 + one leaf p=$p: ranks',r,'b0..b2',[c[k]-r[k]-r[k+1] for k in range(3)], flush=True)"; done
