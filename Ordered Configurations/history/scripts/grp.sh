S=/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad
run() { for p in 13 3; do ./oconf13 -g "$2" $4 -p $p -o $S/grp < $1 2>/dev/null | python3 -c "
import json,sys; d=json.load(sys.stdin); c=d['cells']; r=[0]+[m['pivots'] for m in d['maps']]+[0]
print('$3 p=$p: betti',[c[k]-r[k]-r[k+1] for k in range(len(c))], flush=True)"; done; }
run c4t5.in "1,2,3,4,0;0,4,3,2,1" "tripled C4 n=5, D5 (order 10)"
run c4t5.in "0,2,4,1,3" "tripled C4 n=5, Z4 (order 4)"
run c4t5.in "1,0,2,3,4;0,1,3,2,4;2,3,0,1,4" "tripled C4 n=5, D4 on 0-3 (order 8)"
run c4t5.in "1,0,2,3,4" "tripled C4 n=5, single transposition (order 2)"
for p in 13 3 2; do ./oconf13 -g "1,2,3,4,0" -p $p -o $S/grp < c4tv5.in 2>/dev/null | python3 -c "
import json,sys; d=json.load(sys.stdin); c=d['cells']; r=[0]+[m['pivots'] for m in d['maps']]+[0]
print('CONTROL tripled C4 + degree-3 vertex n=5, Z5, p=$p: betti',[c[k]-r[k]-r[k+1] for k in range(len(c))], flush=True)"; done
