#!/bin/bash
S=/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad
cd /home/claude/gconf
while pgrep -x q21.sh > /dev/null; do sleep 30; done
run() { name=$1; shift; inp=$1; shift; s=$(date +%s); ( while sleep 20; do ps -o rss= -C oconf20 >> $name.mem; done ) & mp=$!; ./oconf20 "$@" -o $S/$name < $inp > $name.json 2> $name.log; echo "$name rc=$? secs=$(( $(date +%s) - s )) peakKB=$(sort -n $name.mem | tail -1)" >> q17.done; kill $mp; }
F20="1,2,3,4,0;0,2,4,1,3"
AUT="1,2,3,4,0;0,1,2,3,4|0,1,2,3,4,15,6,7,8,9,10,11,12,13,14,5,16,17,18,19"
for tw in n l e le; do run dk5n5_p3_$tw dk5r_5.in -q 3 3 -g "$F20" -a "$AUT" -X $tw -D 40000; done
for tw in n l e le; do run dk5n6_p3_$tw dk5r_6.in -q 3 3 -g "$F20" -a "$AUT" -X $tw -D 40000; done
