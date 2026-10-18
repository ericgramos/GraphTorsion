#!/bin/bash
S=/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad
cd /home/claude/gconf
while pgrep -x q34.sh > /dev/null; do sleep 30; done
run() { name=$1; shift; inp=$1; shift; s=$(date +%s); rm -f $name.mem; ( while sleep 20; do ps -o rss= -C oconf23 >> $name.mem; done ) & mp=$!; ./oconf23 -F "$@" -o $S/$name < $inp > $name.json 2> $name.log; echo "$name rc=$? secs=$(( $(date +%s) - s )) peakKB=$(sort -n $name.mem | tail -1)" >> q17.done; kill $mp; }
Z33="1,2,0,3,4,5;0,1,2,4,5,3"
run t123a_n6_p2 t123a_6.in -q 2 2 -g "$Z33" -d 3 -D 100000
run t222a_n6_p2 t222a_6.in -q 2 2 -g "$Z33" -d 3 -D 100000
run t223a_n6_p2 t223a_6.in -q 2 2 -g "$Z33" -d 3 -D 100000
F20="1,2,3,4,0;0,2,4,1,3"
for tw in l n; do run t223a_n6_p3_$tw t223a_6.in -q 3 2 -g "$F20" -X $tw -d 3 -D 70000; done
L40="1,2,3,4,0,5,6;0,2,4,1,3,5,6;0,1,2,3,4,6,5"
for tw in n e; do run t133a_n7_$tw t133a_7.in -q 3 2 -g "$L40" -a "0,1,2,3|0,2,1,3,4,5,6,7,8,9;0,1,2,3|0,1,2,3,5,4,6,7,8,9" -X $tw -D 70000; done
