#!/bin/bash
S=/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad
cd /home/claude/gconf
while pgrep -f "t133a_n6_F20l" > /dev/null; do sleep 30; done
run() { name=$1; shift; inp=$1; shift; s=$(date +%s); rm -f $name.mem; ( while sleep 20; do ps -o rss= -C oconf23 >> $name.mem; done ) & mp=$!; ./oconf23 -F "$@" -o $S/$name < $inp > $name.json 2> $name.log; echo "$name rc=$? secs=$(( $(date +%s) - s )) peakKB=$(sort -n $name.mem | tail -1)" >> q17.done; kill $mp; }
F20="1,2,3,4,0;0,2,4,1,3"
L40="1,2,3,4,0,5,6;0,2,4,1,3,5,6;0,1,2,3,4,6,5"
for tw in l le; do run k5T133n6_$tw k5T133_6.in -q 3 2 -g "$F20" -a "0,1,2,4,3;0,1,2,3,4|0,1,2,3,4,5,6,7,8,9,11,10,12,13;0,1,2,3,4|0,1,2,3,4,5,6,7,8,9,10,11,13,12" -X $tw -D 70000; done
for tw in l le; do run k33T133n6_$tw k33T133_6.in -q 3 2 -g "$F20" -a "0,1,2,3,5,4;0,1,2,3,4,5|0,1,2,3,4,5,6,7,8,9,11,10,12,13;0,1,2,3,4,5|0,1,2,3,4,5,6,7,8,9,10,11,13,12" -X $tw -D 70000; done
for tw in n e; do run k5T133n6_$tw k5T133_6.in -q 3 2 -g "$F20" -a "0,1,2,4,3;0,1,2,3,4|0,1,2,3,4,5,6,7,8,9,11,10,12,13;0,1,2,3,4|0,1,2,3,4,5,6,7,8,9,10,11,13,12" -X $tw -D 70000; done
for tw in l le; do run t133a_n7_$tw t133a_7.in -q 3 2 -g "$L40" -a "0,1,2,3|0,2,1,3,4,5,6,7,8,9;0,1,2,3|0,1,2,3,5,4,6,7,8,9" -X $tw -D 70000; done
