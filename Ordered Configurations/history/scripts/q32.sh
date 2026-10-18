#!/bin/bash
S=/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad
cd /home/claude/gconf
while pgrep -x q31.sh > /dev/null; do sleep 30; done
run() { name=$1; shift; inp=$1; shift; s=$(date +%s); rm -f $name.mem; ( while sleep 20; do ps -o rss= -C oconf23 >> $name.mem; done ) & mp=$!; ./oconf23 -F "$@" -o $S/$name < $inp > $name.json 2> $name.log; echo "$name rc=$? secs=$(( $(date +%s) - s )) peakKB=$(sort -n $name.mem | tail -1)" >> q17.done; kill $mp; }
F20="1,2,3,4,0;0,2,4,1,3"
for tw in l le n e; do run k5mT133n6_$tw k5mT133_6.in -q 3 2 -g "$F20" -a "0,1,2,3,4|0,1,2,3,4,5,6,7,8,10,9,11,12;0,1,2,3,4|0,1,2,3,4,5,6,7,8,9,10,12,11" -X $tw -D 70000; done
