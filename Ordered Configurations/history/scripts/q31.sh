#!/bin/bash
S=/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad
cd /home/claude/gconf
while pgrep -x q30.sh > /dev/null; do sleep 30; done
run() { name=$1; shift; inp=$1; shift; s=$(date +%s); rm -f $name.mem; ( while sleep 20; do ps -o rss= -C oconf23 >> $name.mem; done ) & mp=$!; ./oconf23 -F "$@" -o $S/$name < $inp > $name.json 2> $name.log; echo "$name rc=$? secs=$(( $(date +%s) - s )) peakKB=$(sort -n $name.mem | tail -1)" >> q17.done; kill $mp; }
L40="1,2,3,4,0,5,6;0,2,4,1,3,5,6;0,1,2,3,4,6,5"
run t133a_n7_l t133a_7.in -q 3 2 -g "$L40" -a "0,1,2,3|0,2,1,3,4,5,6,7,8,9;0,1,2,3|0,1,2,3,5,4,6,7,8,9" -X l -D 70000
