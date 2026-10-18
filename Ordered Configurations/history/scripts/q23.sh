#!/bin/bash
S=/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad
cd /home/claude/gconf
while pgrep -x q22.sh > /dev/null || pgrep -x q21.sh > /dev/null; do sleep 30; done
run() { name=$1; shift; inp=$1; shift; s=$(date +%s); rm -f $name.mem; ( while sleep 20; do ps -o rss= -C oconf21 >> $name.mem; done ) & mp=$!; ./oconf21 "$@" -o $S/$name < $inp > $name.json 2> $name.log; echo "$name rc=$? secs=$(( $(date +%s) - s )) peakKB=$(sort -n $name.mem | tail -1)" >> q17.done; kill $mp; }
F21="1,2,3,4,5,6,0;0,2,4,6,1,3,5"
run k5n7_p2z4_d3 k5r_7.in -q 2 2 -g "$F21" -a "1,2,3,4,0" -d 3 -D 100000
