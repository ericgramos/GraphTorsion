#!/bin/bash
S=/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad
cd /home/claude/gconf
run() { name=$1; shift; inp=$1; shift; s=$(date +%s); ( while sleep 20; do ps -o rss= -C oconf19 >> $name.mem; done ) & mp=$!; ./oconf19 "$@" -o $S/$name < $inp > $name.json 2> $name.log; echo "$name rc=$? secs=$(( $(date +%s) - s )) peakKB=$(sort -n $name.mem | tail -1)" >> q17.done; kill $mp; }
F21="1,2,3,4,5,6,0;0,2,4,6,1,3,5"
SYL2="1,2,3,0,4,5,6;2,1,0,3,4,5,6;0,1,2,3,5,4,6"
run k33n7_p2 k33_7.in -q 2 4 -g "$F21" -a "1,2,0,4,5,3" -D 40000
run k33n7_p3 k33_7.in -q 3 3 -g "$SYL2" -a "1,0,2,3,4,5;0,1,2,4,3,5;3,4,5,0,1,2" -D 40000
run k33n7_p3_tw k33_7.in -q 3 3 -g "$SYL2" -a "1,0,2,3,4,5;0,1,2,4,3,5;3,4,5,0,1,2" -X le -D 40000
run k5n7_p2_d3 k5r_7.in -q 2 4 -g "$F21" -a "1,2,3,4,0" -d 3 -D 40000
run dk5n5_p2_d3 dk5r_5.in -q 2 4 -g "1,2,3,4,0" -a "1,2,3,4,0" -d 3 -D 40000
