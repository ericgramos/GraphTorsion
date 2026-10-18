#!/bin/bash
S=/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad
cd /home/claude/gconf
run() { name=$1; shift; inp=$1; shift; s=$(date +%s); rm -f $name.mem; ( while sleep 20; do ps -o rss= -C oconf22 >> $name.mem; done ) & mp=$!; ./oconf22 -F "$@" -o $S/$name < $inp > $name.json 2> $name.log; echo "$name rc=$? secs=$(( $(date +%s) - s )) peakKB=$(sort -n $name.mem | tail -1)" >> q17.done; kill $mp; }
F21="1,2,3,4,5,6,0;0,2,4,6,1,3,5"
F20="1,2,3,4,0;0,2,4,1,3"
AUT="1,2,3,4,0;0,1,2,3,4|0,1,2,3,4,15,6,7,8,9,10,11,12,13,14,5,16,17,18,19"
run dk5n6_p3_n dk5r_6.in -q 3 3 -g "$F20" -a "$AUT" -D 40000
run dk5n6_p3_l dk5r_6.in -q 3 3 -g "$F20" -a "$AUT" -X l -D 40000
run k5n7_p2z4_d3 k5r_7.in -q 2 2 -g "$F21" -a "1,2,3,4,0" -d 3 -D 100000
run dk5n5_p3_n dk5r_5.in -q 3 3 -g "$F20" -a "$AUT" -D 40000
run dk5n6_p3_e dk5r_6.in -q 3 3 -g "$F20" -a "$AUT" -X e -D 40000
