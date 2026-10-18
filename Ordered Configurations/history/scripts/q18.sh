#!/bin/bash
S=/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad
cd /home/claude/gconf
while pgrep -x q17.sh > /dev/null || pgrep oconf17 > /dev/null; do sleep 30; done
run() { name=$1; shift; inp=$1; shift; s=$(date +%s); ./oconf18 "$@" -o $S/$name < $inp > $name.json 2> $name.log; echo "$name rc=$? secs=$(( $(date +%s) - s ))" >> q17.done; }
run dk33_z5_p2_d3 dk33_5.in -q 2 4 -g "1,2,3,4,0" -d 3 -D 40000
run dk33_d5_p3_d3 dk33_5.in -q 3 3 -g "1,2,3,4,0;0,4,3,2,1" -d 3 -D 40000
run k5_z5_p2_full k5_6.in -q 2 4 -g "1,2,3,4,0,5" -D 40000
