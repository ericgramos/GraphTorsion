#!/bin/bash
S=/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad
cd /home/claude/gconf
while pgrep -x q18.sh > /dev/null; do sleep 30; done
run() { name=$1; shift; inp=$1; shift; s=$(date +%s); ( while sleep 20; do ps -o rss= -C oconf18 >> $name.mem; done ) & mp=$!; ./oconf18 "$@" -o $S/$name < $inp > $name.json 2> $name.log; echo "$name rc=$? secs=$(( $(date +%s) - s )) peakKB=$(sort -n $name.mem | tail -1)" >> q17.done; kill $mp; }
run dk5_z5_p2_d3 dk5_5.in -q 2 4 -g "1,2,3,4,0" -d 3 -D 40000
run k5de_z5_p2_d3 k5de_6.in -q 2 4 -g "1,2,3,4,0,5" -d 3 -D 40000
