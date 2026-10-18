#!/bin/bash
S=/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad
cd /home/claude/gconf
while pgrep -x q32.sh > /dev/null; do sleep 30; done
run() { name=$1; shift; inp=$1; shift; s=$(date +%s); rm -f $name.mem; ( while sleep 20; do ps -o rss= -C oconf23 >> $name.mem; done ) & mp=$!; ./oconf23 -F "$@" -o $S/$name < $inp > $name.json 2> $name.log; echo "$name rc=$? secs=$(( $(date +%s) - s )) peakKB=$(sort -n $name.mem | tail -1)" >> q17.done; kill $mp; }
run dk4a3n5_p2 dk4a3_5.in -q 2 2 -g "1,2,3,4,0" -d 3 -D 100000
run dk4a4n5_p2 dk4a4_5.in -q 2 2 -g "1,2,3,4,0" -d 3 -D 100000
