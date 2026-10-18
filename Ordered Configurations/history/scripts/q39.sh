#!/bin/bash
S=/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad
cd /home/claude/gconf
while pgrep -x q38.sh > /dev/null || pgrep -x q35.sh > /dev/null; do sleep 30; done
run() { name=$1; shift; inp=$1; shift; s=$(date +%s); rm -f $name.mem; ( while sleep 20; do ps -o rss= -C oconf26 >> $name.mem; done ) & mp=$!; ./oconf26 -F "$@" -o $S/$name < $inp > $name.json 2> $name.log; echo "$name rc=$? secs=$(( $(date +%s) - s )) peakKB=$(sort -n $name.mem | tail -1)" >> q17.done; kill $mp; }
run k5T222r_n6_w k5T222r_6.in -g "1,2,0,3,4,5;0,1,2,4,5,3" -a "1,2,0,3,4" -W 1 -d 3
