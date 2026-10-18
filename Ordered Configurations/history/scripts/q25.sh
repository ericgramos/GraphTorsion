#!/bin/bash
S=/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad
cd /home/claude/gconf
while pgrep -x q24.sh > /dev/null; do sleep 30; done
run() { name=$1; shift; inp=$1; shift; s=$(date +%s); rm -f $name.mem; ( while sleep 20; do ps -o rss= -C oconf22 >> $name.mem; done ) & mp=$!; ./oconf22 -F "$@" -o $S/$name < $inp > $name.json 2> $name.log; echo "$name rc=$? secs=$(( $(date +%s) - s )) peakKB=$(sort -n $name.mem | tail -1)" >> q17.done; kill $mp; }
AUT="1,0,2,3,4,5;0,1,2,4,3,5;3,4,5,0,1,2;0,1,2,3,4,5|9,1,2,3,4,5,6,7,8,0,10,11,12,13,14,15,16,17"
run dk33n6_p3_n dk33_6.in -q 3 3 -g "1,2,3,4,0;0,2,4,1,3" -a "$AUT" -D 40000
run dk33n6_p3_l dk33_6.in -q 3 3 -g "1,2,3,4,0;0,2,4,1,3" -a "$AUT" -X l -D 40000
