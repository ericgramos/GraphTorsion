S=/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad
for p in 3 13 2; do
  /usr/bin/env python3 - $p <<'PY'
import subprocess, json, sys, time, threading
p = sys.argv[1]; S = "/tmp/claude-0/-home-claude/6317a15a-d6c4-555d-857c-d8fa76cb1885/scratchpad"
proc = subprocess.Popen(["./oconf12", "-g", "1,2,3,4,0", "-d", "3", "-p", p, "-o", S + "/k33c"], stdin=open("k33c.in"), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
peak = [0]
def w():
    while proc.poll() is None:
        try: peak[0] = max(peak[0], int(open(f"/proc/{proc.pid}/status").read().split("VmRSS:")[1].split()[0]) // 1024)
        except Exception: pass
        time.sleep(2)
threading.Thread(target=w, daemon=True).start()
out, _ = proc.communicate(); t = time.time()
if proc.returncode: print("p", p, "FAILED", proc.returncode, "peak", peak[0], flush=True)
else:
    d = json.loads(out); c = d["cells"]; r = [0] + [m["pivots"] for m in d["maps"]]
    print("p", p, "cells", c, "ranks", r[1:], "b0..b2", [c[k] - r[k] - r[k+1] for k in range(3)], "peak", peak[0], flush=True)
PY
done
