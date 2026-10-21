"""_common.py -- shared helpers for the scripts in surveys/.

* puts ../python on the module path (so ``from gconf import homology`` works) and defines ENGINE, the path
  of the compiled engine (environment variable OCONF_BIN overrides it);
* graph constructors used throughout (complete graphs, complete bipartite graphs, triangle multigraphs
  T(a,b,c), cycles, Petersen graph);
* graph_input(edges, n): the engine's input text;
* run_engine(options, edges, n): runs the engine, returns (parsed JSON or None, seconds, peak memory in MB,
  return code), measuring the peak resident memory from /proc (Linux).
"""
import itertools, json, os, subprocess, sys, tempfile, threading, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python"))
ENGINE = os.environ.get("OCONF_BIN", os.path.join(HERE, "..", "engine", "oconf"))
SCRATCH = os.environ.get("GCONF_SCRATCH", tempfile.gettempdir())

K = lambda m: list(itertools.combinations(range(m), 2))                       # complete graph K_m
Kb = lambda a, b: [(i, a + j) for i in range(a) for j in range(b)]           # complete bipartite K_{a,b}
T = lambda a, b, c: [(0, 1)] * a + [(0, 2)] * b + [(1, 2)] * c                # triangle multigraph T(a,b,c)
C = lambda m: [(i, (i + 1) % m) for i in range(m)]                            # cycle of length m
PETERSEN = [(i, (i + 1) % 5) for i in range(5)] + [(i, i + 5) for i in range(5)] + \
           [(5 + i, 5 + (i + 2) % 5) for i in range(5)]


def graph_input(edges, n):
    """Engine input: first line 'V E n', then one 'u v' line per edge (vertices must be 0..V-1)."""
    V = 1 + max(max(e) for e in edges)
    return f"{V} {len(edges)} {n}\n" + "".join(f"{a} {b}\n" for a, b in edges)


def run_engine(options, edges, n, tag="run"):
    """Run the engine with the given option list on (edges, n).

    Returns (data, seconds, peak_MB, returncode); data is the parsed JSON (None on failure)."""
    args = [ENGINE, *options, "-o", os.path.join(SCRATCH, f"{tag}{os.getpid()}")]
    t = time.time()
    p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    peak = [0]

    def watch():                                   # sample VmRSS twice a second
        """Background thread: sample the resident memory of the engine process twice a second and keep the peak."""
        while p.poll() is None:
            try:
                peak[0] = max(peak[0], int(open(f"/proc/{p.pid}/status").read().split("VmRSS:")[1].split()[0]) // 1024)
            except Exception:
                pass
            time.sleep(0.5)

    threading.Thread(target=watch, daemon=True).start()
    out, err = p.communicate(graph_input(edges, n))
    data = json.loads(out) if p.returncode == 0 else None
    return data, time.time() - t, peak[0], p.returncode
