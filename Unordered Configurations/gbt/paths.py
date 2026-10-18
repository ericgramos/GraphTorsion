"""
Where the repository keeps its data, resolved from this file's location.

Every experiment writes its results into `<repo root>/data/` and reads them back from
there, so the scripts can be launched from any working directory:

    python3 experiments/sweep_complete.py                 # from the repo root
    python3 /path/to/repo/experiments/sweep_complete.py    # from anywhere else

both put the JSON in the same place.
"""

import json
import os

#: absolute path of the repository root (the directory containing `gbt/` and `data/`)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: absolute path of the data directory; created on first use
DATA = os.path.join(ROOT, "data")


def data_path(name):
    """Absolute path of `name` inside the data directory, creating the directory."""
    os.makedirs(DATA, exist_ok=True)
    return os.path.join(DATA, name)


def load_json(name, default=None):
    """Read a JSON file from the data directory, returning `default` if it is absent."""
    p = os.path.join(DATA, name)
    if not os.path.exists(p):
        return {} if default is None else default
    with open(p) as fh:
        return json.load(fh)


def save_json(name, obj):
    """Write `obj` to the data directory as JSON, with a stable key order.

    Sweeps call this after every completed job, so that an interrupted run can be resumed
    from the file rather than restarted.
    """
    with open(data_path(name), "w") as fh:
        json.dump(obj, fh, sort_keys=True, indent=1)


def record(H):
    """Serialise a homology dict {d: (free, invariant factors)} for JSON."""
    return {str(d): [H[d][0], list(H[d][1])] for d in H}
