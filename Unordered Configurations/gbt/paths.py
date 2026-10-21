"""
Where the repository keeps its data, resolved from this file's location.

Every experiment writes its results into `<repo root>/data/` and reads them back from
there, so the scripts can be launched from any working directory:

    python3 experiments/sweep_complete.py                 # from the repo root
    python3 /path/to/repo/experiments/sweep_complete.py    # from anywhere else

both put the JSON in the same place.

FILE FORMAT
-----------
Every data file is one JSON document written by `save_json` with sorted keys and one-space
indentation, so successive versions diff cleanly.  A homology dict is stored with string
degrees and lists instead of tuples (JSON has no tuples and no integer keys); `record`
does that conversion:

    >>> from gbt.paths import record
    >>> record({1: (6, [2]), 2: (0, [])})
    {'1': [6, [2]], '2': [0, []]}

To read it back as a homology dict:  {int(d): (v[0], v[1]) for d, v in stored.items()}.
"""

import json
import os

#: absolute path of the repository root (the directory containing `gbt/` and `data/`)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: absolute path of the data directory; created on first use
DATA = os.path.join(ROOT, "data")


def data_path(name):
    """Absolute path of `name` (a file name such as 'sweep_complete.json') inside the data
    directory.  Creates the data directory if it does not exist yet."""
    os.makedirs(DATA, exist_ok=True)
    return os.path.join(DATA, name)


def load_json(name, default=None):
    """Read a JSON file from the data directory, returning `default` if it is absent.

    `name` is a file name such as 'sweep_complete.json'.  If the file does not exist the
    result is `default`, or a new empty dict when `default` is None -- which is what lets
    every sweep start from scratch or resume with the same line of code.
    """
    p = os.path.join(DATA, name)
    if not os.path.exists(p):
        return {} if default is None else default
    with open(p) as fh:
        return json.load(fh)


def save_json(name, obj):
    """Write `obj` to the data directory as JSON, with a stable key order.

    Sweeps call this after every completed job, so that an interrupted run can be resumed
    from the file rather than restarted.  The whole file is rewritten each time (it is
    small), so a script that saves a dict it did not first load will REPLACE the stored
    results: always load, update, save.
    """
    with open(data_path(name), "w") as fh:
        json.dump(obj, fh, sort_keys=True, indent=1)


def record(H):
    """Serialise a homology dict {d: (free, invariant factors)} for JSON.

    Returns {str(d): [free, list(invariant factors)]}; see the module docstring.
    """
    return {str(d): [H[d][0], list(H[d][1])] for d in H}
