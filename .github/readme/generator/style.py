"""Shared colours and helpers for the README figures."""
import os
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.abspath(os.path.join(HERE, ".."))
PY = os.path.join(ROOT, "python")

ACTUAL, EST = "#2a78d6", "#eb6834"            # actual (blue) vs estimated (orange)
INK, MUTED, GRID, SOFT = "#2b2b2b", "#8a8a85", "#e6e6e3", "#f4f3ef"
C_A, C_B, C_P = "#2a78d6", "#14875e", "#d2541f"  # A(y), beta(y), phi(y)theta

plt.rcParams.update({
    "font.size": 12, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "lines.linewidth": 2.2,
    "axes.titlesize": 13, "axes.titlecolor": INK,
})


def run(script, **subs):
    """Run one of the course examples and return its variables."""
    src = open(os.path.join(PY, script)).read().replace("plt.show()", "plt.close('all')")
    src = src.replace("print(", "(lambda *a, **k: None)(")
    for k, v in subs.items():
        src, n = re.subn(rf"^{k}\s*=.*$", f"{k} = {v}", src, flags=re.M)
        assert n == 1, k
    ns = {"__name__": "readme"}
    exec(compile(src, script, "exec"), ns)
    return ns
