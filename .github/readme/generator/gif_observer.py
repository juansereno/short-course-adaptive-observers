"""Animation: the observer learns c_T and reconstructs T cells as data arrive (Example 1)."""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from style import OUT, ACTUAL, EST, INK, MUTED, run

ns = run("ex1_cT.py")
t, V, T, T_hat = ns["t"], ns["V"], ns["T"], ns["T_hat"]
cT, cT_hat = ns["c_T"], ns["c_T_hat"]
lam_min = ns["lam_min_Gamma"]

fig, ax = plt.subplots(1, 3, figsize=(15, 4.0))
fig.subplots_adjust(left=0.05, right=0.99, bottom=0.17, top=0.86, wspace=0.28)

ax[0].set_title("Measured: viral load V(t)")
ax[0].set_yscale("log"); ax[0].set_xlim(0, 21); ax[0].set_ylim(1, 3e9)
ax[1].set_title("Not measured: T cells (×10⁶)")
ax[1].set_xlim(0, 21); ax[1].set_ylim(0, 6)
ax[2].set_title("Unknown parameter $c_T$ (×10⁻⁶)")
ax[2].set_xlim(0, 21); ax[2].set_ylim(0, 2.6)
for a in ax:
    a.set_xlabel("days after symptom onset")

ax[1].plot(t, T / 1e6, color=ACTUAL, alpha=0.25)
ax[2].axhline(cT * 1e6, color=ACTUAL, lw=2, label="actual")
lV, = ax[0].plot([], [], color=ACTUAL, label="data")
lT, = ax[1].plot([], [], color=ACTUAL, label="actual")
lTh, = ax[1].plot([], [], color=EST, ls="--", label="observer")
lc, = ax[2].plot([], [], color=EST, ls="--", label="observer")
ax[0].legend(loc="upper right"); ax[1].legend(loc="upper left"); ax[2].legend(loc="lower right")
txt = fig.text(0.5, 0.95, "", ha="center", color=INK, fontsize=13)

frames = np.r_[np.linspace(0, 8, 50), np.linspace(8.2, 21, 25), [21] * 12]


def update(day):
    k = np.searchsorted(t, day, side="right")
    lV.set_data(t[:k], V[:k])
    lT.set_data(t[:k], T[:k] / 1e6)
    lTh.set_data(t[:k], T_hat[:k] / 1e6)
    lc.set_data(t[:k], cT_hat[:k] * 1e6)
    pe = "no excitation yet: the estimate waits" if day < 2 else (
         "the infection excites the system: learning" if day < 4.5 else "converged: T cells reconstructed")
    txt.set_text(f"day {day:4.1f}   ·   {pe}")
    return lV, lT, lTh, lc, txt


anim = FuncAnimation(fig, update, frames=frames, blit=False)
anim.save(os.path.join(OUT, "observer.gif"), writer=PillowWriter(fps=12), dpi=72)
update(21)
fig.savefig(os.path.join(OUT, "observer_still.png"), dpi=100)
print("observer.gif, observer_still.png")
