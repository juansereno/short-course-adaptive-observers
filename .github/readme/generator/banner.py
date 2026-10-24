"""Banner: title on the left, the observer at work on the right."""
import os
import numpy as np
import matplotlib.pyplot as plt
from style import OUT, EST, INK, run

ns = run("ex1_cT.py")
t, T, T_hat = ns["t"], ns["T"], ns["T_hat"]

fig = plt.figure(figsize=(16, 4.4), facecolor=INK)
ax = fig.add_axes([0, 0, 1, 1], facecolor=INK)
ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)

ax.text(0.04, 0.70, "Parameter Monitoring with", color="white", fontsize=34, fontweight="bold", va="center")
ax.text(0.04, 0.53, "Adaptive Observers", color="white", fontsize=34, fontweight="bold", va="center")
ax.plot([0.04, 0.50], [0.40, 0.40], color=EST, lw=3)
ax.text(0.04, 0.30, "Estimating T cells and patient-specific parameters in SARS-CoV-2",
        color="#d8d7d2", fontsize=15, va="center")
ax.text(0.04, 0.17, "Mini-course  ·  Juan E. Sereno  ·  MATLAB & Python  ·  2026",
        color="#a9a8a2", fontsize=13, va="center")

# right: actual vs estimated T cells
inset = fig.add_axes([0.60, 0.16, 0.36, 0.70], facecolor=INK)
inset.plot(t, T / 1e6, color="#86b6ef", lw=3, label="T cells (actual, not measured)")
inset.plot(t, T_hat / 1e6, color=EST, lw=3, ls="--", label="adaptive observer")
inset.set_ylim(0, 6.2); inset.set_xlim(0, 21)
inset.grid(False)
for s in ["left", "bottom"]:
    inset.spines[s].set_color("#6b6b66")
inset.tick_params(colors="#a9a8a2", labelsize=11)
inset.set_xlabel("days after symptom onset", color="#a9a8a2", fontsize=11)
inset.set_ylabel(r"T  ($\times 10^6$)", color="#a9a8a2", fontsize=11)
leg = inset.legend(loc="upper right", fontsize=11, labelcolor="#e8e7e2")

fig.savefig(os.path.join(OUT, "banner.png"), dpi=110, facecolor=INK)
print("banner.png")
