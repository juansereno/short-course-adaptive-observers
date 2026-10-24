"""Roadmap: the 4-step recipe of the course."""
import os
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from style import OUT, INK, MUTED, SOFT, EST

steps = [
    ("1", "Model", "What is measured (y)?\nWhat is unknown (x, θ)?", "virus–T-cell model"),
    ("2", "Can we?", "Observability and\nidentifiability", "rank conditions"),
    ("3", "Rewrite", "State-affine form\nẋ = A(y)x + β(y) + φ(y)θ", "sort terms into A, β, φ"),
    ("4", "Design & code", "Adaptive observer (17),\nassumptions A1–A3, tuning", "MATLAB · Python"),
]

fig = plt.figure(figsize=(16, 3.2), facecolor="white")
ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off(); ax.set_xlim(0, 16); ax.set_ylim(0, 3.2)
w, h, gap, x0, y0 = 3.3, 2.2, 0.75, 0.45, 0.55
for i, (num, title, body, foot) in enumerate(steps):
    x = x0 + i * (w + gap)
    ax.add_patch(FancyBboxPatch((x, y0), w, h, boxstyle="round,pad=0.02,rounding_size=0.15",
                                fc=SOFT, ec=INK, lw=1.4))
    ax.add_patch(FancyBboxPatch((x, y0 + h - 0.6), w, 0.6, boxstyle="round,pad=0.02,rounding_size=0.15",
                                fc=INK, ec=INK, lw=1.4))
    ax.text(x + 0.2, y0 + h - 0.3, f"{num} · {title}", color="white", fontsize=15, fontweight="bold", va="center")
    ax.text(x + w / 2, y0 + 0.95, body, color=INK, fontsize=12.5, ha="center", va="center")
    ax.text(x + w / 2, y0 + 0.25, foot, color=MUTED, fontsize=11, ha="center", va="center", style="italic")
    if i < len(steps) - 1:
        ax.annotate("", xy=(x + w + gap - 0.05, y0 + h / 2), xytext=(x + w + 0.05, y0 + h / 2),
                    arrowprops=dict(arrowstyle="-|>", color=EST, lw=2.2, mutation_scale=18))
fig.savefig(os.path.join(OUT, "roadmap.png"), dpi=100)
print("roadmap.png")
