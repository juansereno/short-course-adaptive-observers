"""Highlights strip: two parameters, the PE meter, and noise robustness."""
import os
import numpy as np
import matplotlib.pyplot as plt
from style import OUT, ACTUAL, EST, INK, MUTED, run

e2 = run("ex2_cT_r.py")
n_def = run("ex2_cT_r.py", ADD_NOISE="True")
n_tun = run("ex2_cT_r.py", ADD_NOISE="True", rho_theta="1.0", rho_x="10.0")

fig, ax = plt.subplots(1, 3, figsize=(15, 3.9))
fig.subplots_adjust(left=0.05, right=0.99, bottom=0.17, top=0.86, wspace=0.28)

t = e2["t"]
ax[0].axhline(1, color=MUTED, lw=1)
ax[0].plot(t, e2["c_T_hat"] / e2["c_T"], color=EST, label=r"$\hat c_T / c_T$")
ax[0].plot(t, e2["r_hat"] / e2["r"], color=ACTUAL, ls="--", label=r"$\hat r / r$")
ax[0].set_ylim(0, 1.6); ax[0].legend(loc="lower right")
ax[0].set_title("Example 2: two parameters at once")

ax[1].semilogy(t, e2["lam_min_Gamma"], color=EST, label=r"$\lambda_{\min}(\Gamma(t))$")
ax[1].semilogy(t, e2["no_info"], color=MUTED, ls=":", label="no new data")
ax[1].set_ylim(1e-24, 10); ax[1].legend(loc="lower left")
ax[1].set_title("PE meter: smallest eigenvalue of Γ")

ax[2].plot(n_def["t"], n_def["c_T_hat"] * 1e6, color=MUTED, lw=1.4, label=r"$\rho_x=100,\ \rho_\theta=4$")
ax[2].plot(n_tun["t"], n_tun["c_T_hat"] * 1e6, color=EST, label=r"$\rho_x=10,\ \rho_\theta=1$")
ax[2].axhline(e2["c_T"] * 1e6, color=ACTUAL, lw=2, label="actual")
ax[2].set_ylim(-5, 6); ax[2].legend(loc="lower left", fontsize=10)
ax[2].set_title("5 % noise: re-tuning the forgetting factors")
ax[2].set_ylabel(r"$\hat c_T$ (×10⁻⁶)")
for a in ax:
    a.set_xlabel("days after symptom onset")
fig.savefig(os.path.join(OUT, "highlights.png"), dpi=100)
print("highlights.png")
