"""Generate the figures used in the slides by running the two example scripts."""
import os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
PY = os.path.join(HERE, "..", "python")
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

ACTUAL, EST, INK, MUTED = "#2a78d6", "#eb6834", "#2b2b2b", "#8a8a85"
plt.rcParams.update({
    "font.size": 11, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
    "grid.color": "#e6e6e3", "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "lines.linewidth": 2,
})


def run(script, **subs):
    src = open(os.path.join(PY, script)).read().replace("plt.show()", "plt.close('all')")
    for k, v in subs.items():
        src, n = re.subn(rf"^{k}\s*=.*$", f"{k} = {v}", src, flags=re.M)
        assert n == 1, k
    ns = {"__name__": "fig"}
    exec(compile(src, script, "exec"), ns)
    return ns


def mark_pe(ax):
    ax.axvspan(0, 2, color="#f1f0ec", zorder=0)


# ---------------------------------------------------------------- patient
ns = run("ex1_cT.py")
f = lambda t, z: ns["patient"](t, z[0], z[1])
s = solve_ivp(f, [-3, 21], [ns["V0"], ns["T0"]], rtol=1e-8, atol=1e-8, max_step=0.02)
fig, ax = plt.subplots(1, 2, figsize=(9, 3.0))
ax[0].semilogy(s.t, s.y[0], color=ACTUAL)
ax[0].axhline(50, color=MUTED, ls=":", lw=1)
ax[0].text(20.5, 70, "detection limit (50 copies/ml)", ha="right", va="bottom", fontsize=9, color=MUTED)
ax[0].set_title("Virus V(t)  (measured)", fontsize=11, color=INK)
ax[1].plot(s.t, s.y[1] / 1e6, color=ACTUAL)
ax[1].set_title(r"T cells T(t) $\times 10^6$  (NOT measured)", fontsize=11, color=INK)
for a in ax:
    a.axvspan(-3, 0, color="#f1f0ec", zorder=0)
    a.set_xlabel("t [days after symptom onset]")
ax[0].text(-1.5, 1e6, "no data", ha="center", fontsize=9, color=MUTED, rotation=90)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "patient.pdf")); plt.close()

# ---------------------------------------------------------------- example 1
t, T, That, cT, cTh, V, lam = ns["t"], ns["T"], ns["T_hat"], ns["c_T"], ns["c_T_hat"], ns["V"], ns["lam"]
fig, ax = plt.subplots(3, 1, figsize=(6.2, 5.6), sharex=True)
ax[0].plot(t, cT * np.ones_like(t) * 1e6, color=ACTUAL, label="actual")
ax[0].plot(t, cTh * 1e6, color=EST, ls="--", label="estimated")
ax[0].set_ylabel(r"$c_T\ (\times10^{-6})$"); ax[0].legend(loc="lower right")
ax[1].plot(t, T / 1e6, color=ACTUAL, label="actual")
ax[1].plot(t, That / 1e6, color=EST, ls="--", label="estimated")
ax[1].set_ylabel(r"$T\ (\times10^{6})$"); ax[1].set_ylim(0, 6); ax[1].legend(loc="upper right")
ax[2].plot(t, -V / lam, color=INK)
ax[2].set_ylabel(r"$\varphi_{21}=-y_1$"); ax[2].set_xlabel("t [days]")
ax[2].annotate("regressor ~ 0:\nnothing to learn yet", xy=(1.0, -0.1), xytext=(12, -6),
               fontsize=9, color=MUTED, arrowprops=dict(arrowstyle="->", color=MUTED))
for a in ax: mark_pe(a)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "ex1.pdf")); plt.close()

# ---------------------------------------------------------------- example 2
n2 = run("ex2_cT_r.py")
t, T, That = n2["t"], n2["T"], n2["T_hat"]
fig, ax = plt.subplots(3, 1, figsize=(6.2, 5.6), sharex=True)
ax[0].plot(t, n2["c_T"] * np.ones_like(t) * 1e6, color=ACTUAL, label="actual")
ax[0].plot(t, n2["c_T_hat"] * 1e6, color=EST, ls="--", label="estimated")
ax[0].set_ylabel(r"$c_T\ (\times10^{-6})$"); ax[0].legend(loc="lower right")
ax[1].plot(t, n2["r"] * np.ones_like(t), color=ACTUAL, label="actual")
ax[1].plot(t, n2["r_hat"], color=EST, ls="--", label="estimated")
ax[1].set_ylabel(r"$r$"); ax[1].legend(loc="lower right")
ax[2].plot(t, T / 1e6, color=ACTUAL, label="actual")
ax[2].plot(t, That / 1e6, color=EST, ls="--", label="estimated")
ax[2].set_ylabel(r"$T\ (\times10^{6})$"); ax[2].set_ylim(0, 6); ax[2].set_xlabel("t [days]")
ax[2].legend(loc="upper right")
plt.tight_layout(); plt.savefig(os.path.join(OUT, "ex2.pdf")); plt.close()

# ---------------------------------------------------------------- noise
a = run("ex2_cT_r.py", ADD_NOISE="True")
b = run("ex2_cT_r.py", ADD_NOISE="True", rho_theta="1.0", rho_x="10.0")
fig, ax = plt.subplots(1, 2, figsize=(9, 2.9), sharey=True)
for axi, nn, ttl in [(ax[0], a, r"$\rho_x=100,\ \rho_\theta=4$ (default)"),
                     (ax[1], b, r"$\rho_x=10,\ \rho_\theta=1$ (slower, smoother)")]:
    axi.plot(nn["t"], nn["c_T"] * np.ones_like(nn["t"]) * 1e6, color=ACTUAL, label="actual")
    axi.plot(nn["t"], nn["c_T_hat"] * 1e6, color=EST, ls="--", label="estimated")
    axi.set_title(ttl, fontsize=11, color=INK); axi.set_xlabel("t [days]")
    axi.set_ylim(-5, 6)
ax[0].set_ylabel(r"$c_T\ (\times10^{-6})$"); ax[0].legend(loc="lower left")
plt.tight_layout(); plt.savefig(os.path.join(OUT, "noise.pdf")); plt.close()

# ---------------------------------------------------------------- PE meter
e1, e2 = run("ex1_cT.py"), run("ex2_cT_r.py")
fig, ax = plt.subplots(1, 2, figsize=(9, 3.1), sharey=True)
for axi, nn, ttl, lab in [(ax[0], e1, r"Ex. 1  ($q=1$: $\Gamma$ is a scalar)", r"$\Gamma(t)$"),
                          (ax[1], e2, r"Ex. 2  ($q=2$: smallest eigenvalue)", r"$\lambda_{\min}(\Gamma(t))$")]:
    tt = nn["t"]
    axi.axvspan(1.8, 4.0, color="#f1f0ec", zorder=0)
    axi.semilogy(tt, nn["lam_min_Gamma"], color=EST, label=lab)
    axi.semilogy(tt, nn["no_info"], color=MUTED, ls=":", lw=1.5, label=r"no new data: $e^{-\rho_\theta t}\,\Gamma(0)$")
    axi.set_title(ttl, fontsize=11, color=INK); axi.set_xlabel("t [days]")
    axi.set_ylim(1e-24, 10); axi.legend(loc="lower left", fontsize=9)
ax[0].annotate("rises: data excite $\\theta$", xy=(3.0, 1e-3), xytext=(6.5, 1e-15),
               fontsize=9, color=INK, arrowprops=dict(arrowstyle="->", color=MUTED))
ax[0].annotate("parallel to dotted line:\nonly forgetting, PE lost", xy=(13, 1e-12), xytext=(15.5, 1e-4),
               fontsize=9, color=INK, ha="center", arrowprops=dict(arrowstyle="->", color=MUTED))
ax[1].annotate("weaker excitation of\nthe 2nd direction", xy=(5, 1.5e-6), xytext=(13, 1e-2),
               fontsize=9, color=INK, ha="center", arrowprops=dict(arrowstyle="->", color=MUTED))
ax[0].set_ylabel(r"$\lambda_{\min}(\Gamma)$")
plt.tight_layout(); plt.savefig(os.path.join(OUT, "pe_meter.pdf")); plt.close()
print("figures written to", OUT)
