"""
===========================================================================
 Mini-course: Adaptive observers for immune digital twins
 EXAMPLE 2 -- Estimate TWO parameters: the T-cell killing rate c_T and
             the T-cell activation rate r (and reconstruct T cells)

 Based on: Rodriguez-Linan, Sereno, Quiroz-Compean, Hernandez-Vargas,
           "Adaptive Observers for Calibrating Immune Digital Twins in
            Viral Infections", IEEE CDC 2025.  Eqs. (1)-(2), (13)-(17).

 Same structure as Example 1. What changes (look for  <-- NEW):
   * theta = [gamma; r]  with gamma = c_T * S_T   (q = 2 unknowns)
   * the r-terms move from A(y) into phi(y)
   * Lambda becomes 2x2, Gamma becomes 2x2, theta_hat has 2 entries
===========================================================================
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# =========================================================================
# STEP 1. Parameters of the virtual patient ("Mean" patient)
# =========================================================================
p     = 8.57      # virus replication rate                 [1/day]
K     = 2.35e9    # carrying capacity of the virus         [copies/ml]
c     = 2.4       # viral clearance rate                   [1/day]
c_T   = 1.89e-6   # killing rate of virus by T cells       <-- UNKNOWN
delta_T = 0.1     # death rate of T cells                  [1/day]
k_T   = 1.26e5    # half-saturation constant of T cells
m     = 2         # Hill coefficient (dimensionless)
r     = 0.194     # T-cell activation/proliferation rate   <-- UNKNOWN (NEW)
T0    = 1e6       # initial number of T cells
S_T   = delta_T * T0   # homeostatic replenishment rate (= 1e5)
V0    = 0.31      # initial viral load at t = -3 days

alpha = p - c          # net viral growth rate
sigma = p / K          # logistic term
gamma = c_T * S_T      # unknown #1
lam   = 1e8            # scaling factor lambda


# =========================================================================
# STEP 2. The in-host model (Eqs. 1-2): our "virtual patient"
# =========================================================================
def G(V):
    """T-cell activation function G(V) = V^m / (V^m + k_T^m)."""
    return V**m / (V**m + k_T**m)


def patient(t, V, T):
    """Right-hand side of Eqs. (1)-(2). Returns dV/dt and dT/dt."""
    dV = p * V * (1 - V / K) - c_T * V * T - c * V
    dT = S_T + r * T * G(V) - delta_T * T
    return dV, dT


# =========================================================================
# STEP 3. Before symptoms (t = -3 ... 0 days): no measurements
# =========================================================================
t_start, t_obs, t_end = -3.0, 0.0, 21.0

pre = solve_ivp(lambda t, z: patient(t, z[0], z[1]),
                [t_start, t_obs], [V0, T0], rtol=1e-8, atol=1e-8)
V_at_0 = pre.y[0, -1]
T_at_0 = pre.y[1, -1]
print(f"At t = {t_obs:g} (first measurement): V = {V_at_0:.3e}, T = {T_at_0:.3e}")


# =========================================================================
# STEP 4. Measurements y1 = V/lambda, y2 = (dV/dt)/lambda (+ optional noise)
# =========================================================================
ADD_NOISE   = False   # <-- set to True for the noise exercise
NOISE_LEVEL = 0.05    # 0.05 = 5 %

rng = np.random.default_rng(seed=1)
t_noise = np.arange(t_obs, t_end + 1, 1 / 24)            # one value per hour
eta     = 1 + NOISE_LEVEL * rng.standard_normal(t_noise.size)


def measurements(t, V, T):
    dV, _ = patient(t, V, T)
    V_meas = V
    if ADD_NOISE:
        k = int((t - t_obs) * 24)
        V_meas = V * eta[k]
    y1 = V_meas / lam
    y2 = dV / lam
    return y1, y2


# =========================================================================
# STEP 5. State-affine form (Eq. 15) for EXAMPLE 2.             <-- NEW
#
#   dx2/dt = x2^2/x1 + (alpha*delta_T - gamma) x1 - lambda*delta_T*sigma x1^2
#            - delta_T x2 - lambda*sigma x1 x2
#            + r * G * ( x2 - alpha x1 + lambda*sigma x1^2 )
#
#   Group the unknowns:
#     -gamma * x1                              -> phi column 1:  -y1
#     r * G*(x2 - alpha x1 + lambda sigma x1^2) -> phi column 2:  G(y)*(y2 - alpha y1 + lambda sigma y1^2)
#   Known coefficients times x1, x2           -> A(y)
#   Pure functions of y                       -> beta(y)
# =========================================================================
def state_affine_matrices(y1, y2):
    Gy = G(lam * y1)

    a21 = alpha * delta_T - lam * delta_T * sigma * y1     # no r here anymore
    a22 = -delta_T - lam * sigma * y1                      # no r here anymore

    A = np.array([[0.0, 1.0],
                  [a21, a22]])

    beta = np.array([[0.0],
                     [y2**2 / y1]])

    phi = np.array([[0.0, 0.0],
                    [-y1, Gy * (y2 - alpha * y1 + lam * sigma * y1**2)]])  # 2 x 2

    C = np.eye(2)
    return A, beta, phi, C


# =========================================================================
# STEP 6. Observer tuning knobs (Eq. 17)
# =========================================================================
rho_x     = 100.0
rho_theta = 4.0
Q         = np.eye(2)


# =========================================================================
# STEP 7. Patient + adaptive observer in ONE ODE system (Eq. 17)
#
#   Layout of z (16 entries):                                    <-- NEW
#     z[0]  V           z[1]  T
#     z[2]  x1_hat      z[3]  x2_hat
#     z[4]  S11         z[5]  S12         z[6]  S22
#     z[7]  Lambda11    z[8]  Lambda12    z[9]  Lambda21   z[10] Lambda22
#     z[11] Gamma11     z[12] Gamma12     z[13] Gamma22    (Gamma is symmetric)
#     z[14] gamma_hat   z[15] r_hat
# =========================================================================
def patient_and_observer(t, z):
    # --- unpack ----------------------------------------------------------
    V, T      = z[0], z[1]
    x_hat     = np.array([[z[2]], [z[3]]])
    S         = np.array([[z[4], z[5]],
                          [z[5], z[6]]])
    Lambda    = np.array([[z[7], z[8]],
                          [z[9], z[10]]])
    Gamma     = np.array([[z[11], z[12]],
                          [z[12], z[13]]])
    theta_hat = np.array([[z[14]], [z[15]]])

    # --- the patient ------------------------------------------------------
    dV, dT = patient(t, V, T)

    # --- measurements -----------------------------------------------------
    y1, y2 = measurements(t, V, T)
    y = np.array([[y1], [y2]])

    # --- state-affine matrices from measured data -------------------------
    A, beta, phi, C = state_affine_matrices(y1, y2)

    # --- adaptive observer, Eq. (17): IDENTICAL to Example 1 --------------
    S_inv     = np.linalg.inv(S)
    Gamma_inv = np.linalg.inv(Gamma)
    error     = y - C @ x_hat

    dx_hat     = A @ x_hat + beta + phi @ theta_hat \
                 + (S_inv @ C.T + Lambda @ Gamma_inv @ Lambda.T @ C.T) @ Q @ error
    dS         = -rho_x * S - A.T @ S - S @ A + C.T @ Q @ C
    dLambda    = (A - S_inv @ C.T @ Q @ C) @ Lambda + phi
    dGamma     = -rho_theta * Gamma + Lambda.T @ C.T @ Q @ C @ Lambda
    dtheta_hat = Gamma_inv @ Lambda.T @ C.T @ Q @ error

    # --- pack -------------------------------------------------------------
    return [dV, dT,
            dx_hat[0, 0], dx_hat[1, 0],
            dS[0, 0], dS[0, 1], dS[1, 1],
            dLambda[0, 0], dLambda[0, 1], dLambda[1, 0], dLambda[1, 1],
            dGamma[0, 0], dGamma[0, 1], dGamma[1, 1],
            dtheta_hat[0, 0], dtheta_hat[1, 0]]


# =========================================================================
# STEP 8. Initial conditions of the observer
# =========================================================================
y1_0, y2_0 = measurements(t_obs, V_at_0, T_at_0)

x_hat_0     = [0.5 * y1_0, 0.5 * y2_0]   # wrong guess of the states
S_0         = [1.0, 0.0, 1.0]            # S(0) = I
Lambda_0    = [0.0, 0.0, 0.0, 0.0]       # Lambda(0) = 0
Gamma_0     = [1.0, 0.0, 1.0]            # Gamma(0) = I
theta_hat_0 = [0.05, 0.05]               # guesses: gamma (true 0.189), r (true 0.194)

z0 = [V_at_0, T_at_0] + x_hat_0 + S_0 + Lambda_0 + Gamma_0 + theta_hat_0


# =========================================================================
# STEP 9. Solve
# =========================================================================
t_eval = np.linspace(t_obs, t_end, 2001)
sol = solve_ivp(patient_and_observer, [t_obs, t_end], z0, t_eval=t_eval,
                rtol=1e-6, atol=1e-9, max_step=1 / 24)
print(sol.message)

t         = sol.t
V         = sol.y[0]
T         = sol.y[1]
V_hat     = lam * sol.y[2]
dV_hat    = lam * sol.y[3]
gamma_hat = sol.y[14]
r_hat     = sol.y[15]


# =========================================================================
# STEP 10. Physical parameters and the unmeasured state T
# =========================================================================
c_T_hat = gamma_hat / S_T
c_T_hat_safe = np.where(c_T_hat > 0, c_T_hat, np.nan)
T_hat = (alpha - sigma * V_hat - dV_hat / V_hat) / c_T_hat_safe

print(f"True c_T = {c_T:.4e}   final c_T_hat = {c_T_hat[-1]:.4e}")
print(f"True r   = {r:.4f}       final r_hat   = {r_hat[-1]:.4f}")
for day in [2, 4, 6, 8, 10]:
    k = np.argmin(abs(t - day))
    print(f"  day {day:2d}: c_T_hat = {c_T_hat[k]:.3e}, r_hat = {r_hat[k]:.3f}")


# =========================================================================
# STEP 11. Persistent-excitation (PE) meter
#
#   PE condition (Aranovskiy et al. 2017, Eq. 4), written for our observer:
#       lambda_min( integral_t^{t+T}  Lambda' C' Q C Lambda ds ) >= delta > 0
#   Gamma(t) is this integral with exponential forgetting, so
#   lambda_min(Gamma(t)) measures ONLINE how much information we have about
#   the LEAST excited direction of theta = [gamma; r].
#   Reference: pure forgetting e^{-rho_theta t} Gamma(0) = "no new information".
# =========================================================================
lam_min_Gamma = np.zeros_like(t)
for k in range(t.size):
    Gamma_k = np.array([[sol.y[11, k], sol.y[12, k]],
                        [sol.y[12, k], sol.y[13, k]]])
    lam_min_Gamma[k] = np.min(np.linalg.eigvalsh(Gamma_k))  # smallest eigenvalue
no_info = np.exp(-rho_theta * (t - t_obs)) * 1.0               # Gamma(0) = I


# =========================================================================
# STEP 12. Plots
# =========================================================================
fig, ax = plt.subplots(4, 1, figsize=(8, 10), sharex=True)

ax[0].semilogy(t, V, 'b-', lw=2, label='Actual')
ax[0].semilogy(t, V_hat, 'r--', lw=2, label='Observer')
ax[0].set_ylabel('V(t)  [copies/ml]')
ax[0].legend()
ax[0].set_title('Example 2: estimating $c_T$ and $r$ from V and dV/dt')

ax[1].plot(t, T, 'b-', lw=2, label='Actual')
ax[1].plot(t, T_hat, 'r--', lw=2, label='Estimated')
ax[1].set_ylabel('T(t)  [cells]')
ax[1].set_ylim(0, 1.5 * T.max())
ax[1].legend()

ax[2].plot(t, c_T * np.ones_like(t), 'b-', lw=2, label='Actual')
ax[2].plot(t, c_T_hat, 'r--', lw=2, label='Estimated')
ax[2].set_ylabel('$c_T$')
ax[2].legend()

ax[3].plot(t, r * np.ones_like(t), 'b-', lw=2, label='Actual')
ax[3].plot(t, r_hat, 'r--', lw=2, label='Estimated')
ax[3].set_ylabel('$r$')
ax[3].set_xlabel('t [days]')
ax[3].legend()

for a in ax:
    a.grid(True)
plt.tight_layout()
plt.show()

# --- PE meter -------------------------------------------------------------
plt.figure(figsize=(8, 3.5))
plt.semilogy(t, lam_min_Gamma, 'k-', lw=2, label=r'$\lambda_{\min}(\Gamma(t))$')
plt.semilogy(t, no_info, 'k:', lw=1.5, label=r'pure forgetting $e^{-\rho_\theta t}\Gamma(0)$')
plt.xlabel('t [days]')
plt.ylabel(r'$\lambda_{\min}(\Gamma)$')
plt.title('PE meter: far above the dotted line = data are exciting the parameters')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
