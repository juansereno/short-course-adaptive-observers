# Exercises and open projects

All exercises start from the two examples in [`matlab/`](../matlab) or [`python/`](../python).
Before running each change, **write down what you expect to happen**. Then run it and
explain the difference.

## Hands-on exercises

| # | Task | What to change | What to look at |
|:-:|---|---|---|
| 1 | **Run** Example 1 and Example 2. Change the initial guesses $\hat\theta(0)$ and $\Gamma(0)$. | `gamma_hat_0` / `theta_hat_0`, `Gamma_0` (STEP 5 in MATLAB, STEP 8 in Python) | When does $\hat c_T$ start to move? Does a larger $\Gamma(0)$ speed it up or slow it down? |
| 2 | **Other patients.** Use the parameter sets A–I of Hernandez-Vargas & Velasco-Hernandez (2020), table below. Which converge slowest, and why? | `p`, `k_T`, `r`, `c_T` in STEP 1 | Plot $\varphi$ and the PE meter $\lambda_{\min}(\Gamma)$. |
| 3 | **Noise.** Turn on 5 % multiplicative noise on $V$, then re-tune the forgetting factors. | `add_noise = true`, then `rho_x`, `rho_theta` | Try $\rho_x = 10$, $\rho_\theta = 1$. What do you gain, and what do you lose? |
| 4 | **Start late.** Start the observer at day 12 instead of day 0. Is there still persistent excitation? | `t_obs` | The PE meter: does $\lambda_{\min}(\Gamma)$ ever rise above the dotted line? |
| 5 | **Your turn.** Estimate $\theta = [\gamma,\ \delta_T]$ instead of $[\gamma,\ r]$. Find $A(y)$, $\beta(y)$ and $\varphi(y)$, then code them. | a copy of Example 2 | Hint: the $\delta_T$-column of $\varphi$ is $\alpha y_1 - \lambda\sigma y_1^2 - y_2$. Remember that $S_T = \delta_T T_0$. |

### Patient parameter sets (Hernandez-Vargas & Velasco-Hernandez, 2020)

All patients share $K = 2.35\times10^9$, $c = 2.4$, $\delta_T = 0.1$, $m = 2$, $T_0 = 10^6$, $S_T = 10^5$ and $V_0 = 0.31$.

| Patient | $p$ | $k_T$ | $r$ | $c_T$ |
|:-:|---:|---:|---:|---:|
| A | 6.31 | 7.94e7 | 0.794 | 1.58e-6 |
| B | 12.58 | 1.99e6 | 0.126 | 7.94e-6 |
| C | 19.95 | 1.58e3 | 0.020 | 1.58e-5 |
| D | 6.31 | 3.16e4 | 0.251 | 1.58e-6 |
| E | 5.01 | 5.01e5 | 0.316 | 1.00e-6 |
| F | 6.31 | 1.00e7 | 0.398 | 1.26e-6 |
| G | 6.31 | 7.94e4 | 0.158 | 1.99e-6 |
| H | 31.62 | 2.52e3 | 0.050 | 2.52e-5 |
| I | 3.98 | 1.00e3 | 0.794 | 1.58e-8 |
| Mean (used in the examples) | 8.57 | 1.26e5 | 0.194 | 1.89e-6 |

With exact data the observer converges for all nine patients; some take much longer than others.

## Open projects: PAE · undergrad thesis · master thesis

### 1 · Circadian insulin sensitivity in type 1 diabetes

**Objective.** Design an adaptive observer that tracks the time-varying insulin sensitivity
$S_I(t)$ of a person with T1D from continuous glucose monitoring (CGM) and insulin-pump
data. The target is the breakfast/lunch/dinner pattern of Fig. 1 in Visentin et al. (2015).

**Methodology.**
1. Glucose–insulin minimal model with $S_I$ as the unknown parameter.
2. Rewrite it in state-affine form; check observability, identifiability and PE (meals as excitation).
3. Generate in silico data with a circadian $S_I$ profile.
4. Tune $\rho_\theta$ to follow $S_I$ changes; test noise and 5-min CGM sampling.

### 2 · Menstrual cycle monitoring in women with diabetes

**Objective.** Estimate the cycle-dependent insulin sensitivity $S_I(E_2, P_4)$ along the
menstrual cycle of women with diabetes, based on the model of Ramírez Mazo et al. (2026).

**Methodology.**
1. Adapt the glucose–insulin part of the model to diabetes (exogenous insulin input).
2. Choose the unknown parameters of $S_I(E_2, P_4)$ and write the state-affine form.
3. Analyse observability, identifiability and PE over a 28-day cycle.
4. Validate the observer in silico with $E_2$, $P_4$ profiles and CGM data.

### References for the projects

- R. Visentin, C. Dalla Man, Y. C. Kudva, A. Basu, C. Cobelli. *Circadian variability of insulin sensitivity: physiological input for in silico artificial pancreas.* Diabetes Technology & Therapeutics 17(1):1–7, 2015.
- C. Ramírez Mazo, L. Gómez-Echavarría, L. Lema-Pérez, S. Röblitz. *A physiologically based mathematical model describing estradiol and progesterone actions on glucose-insulin dynamics throughout the menstrual cycle.* Journal of Theoretical Biology 631:112518, 2026.
