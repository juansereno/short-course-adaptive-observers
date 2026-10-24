# MATLAB code

| File | What it does |
|---|---|
| [`Ex1_main_cT.m`](Ex1_main_cT.m) | **Run this.** Example 1: estimates $\gamma = c_T S_T$ and reconstructs the T cells |
| [`Ex1_AdapObs_cT.m`](Ex1_AdapObs_cT.m) | patient + adaptive observer ODEs (Eq. 17) for Example 1 |
| [`Ex2_main_cT_r.m`](Ex2_main_cT_r.m) | **Run this.** Example 2: estimates $\gamma$ and $r$ |
| [`Ex2_AdapObs_cT_r.m`](Ex2_AdapObs_cT_r.m) | patient + adaptive observer ODEs for Example 2 |
| [`patient_model.m`](patient_model.m) | in-host model, Eqs. (1)–(2): the "virtual patient" |

The main scripts follow the same steps as the slides:

1. Parameters of the virtual patient
2. Simulate the days before symptoms (no data)
3. Measurements, with an optional noise switch `param.add_noise`
4. Observer tuning: `rho_x`, `rho_theta`, `Q`
5. Initial conditions of the observer
6. Solve patient + observer together with `ode45`
7. Recover $\hat c_T$ and $\hat T$
8. PE meter: $\lambda_{\min}(\Gamma(t))$
9. Plots

In the `*_AdapObs_*.m` files the state vector `z` is unpacked by index, the matrices
$A(y)$, $\beta(y)$, $\varphi(y)$, $C$ are written out, and Eq. (17) takes five lines.
The only differences between the two examples are $A$, $\varphi$ and the sizes of
$\Lambda$, $\Gamma$ and $\hat\theta$.

Needs MATLAB R2016b or newer. The numerical results were checked in GNU Octave 8.4 and match the Python version.
