%==========================================================================
% Mini-course: Adaptive observers for immune digital twins
% EXAMPLE 1 -- Estimate ONE parameter: the T-cell killing rate c_T
%             (and reconstruct the unmeasured CD8+ T cells)
%
% Based on: Rodriguez-Linan, Sereno, Quiroz-Compean, Hernandez-Vargas,
%           "Adaptive Observers for Calibrating Immune Digital Twins in
%            Viral Infections", IEEE CDC 2025.  Eqs. (1)-(2), (13)-(17).
%
% Files needed in the same folder:
%   Ex1_main_cT.m      (this script)
%   Ex1_AdapObs_cT.m   (patient + adaptive observer ODEs)
%   patient_model.m    (in-host model, Eqs. 1-2)
%
% This code is written to be READ, not to be fast.
%==========================================================================
% Change of coordinates (Eq. 13):  x1 = V/lambda,  x2 = (dV/dt)/lambda
% Measured outputs:                y1 = x1,        y2 = x2   (C = I)
% State-affine form (Eq. 15):      dx/dt = A(y) x + beta(y) + phi(y) theta
% Unknown parameter:               theta = gamma = c_T * S_T
%   (gamma ~ 0.19 is much friendlier for the solver than c_T ~ 2e-6;
%    afterwards c_T = gamma / S_T)
%==========================================================================
clear; close all; clc;

%% STEP 1. Parameters of the virtual patient ("Mean" patient of
%          Hernandez-Vargas & Velasco-Hernandez, Annu. Rev. Control 2020)
param.p      = 8.57;      % virus replication rate               [1/day]
param.K      = 2.35e9;    % carrying capacity of the virus       [copies/ml]
param.c      = 2.4;       % viral clearance rate                 [1/day]
param.cT     = 1.89e-6;   % killing rate of virus by T cells     <-- UNKNOWN
param.deltaT = 0.1;       % death rate of T cells                [1/day]
param.kT     = 1.26e5;    % half-saturation constant of T cells
param.m      = 2;         % Hill coefficient
param.r      = 0.194;     % T-cell activation rate               [1/day]
param.T0     = 1e6;       % initial number of T cells
param.ST     = param.deltaT*param.T0;  % homeostatic replenishment (= 1e5)
param.V0     = 0.31;      % initial viral load at t = -3 days

% Lumped parameters used in the paper
param.alpha  = param.p - param.c;     % net viral growth rate
param.sigma  = param.p/param.K;       % logistic term
gamma_true   = param.cT*param.ST;     % <-- what the observer estimates (theta)

param.lambda = 1e8;                   % scaling factor lambda

%% STEP 2. Before symptoms (t = -3 ... 0 days) there are NO measurements.
%          We only simulate the patient to know the state at t = 0.
t_start = -3;  param.t_obs = 0;  t_end = 21;

[~, x_pre] = ode45(@(t,x) patient_model(t, x, param), [t_start param.t_obs], ...
                   [param.V0; param.T0]);
V_at_0 = x_pre(end,1);
T_at_0 = x_pre(end,2);
fprintf('At t = %g (first measurement): V = %.3e, T = %.3e\n', param.t_obs, V_at_0, T_at_0);

%% STEP 3. Measurements: optional multiplicative noise on V
%          V_meas = V*eta,  eta ~ N(1, noise_level), one value per hour
param.add_noise   = false;   % <-- set to true for the noise exercise
param.noise_level = 0.05;    % 0.05 = 5 %
rng(1);                      % same random numbers every run
param.eta = 1 + param.noise_level*randn(24*(t_end - param.t_obs) + 24, 1);

%% STEP 4. Observer tuning knobs (Eq. 17)
param.rho_x     = 100;       % forgetting factor of S     (state correction speed)
param.rho_theta = 4;         % forgetting factor of Gamma (parameter adaptation speed)
param.Q         = eye(2);    % weight on the output error

%% STEP 5. Initial conditions of the observer (we do NOT know the truth!)
dx_at_0 = patient_model(param.t_obs, [V_at_0; T_at_0], param);
y1_0 = V_at_0/param.lambda;
y2_0 = dx_at_0(1)/param.lambda;

x_hat_0     = [0.5*y1_0; 0.5*y2_0];   % a wrong guess of the states (50 % off)
S_0         = [1; 0; 1];              % S(0) = I      -> S11, S12, S22
Lambda_0    = [0; 0];                 % Lambda(0) = 0
Gamma_0     = 1;                      % Gamma(0) = 1
gamma_hat_0 = 0.05;                   % initial guess of gamma (true value 0.189)

z0 = [V_at_0; T_at_0; x_hat_0; S_0; Lambda_0; Gamma_0; gamma_hat_0];

%% STEP 6. Solve patient + observer together
opts   = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'MaxStep', 1/24);
t_span = linspace(param.t_obs, t_end, 2001);
[t, z] = ode45(@(t,z) Ex1_AdapObs_cT(t, z, param), t_span, z0, opts);

V         = z(:,1);
T         = z(:,2);
V_hat     = param.lambda*z(:,3);     % back to original units
dV_hat    = param.lambda*z(:,4);
gamma_hat = z(:,11);

%% STEP 7. Back to the physical parameter and the unmeasured state T
%   c_T = gamma / S_T                                  (S_T known)
%   From Eq. (1):  c_T*T = alpha - sigma*V - (dV/dt)/V
%            =>    T_hat = (alpha - sigma*V_hat - dV_hat/V_hat) / c_T_hat
cT_hat = gamma_hat/param.ST;
cT_hat_safe = cT_hat;
cT_hat_safe(cT_hat <= 0) = NaN;      % T_hat is undefined if c_T_hat <= 0
T_hat = (param.alpha - param.sigma*V_hat - dV_hat./V_hat)./cT_hat_safe;

fprintf('True  c_T     = %.4e\n', param.cT);
fprintf('Final c_T_hat = %.4e  (relative error %.2f %%)\n', ...
        cT_hat(end), abs(cT_hat(end) - param.cT)/param.cT*100);
for day = [2 4 6 8]
    [~, k] = min(abs(t - day));
    fprintf('  day %2d: c_T_hat = %.3e\n', day, cT_hat(k));
end

%% STEP 8. Persistent-excitation (PE) meter
%   PE condition (Aranovskiy et al. 2017, Eq. 4), written for our observer:
%       lambda_min( integral_t^{t+T} Lambda'*C'*Q*C*Lambda ds ) >= delta > 0
%   Gamma(t) is exactly this integral with exponential forgetting:
%       Gamma(t) = e^{-rho_theta t} Gamma(0) + int_0^t e^{-rho_theta(t-s)} Lambda'*C'*Q*C*Lambda ds
%   so lambda_min(Gamma(t)) tells us ONLINE how much information we have.
%   Here q = 1: Gamma is a scalar and lambda_min(Gamma) = Gamma.
lam_min_Gamma = z(:,10);
no_info = exp(-param.rho_theta*(t - param.t_obs))*Gamma_0;   % Gamma if Lambda = 0

%% STEP 9. Plots
figure('Position', [100 50 700 900]);

subplot(4,1,1)
semilogy(t, V, 'b-', t, V_hat, 'r--', 'LineWidth', 2); grid on
ylabel('V(t) [copies/ml]'); legend('Actual', 'Observer')
title('Example 1: estimating c_T from V and dV/dt')

subplot(4,1,2)
plot(t, T, 'b-', t, T_hat, 'r--', 'LineWidth', 2); grid on
ylabel('T(t) [cells]'); legend('Actual', 'Estimated')
ylim([0 1.5*max(T)])

subplot(4,1,3)
plot(t, param.cT*ones(size(t)), 'b-', t, cT_hat, 'r--', 'LineWidth', 2); grid on
ylabel('c_T'); legend('Actual', 'Estimated')

% Persistent excitation: the regressor phi_21 = -y1 must "move" to learn theta
subplot(4,1,4)
plot(t, -V/param.lambda, 'k-', 'LineWidth', 2); grid on
ylabel('\phi_{21} = -y_1'); xlabel('t [days]')

%% PE meter figure
figure('Position', [100 50 700 300]);
semilogy(t, lam_min_Gamma, 'k-', t, no_info, 'k:', 'LineWidth', 2); grid on
xlabel('t [days]'); ylabel('\lambda_{min}(\Gamma)')
legend('\lambda_{min}(\Gamma(t))', 'pure forgetting e^{-\rho_\theta t}\Gamma(0)')
title('PE meter: far above the dotted line = data are exciting the parameters')
