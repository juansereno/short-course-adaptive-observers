%==========================================================================
% Mini-course: Adaptive observers for immune digital twins
% EXAMPLE 2 -- Estimate TWO parameters: the T-cell killing rate c_T and
%             the T-cell activation rate r (and reconstruct T cells)
%
% Based on: Rodriguez-Linan, Sereno, Quiroz-Compean, Hernandez-Vargas,
%           "Adaptive Observers for Calibrating Immune Digital Twins in
%            Viral Infections", IEEE CDC 2025.  Eqs. (1)-(2), (13)-(17).
%
% Files needed in the same folder:
%   Ex2_main_cT_r.m      (this script)
%   Ex2_AdapObs_cT_r.m   (patient + adaptive observer ODEs)
%   patient_model.m      (in-host model, Eqs. 1-2)
%
% Same structure as Example 1. What changes (look for  <-- NEW):
%   * theta = [gamma; r]  with gamma = c_T * S_T   (q = 2 unknowns)
%   * the r-terms move from A(y) into phi(y)
%   * Lambda becomes 2x2, Gamma becomes 2x2, theta_hat has 2 entries
%==========================================================================
clear; close all; clc;

%% STEP 1. Parameters of the virtual patient ("Mean" patient)
param.p      = 8.57;      % virus replication rate               [1/day]
param.K      = 2.35e9;    % carrying capacity of the virus       [copies/ml]
param.c      = 2.4;       % viral clearance rate                 [1/day]
param.cT     = 1.89e-6;   % killing rate of virus by T cells     <-- UNKNOWN
param.deltaT = 0.1;       % death rate of T cells                [1/day]
param.kT     = 1.26e5;    % half-saturation constant of T cells
param.m      = 2;         % Hill coefficient
param.r      = 0.194;     % T-cell activation rate               <-- UNKNOWN (NEW)
param.T0     = 1e6;       % initial number of T cells
param.ST     = param.deltaT*param.T0;  % homeostatic replenishment (= 1e5)
param.V0     = 0.31;      % initial viral load at t = -3 days

param.alpha  = param.p - param.c;     % net viral growth rate
param.sigma  = param.p/param.K;       % logistic term
param.lambda = 1e8;                   % scaling factor lambda
% NOTE: param.cT and param.r are used ONLY by the virtual patient
%       (patient_model.m). The observer in Ex2_AdapObs_cT_r.m never uses them.

%% STEP 2. Before symptoms (t = -3 ... 0 days): no measurements
t_start = -3;  param.t_obs = 0;  t_end = 21;

[~, x_pre] = ode45(@(t,x) patient_model(t, x, param), [t_start param.t_obs], ...
                   [param.V0; param.T0]);
V_at_0 = x_pre(end,1);
T_at_0 = x_pre(end,2);
fprintf('At t = %g (first measurement): V = %.3e, T = %.3e\n', param.t_obs, V_at_0, T_at_0);

%% STEP 3. Measurements: optional multiplicative noise on V
param.add_noise   = false;   % <-- set to true for the noise exercise
param.noise_level = 0.05;    % 0.05 = 5 %
rng(1);
param.eta = 1 + param.noise_level*randn(24*(t_end - param.t_obs) + 24, 1);

%% STEP 4. Observer tuning knobs (Eq. 17)
param.rho_x     = 100;
param.rho_theta = 4;
param.Q         = eye(2);

%% STEP 5. Initial conditions of the observer
dx_at_0 = patient_model(param.t_obs, [V_at_0; T_at_0], param);
y1_0 = V_at_0/param.lambda;
y2_0 = dx_at_0(1)/param.lambda;

x_hat_0     = [0.5*y1_0; 0.5*y2_0];   % wrong guess of the states
S_0         = [1; 0; 1];              % S(0) = I
Lambda_0    = [0; 0; 0; 0];           % Lambda(0) = 0   (2x2)       <-- NEW
Gamma_0     = [1; 0; 1];              % Gamma(0) = I    (2x2)       <-- NEW
theta_hat_0 = [0.05; 0.05];           % guesses: gamma (true 0.189), r (true 0.194)

z0 = [V_at_0; T_at_0; x_hat_0; S_0; Lambda_0; Gamma_0; theta_hat_0];

%% STEP 6. Solve patient + observer together
opts   = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'MaxStep', 1/24);
t_span = linspace(param.t_obs, t_end, 2001);
[t, z] = ode45(@(t,z) Ex2_AdapObs_cT_r(t, z, param), t_span, z0, opts);

V         = z(:,1);
T         = z(:,2);
V_hat     = param.lambda*z(:,3);
dV_hat    = param.lambda*z(:,4);
gamma_hat = z(:,15);
r_hat     = z(:,16);

%% STEP 7. Physical parameters and the unmeasured state T
cT_hat = gamma_hat/param.ST;
cT_hat_safe = cT_hat;
cT_hat_safe(cT_hat <= 0) = NaN;      % T_hat is undefined if c_T_hat <= 0
T_hat = (param.alpha - param.sigma*V_hat - dV_hat./V_hat)./cT_hat_safe;

fprintf('True c_T = %.4e   final c_T_hat = %.4e\n', param.cT, cT_hat(end));
fprintf('True r   = %.4f       final r_hat   = %.4f\n', param.r, r_hat(end));
for day = [2 4 6 8 10]
    [~, k] = min(abs(t - day));
    fprintf('  day %2d: c_T_hat = %.3e, r_hat = %.3f\n', day, cT_hat(k), r_hat(k));
end

%% STEP 8. Persistent-excitation (PE) meter
%   PE condition (Aranovskiy et al. 2017, Eq. 4), written for our observer:
%       lambda_min( integral_t^{t+T} Lambda'*C'*Q*C*Lambda ds ) >= delta > 0
%   Gamma(t) is this integral with exponential forgetting, so
%   lambda_min(Gamma(t)) measures ONLINE the information about the LEAST
%   excited direction of theta = [gamma; r].
lam_min_Gamma = zeros(size(t));
for k = 1:length(t)
    Gamma_k = [z(k,12), z(k,13);
               z(k,13), z(k,14)];
    lam_min_Gamma(k) = min(eig(Gamma_k));                    % smallest eigenvalue
end
no_info = exp(-param.rho_theta*(t - param.t_obs))*1;         % Gamma(0) = I

%% STEP 9. Plots
figure('Position', [100 50 700 900]);

subplot(4,1,1)
semilogy(t, V, 'b-', t, V_hat, 'r--', 'LineWidth', 2); grid on
ylabel('V(t) [copies/ml]'); legend('Actual', 'Observer')
title('Example 2: estimating c_T and r from V and dV/dt')

subplot(4,1,2)
plot(t, T, 'b-', t, T_hat, 'r--', 'LineWidth', 2); grid on
ylabel('T(t) [cells]'); legend('Actual', 'Estimated')
ylim([0 1.5*max(T)])

subplot(4,1,3)
plot(t, param.cT*ones(size(t)), 'b-', t, cT_hat, 'r--', 'LineWidth', 2); grid on
ylabel('c_T'); legend('Actual', 'Estimated')

subplot(4,1,4)
plot(t, param.r*ones(size(t)), 'b-', t, r_hat, 'r--', 'LineWidth', 2); grid on
ylabel('r'); xlabel('t [days]'); legend('Actual', 'Estimated')

%% PE meter figure
figure('Position', [100 50 700 300]);
semilogy(t, lam_min_Gamma, 'k-', t, no_info, 'k:', 'LineWidth', 2); grid on
xlabel('t [days]'); ylabel('\lambda_{min}(\Gamma)')
legend('\lambda_{min}(\Gamma(t))', 'pure forgetting e^{-\rho_\theta t}\Gamma(0)')
title('PE meter: far above the dotted line = data are exciting the parameters')
