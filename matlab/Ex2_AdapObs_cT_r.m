%==========================================================================
% Ex2_AdapObs_cT_r.m
% Patient + adaptive observer (Eq. 17) for EXAMPLE 2:
% unknown theta = [gamma; r]   (gamma = c_T * S_T)
%
% Layout of the big state vector z (16 entries):                  <-- NEW
%   z(1)  V           z(2)  T
%   z(3)  x1_hat      z(4)  x2_hat
%   z(5)  S11         z(6)  S12         z(7)  S22
%   z(8)  Lambda11    z(9)  Lambda12    z(10) Lambda21   z(11) Lambda22
%   z(12) Gamma11     z(13) Gamma12     z(14) Gamma22    (Gamma is symmetric)
%   z(15) gamma_hat   z(16) r_hat
%==========================================================================
function dzdt = Ex2_AdapObs_cT_r(t, z, param)

% ---- 1. Unpack ----------------------------------------------------------
V = z(1);
T = z(2);
x_hat     = [z(3);
             z(4)];
S         = [z(5), z(6);
             z(6), z(7)];
Lambda    = [z(8),  z(9);
             z(10), z(11)];
Gamma     = [z(12), z(13);
             z(13), z(14)];
theta_hat = [z(15);
             z(16)];

% ---- 2. The patient evolves ---------------------------------------------
dpatient = patient_model(t, [V; T], param);
dV = dpatient(1);
dT = dpatient(2);

% ---- 3. Measurements y1 = V/lambda, y2 = (dV/dt)/lambda (+ optional noise)
V_meas = V;
if param.add_noise
    k = floor((t - param.t_obs)*24) + 1;
    V_meas = V*param.eta(k);
end
y1 = V_meas/param.lambda;
y2 = dV/param.lambda;
y  = [y1;
      y2];

% ---- 4. State-affine form (Eq. 15) for EXAMPLE 2                 <-- NEW
%   dx2/dt = x2^2/x1 + (alpha*deltaT - gamma) x1 - lambda*deltaT*sigma x1^2
%            - deltaT x2 - lambda*sigma x1 x2
%            + r * G * ( x2 - alpha x1 + lambda*sigma x1^2 )
%   -gamma*x1                              -> phi column 1:  -y1
%   r*G*(x2 - alpha x1 + lambda sigma x1^2) -> phi column 2:  G(y)*(y2 - alpha y1 + lambda sigma y1^2)
Gy = y1^param.m / ( y1^param.m + (param.kT/param.lambda)^param.m );  % G(lambda*y1)

a21 = param.alpha*param.deltaT - param.lambda*param.deltaT*param.sigma*y1;  % no r here anymore
a22 = -param.deltaT - param.lambda*param.sigma*y1;                          % no r here anymore

A    = [ 0,   1;
         a21, a22];
beta = [ 0;
         y2^2/y1];
phi  = [ 0,   0;
        -y1,  Gy*(y2 - param.alpha*y1 + param.lambda*param.sigma*y1^2)];   % 2 x 2
C    = [ 1, 0;
         0, 1];

% ---- 5. Adaptive observer, Eq. (17): IDENTICAL to Example 1 --------------
Q         = param.Q;
rho_x     = param.rho_x;
rho_theta = param.rho_theta;
y_err     = y - C*x_hat;                                % output estimation error

dx_hat     = A*x_hat + beta + phi*theta_hat ...
             + ( S\C' + Lambda*(Gamma\Lambda')*C' )*Q*y_err;
dS         = -rho_x*S - A'*S - S*A + C'*Q*C;
dLambda    = ( A - S\C'*Q*C )*Lambda + phi;
dGamma     = -rho_theta*Gamma + Lambda'*C'*Q*C*Lambda;
dtheta_hat = Gamma\Lambda'*C'*Q*y_err;

% ---- 6. Pack (same order as z) -------------------------------------------
dzdt = [ dV;
         dT;
         dx_hat(1);
         dx_hat(2);
         dS(1,1);
         dS(1,2);
         dS(2,2);
         dLambda(1,1);
         dLambda(1,2);
         dLambda(2,1);
         dLambda(2,2);
         dGamma(1,1);
         dGamma(1,2);
         dGamma(2,2);
         dtheta_hat(1);
         dtheta_hat(2) ];
end
