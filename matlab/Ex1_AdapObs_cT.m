%==========================================================================
% Ex1_AdapObs_cT.m
% Patient + adaptive observer (Eq. 17) for EXAMPLE 1: unknown theta = gamma
% (gamma = c_T * S_T).
%
% Layout of the big state vector z (11 entries):
%   z(1)  V          z(2)  T             <- virtual patient (unknown to observer)
%   z(3)  x1_hat     z(4)  x2_hat        <- estimated states
%   z(5)  S11        z(6)  S12    z(7) S22   <- symmetric matrix S (2x2)
%   z(8)  Lambda1    z(9)  Lambda2       <- Lambda (2x1)
%   z(10) Gamma                          <- Gamma  (1x1)
%   z(11) gamma_hat                      <- theta_hat (1x1)
%==========================================================================
function dzdt = Ex1_AdapObs_cT(t, z, param)

% ---- 1. Unpack ----------------------------------------------------------
V = z(1);
T = z(2);
x_hat     = [z(3);
             z(4)];
S         = [z(5), z(6);
             z(6), z(7)];
Lambda    = [z(8);
             z(9)];
Gamma     = z(10);
theta_hat = z(11);

% ---- 2. The patient evolves (in real life: the real person) -------------
dpatient = patient_model(t, [V; T], param);
dV = dpatient(1);
dT = dpatient(2);

% ---- 3. The clinic measures y1 = V/lambda and y2 = (dV/dt)/lambda --------
V_meas = V;
if param.add_noise
    k = floor((t - param.t_obs)*24) + 1;   % index of the current hour
    V_meas = V*param.eta(k);              % multiplicative noise on V
end
y1 = V_meas/param.lambda;
y2 = dV/param.lambda;
y  = [y1;
      y2];

% ---- 4. State-affine form (Eq. 15), built ONLY from measured data -------
%   dx2/dt = x2^2/x1 + (alpha*deltaT - gamma) x1 - lambda*deltaT*sigma x1^2
%            - r*alpha*x1*G + lambda*r*sigma*x1^2*G
%            - deltaT x2 - lambda*sigma x1 x2 + r x2 G
%   Recipe: unknown gamma -> phi ; known coeff(y)*x -> A ; only-y terms -> beta
Gy = y1^param.m / ( y1^param.m + (param.kT/param.lambda)^param.m );  % G(lambda*y1)

a21 = param.alpha*(param.deltaT - param.r*Gy) ...
      + param.lambda*param.sigma*y1*(param.r*Gy - param.deltaT);
a22 = param.r*Gy - param.lambda*param.sigma*y1 - param.deltaT;

A    = [ 0,   1;
         a21, a22];
beta = [ 0;
         y2^2/y1];
phi  = [ 0;
        -y1];          % theta = gamma enters as  -gamma*x1
C    = [ 1, 0;
         0, 1];        % both x1 and x2 are measured

% ---- 5. Adaptive observer, Eq. (17) --------------------------------------
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
         dLambda(1);
         dLambda(2);
         dGamma;
         dtheta_hat ];
end
