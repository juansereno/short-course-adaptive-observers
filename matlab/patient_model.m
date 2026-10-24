%==========================================================================
% patient_model.m
% In-host SARS-CoV-2 model, Eqs. (1)-(2) of the paper:
%
%   dV/dt = p V (1 - V/K) - c_T V T - c V
%   dT/dt = S_T + r T G(V) - delta_T T,     G(V) = V^m / (V^m + k_T^m)
%
% This is our "virtual patient". In a real application this function is
% replaced by the measurements coming from the clinic.
%==========================================================================
function dxdt = patient_model(t, x, param)

V = x(1);   % viral load
T = x(2);   % CD8+ T cells

G = V^param.m / ( V^param.m + param.kT^param.m );   % T-cell activation

dV = param.p*V*(1 - V/param.K) - param.cT*V*T - param.c*V;
dT = param.ST + param.r*T*G - param.deltaT*T;

dxdt = [dV;
        dT];
end
