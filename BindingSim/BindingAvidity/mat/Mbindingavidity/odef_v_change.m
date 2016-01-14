%user defined ode function
%deltaV/deltaT = kc%dBeta_dV 
function X_dot = odef_v_change(t,x,pars)
global params;
% check whether params is empty 
kc = params.kc;
V = x;
%Beta = get_beta_list([0:params.N_Infect-1]',V, params.p, params.r, params.a, params.b, params.c);
X_dot = kc*dBeta_dV(pars.k,V,params.p, params.r, params.a, params.b, params.c);
end


