function [ rho ] = calTransmissionProb(time_step1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   generate the rate of V change: change of V per day.
%   row number corresponds to #previous infection k
%   column number corresponds to binding avidity from 0:0.01:5;
    global params;
    params.filename = 'params_std_reinfect.t2000';
    params = InitParameters(['dat/' params.filename]);
    
    deltaVMatrix = [];
    v = 0.5;
    n0 = 4;
    params.kc = 20;
    rho = [];
    for k = 0:99
        k
        vini = v;
        vt = getVChange_ode(vini,k,time_step1);
        rho(k+1) = get_rho(k, vt, params.p, params.r, params.a, params.b, 1 , n0);
    end
    disp done;
    save('dat/transmission_meds.mat','rho');
    
function  val = get_rho(k, v, p, r, a, b, c, n0)
   %p, r: parameters for transmission benefit
   %a, b: parameters for replication cost
   %c: contact rate
   P_Ab = exp(-p*(v+1));        % Probability to be recognized by antibodies
   %P_Trans = (1-P_Ab).^(r*k);   % Probability to escape immunity
   P_Trans = (1-P_Ab).^(r*k);   % Probability to escape immunity
   P_Rep = exp(-a*v.^b);        % Probability to replicate
   %val = c.*(P_Trans.*P_Rep);
   R0 = (P_Trans.*P_Rep).*n0;
   rho1 = 1 - R0.^-1;
   loc_r = find(rho1<=0);
   rho1(loc_r) = 0;
   val = rho1;
end
end
   