function [finalV dV]= getVChange_ode(initialV,infectionK,time_step)
global params;

%if isempty(params)
%    params.filename = 'params_std_reinfect.t2000';
%    params = InitParameters(['dat/' params.filename]);
%end

dT = time_step;
v = initialV;
pars.k = infectionK;
t0 = 0;
nsteps = 30;
tspan = linspace(t0,dT,nsteps);
%[v1] = ode2('odef_v_change',tspan,v,k); %If doesn't work, use ode15s
[v1] = ode2(@odef_v_change, tspan, v, pars);
%[t1, v1] = ode45('odef_v_change',t,X0'); %If doesn't work, use ode15s
%[t1, v1] = ode45(@odef_v_change,tspan,v); %If doesn't work, use ode15s
finalV = v1(end,:);
if finalV<0
  error('V < 0')
end
dV = finalV-v;
end



