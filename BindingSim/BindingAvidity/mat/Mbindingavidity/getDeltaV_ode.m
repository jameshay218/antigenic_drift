%calculate deltaV during unit time_step
function [dV finalV]= getDeltaV_ode(initialV,infectionK,time_step)
dT = time_step;
v = initialV;
pars.k = infectionK;
t0 = 0;
nsteps = 20;
tspan = linspace(t0,dT,nsteps);
[v1] = ode2(@odef_v_change, tspan, v, pars);
%[t1 v1] = ode45(@(t,y)odef_v_change(t,y,pars), tspan, v);
finalV = v1(end,:)';
if finalV<0
  error('V < 0')
end
dV = finalV-v;
end



