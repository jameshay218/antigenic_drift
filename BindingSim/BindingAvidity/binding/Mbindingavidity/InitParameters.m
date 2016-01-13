%Read and initialize Parameters
function [Params] = InitParameters(filename)
pa = readparams(filename);
time_step = 1;
idx = 1;
deltaVMatrix = [];
vlist = 0:0.1:6;

Birth = pa.N*ones(1,pa.N_Strains);
for n=1:pa.N_Infect-1
Birth = [Birth; zeros(1,pa.N_Strains)];
%  for i=1:length(vlist)
%      vini=vlist(i);
%      deltaVMatrix(n,i)=getVChange_ode(vini,n,time_step);
%  end
end
si_bar = [];
Params = pa;
%Params.deltaVMatrix = deltaVMatrix;
%Params = struct('gamma',pa.gamma,'mu',pa.mu,'size',pa.N,'time_step',time_step,'mut',pa.mutation,'b_mut',pa.b_mutation,'N',pa.N,'Birth',Birth,'N_Strains',pa.N_Strains,'N_Infect',pa.N_Infect,'N_Binding',pa.N_Binding,'wan',pa.wan,'Vg',pa.Vg,'si_bar',si_bar,'filename',filename,'p',pa.p,'r',pa.r,'v0',pa.v0,'a',pa.a,'b',pa.b,'c',pa.c);
Params.dBdV = zeros(1,pa.N_Infect+1); %[t,v1,v2,...vn]


end
