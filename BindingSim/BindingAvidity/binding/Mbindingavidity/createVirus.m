%Create a virus particle with birth time, death time and source virus.
%Birth time: the time when individual is infected
%Death time
%   Detailed explanation goes here
%Jul 20, 2012

function [ Virus ] = createVirus( vid, birth, death, sourceVirus, infectionK, beta, initialV, currentV )

global params; 
Virus.vid = vid; 
Virus.birth = birth;
Virus.death = death;
Virus.parent = sourceVirus;
Virus.infectionK = infectionK;

%Sk = [0 params.N_Infect];
%beta = get_beta(Sk, v, params.p, params.r, params.a, params.b, params.c);
Virus.beta = beta;
if exist('initialV')
    Virus.initialV = initialV;
end
if exist('currentV')
    Virus.currentV = currentV;
end
end

