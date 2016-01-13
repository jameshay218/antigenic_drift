function [ VChangeMatrix ] = genVChangeMatrix(time_step1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   generate the rate of V change: change of V per day.
%   row number corresponds to #previous infection k
%   column number corresponds to binding avidity from 0:0.01:5;
    global params;
    params.filename = 'params_std_reinfect.t2000';
    params = InitParameters(['dat/' params.filename]);
    
    deltaVMatrix = [];
    v = 0:0.01:5;
    for k = 1:100
        k-1
       for i=1:length(v)    
            vini = v(i);
            deltaV = getVChange_ode(vini,k-1,time_step1)-vini;
            deltaVMatrix(k,i) = deltaV;
       end
    end
    save('dat/deltaVMatrix.mat','deltaVMatrix');
end

