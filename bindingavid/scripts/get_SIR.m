function Y = get_SIR(path)
% get_SIR Summary of this function goes here
% return the lists of infected number for all simulations.
%   path contrains the folder   directory

totalRuns = 3;
M = csvread( [path 'scenario_SIR_1.csv']);
Y = M(:,2);
for i=2:totalRuns
	Mi = csvread( [path 'scenario_SIR_' num2str(i) '.csv']);
    Y = [Y Mi(:,2)];     
end

end

