function void = main_simulate_genealogy_indiv(void)

clear all; close all;

infile = 'bindingAvidity_stochSimulation_I_matrix_alpha1';

filename_infectionTreeData = strcat(infile, '_indiv_infectionTree');

%[births, deaths, parent] = GetInfectionTree(infile);
load(infile)

births = I_matrix(:,1);
deaths = I_matrix(:,2);
parent = I_matrix(:,3);
locs = find(parent == 0); 
parent(locs) = NaN;

save(filename_infectionTreeData, 'births', 'deaths', 'parent');

%n_seqs = 400;
n_seqs = 40;

BuildTree_indiv(infile, filename_infectionTreeData, n_seqs);
