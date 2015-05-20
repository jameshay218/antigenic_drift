%Reconstruct phylogenies
function void = main_simulate_genealogy_indiv_binding(void)

%clear all; close all;

%infile = 'dat/dat_x_trans_tmp';
%infile = 'dat/20130225/dat_x_trans_tmp';
%infile = 'dat/20130226_final/dat_x_trans';
%infile = 'dat/20130228/dat_x_trans_tmp';
%infile = 'dat/20130308/dat_x_trans_tmp';
%infile = 'dat/20130310/dat_x_trans_tmp'; %Job:20130308_2;N=5x10^5
%infile = 'dat/20130308_2/20130314/dat_x_trans_tmp'; %Job:20130308_2_2;N=6x10^5  %Thesis version
%infile = 'dat/20130311/out/20130315/dat_x_trans_tmp'; %Job:20130315_2;N=8x10^5 
%infile = 'dat/20130323/dat_x_trans_tmp'; %N=10^6
%infile = 'dat/20130325/dat_x_trans_tmp'; %N=10^6
%infile = 'dat/20130717/dat_x_trans_tmp'; %N=2x10^6
%infile = 'dat/20130718/dat_x_trans_tmp'; %N=3x10^6; 45 yrs
%proj = '20130801';
%proj = '20130809_med';
%proj = '20130812_med'; %nv=16; r=2
%proj = '20130814_med_fixb';
%proj = '20130816_med_n4_c09'; %nv=4; c=0.9
%proj = '20130816_med_n8_c08'; %nv=4; c=0.9
%proj = '20130819_lar_fixb28'; %beta=0.28
%proj = '20130818_lar_n4_c07';  %nv=4; c=0.7
proj = '20130827_lar_n4_c08_N5';
%proj = '20130828_lar_n4_c08_kc01';

%n_seqs = 400;
%n_seqs = 300;
n_seqs = 3;

infile = ['dat/' proj '/dat_x_trans_tmp']; %N=0.5x10^6; 5yrs
traitfile = ['dat/' proj '/virus_traits'];
filename_infectionTreeData = strcat(infile, '_', num2str(n_seqs), '_indiv_infectionTree');

%[births, deaths, parent] = GetInfectionTree(infile);
load(infile);
load(traitfile);


%parfile = 'dat/params.mat';
%dat = load(parfile)
global epi_params;
%epi_params = dat.epi_params;
%epi_params.tRange_min=1;
%epi_params.tRange_stoch(1,1)=1;
epi_params.tRange_stoch(1,1)=365*20;
%epi_params.tRange_stoch(1,2)=365*29; %thesis version 365*29
epi_params.tRange_stoch(1,2)=365*45; %final version 365*45 (1968-2013)
%births = I_matrix(:,1);
%deaths = I_matrix(:,2);
%parent = I_matrix(:,3);

count = length(dat_viruses(:,1));
births = dat_viruses(:,2);
deaths = dat_viruses(:,3);
parent = dat_viruses(:,4);
binding = dat_VirusesArray(:,8);
%births = dat_viruses(201:50000,2);
%deaths = dat_viruses(201:50000,3);
%parent = dat_viruses(201:50000,4)-200;

locs = find(parent == 0); 
parent(locs) = NaN;

save(filename_infectionTreeData, 'births', 'deaths', 'parent', 'binding');



%BuildTree_indiv_v2(filename_infectionTreeData, n_seqs);
BuildTree_indiv_nexus(filename_infectionTreeData, n_seqs);
