%Reconstruct phylogenies
function void = main_plot_lineage(void)

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
%proj = '20130827_lar_n4_c08_N5';

%optimum
%proj1 = '20130827_lar_n4_c08_N5';
%proj2 = '20130827_lar_fixb_c05';
global metadata;
proj1 = metadata.ibms.proj1;

%kc = 0.1
%proj1 = '20130828_lar_n4_c08_kc01';
%proj2 = '20130829_lar_fixb_c08_kc01';

%kc = 0.5
%proj1 = '20130829_lar_n4_c08_kc05';
%proj2 = '20130829_lar_fixb_c08_kc05';
lin_barray = {};
lin_varray = {};

figure; hold on;
plot_lineage(proj1,'b');
plot_lineage(proj2,'r');
figure;
hist(cell2mat(lin_varray(1)),20)
figure;
hist(cell2mat(lin_varray(2)),20)
disp 'plot the virus traits on the lineage.';


function [] = plot_lineage(proj, color)
%traitfile = ['dat/' proj '/virus_traits'];
traitfile = [metadata.ibms.out_dir '/virus_traits'];
load(traitfile);
global epi_params;
epi_params.tRange_stoch(1,1)=365*20;
%epi_params.tRange_stoch(1,2)=365*29; %thesis version 365*29
epi_params.tRange_stoch(1,2)=365*45; %final version 365*45 (1968-2013)


count = length(dat_VirusesArray(:,1));
vid = dat_VirusesArray(:,1);
births = dat_VirusesArray(:,2);
deaths = dat_VirusesArray(:,3);
parents = dat_VirusesArray(:,4);
initialV = dat_VirusesArray(:,7);
finalV = dat_VirusesArray(:,8);


n_total = length(parents);
%tip = 27030039;
tip = vid(end);
lin = getParentList(parents, tip)';
lin_birth = births(lin(1:end-1));
loc_s = find(lin_birth>365*10);
lin_e = lin(loc_s);

lin_birth = births(lin_e(1:end-1));
lin_initialV = initialV(lin_e(2:end-1));
lin_initialV(end+1) = 0.5;
plot(lin_birth,lin_initialV, color);
lin_barray(end+1)={lin_birth};
lin_varray(end+1)={lin_initialV};
end
%locs = find(parent == 0); 
%parent(locs) = NaN;

%save(filename_infectionTreeData, 'births', 'deaths', 'parent', 'binding');
end
