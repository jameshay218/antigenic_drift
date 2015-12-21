%Reconstruct phylogenies
function void = main_generate_tree_withtraits(filename, smpno, proj, starttime, endtime)
p = path;
p = path(p,'lib/');
%clear all; close all;


if exist('proj','var') 
else
proj = '20130827_lar_n4_c08_N5';
end




n_seqs = str2num(smpno);

%read from matlab
infile = ['dat/' proj '/' filename];
filename_infectionTreeData = strcat(infile, '_', num2str(n_seqs), '_indiv_infectionTree');

%%[births, deaths, parent] = GetInfectionTree(infile);
%load(infile);
%count = length(dat_viruses(:,1));
%births = dat_viruses(:,2);
%deaths = dat_viruses(:,3);
%parent = dat_viruses(:,4);
%infectionK = dat_viruses(:,5);
%binding = dat_VirusesArray(:,7);
%bindingFinal = dat_VirusesArray(:,8); 

%read from csv file
M = csvread([infile '.csv'],2);
dat_viruses = M;
vid	= dat_viruses(:,1);
births = dat_viruses(:,2);
deaths = dat_viruses(:,3);
parent = dat_viruses(:,4);
infectionK = dat_viruses(:,9);
binding = dat_viruses(:,5);
bindingFinal = dat_viruses(:,6);




%Declare global parameters
global epi_params;
if exist('starttime','var') 
    starttime = str2num(starttime);
else
    starttime = 30;
end
if exist('endtime','var') 
    endtime = str2num(endtime);
else
    endtime = max(deaths)-10;
end
epi_params.tRange_stoch(1,1)=starttime;
epi_params.tRange_stoch(1,2)=endtime; %final version 365*45 (1968-2013)




locs = find(parent == 0); 
parent(locs) = NaN;

save(filename_infectionTreeData, 'births', 'deaths', 'parent', 'infectionK', 'binding', 'bindingFinal');



%BuildTree_indiv_v2(filename_infectionTreeData, n_seqs);
BuildTree_indiv_nexus(filename_infectionTreeData, n_seqs);
