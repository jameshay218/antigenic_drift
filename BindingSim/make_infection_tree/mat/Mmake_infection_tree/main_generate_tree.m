<<<<<<< HEAD
function [] = main_generate_tree(infile, smpno, starttime, endtime, display, lang)
%Reconstruct phylogenies
%example: main_generate_tree('voutput_small', '30', '10', '200', '0')
%Hsiang-Yu Yuan
%12/01/2016

p = path;
p = path(p,'lib/');
%clear all; close all;
global epi_params;

if exist('infile','var') 
  if strcmp(infile,'')
      infile = ['dat/mod01/voutput_small'];
  else
  end
end

n_seqs = str2num(smpno);

lg = 'c'; %c: c code; m:matlab; s:simple version in matlab 
epi_params.annotation='annotation'; % the 1st displayed annotation text
epi_params.display=0; % display the matlab figure?
epi_params.savefigure=1; % save the figure?
epi_params.savetree=1; % files to be saved. 0:no files, 1:only nexus tree file, 2:all the tree files
epi_params.lg=lg; % language platform

if exist('lang','var') 
  lg = lang;
end

if exist('display','var') 
  if str2num(display) == 1
    epi_params.display=1;
  end
end

if strcmp(lg,'c')
  epi_params.annotation='AntigenicDrift';
  %infile example: infile = ['dat/mod01/voutput_small']; 
  M = csvread([infile '.csv'],2);
  filename_infectionTreeData = strcat(infile, '_tree_', num2str(n_seqs));
  dat_VirusesArray = M;
  count = length(dat_VirusesArray(:,1));
  births = dat_VirusesArray(:,2);
  deaths = dat_VirusesArray(:,3);
  parent = dat_VirusesArray(:,4);
  infectionK = dat_VirusesArray(:,7);
  binding = dat_VirusesArray(:,5);
  bindingFinal = dat_VirusesArray(:,6); 
  antigenicTot = dat_VirusesArray(:,12);
end

if strcmp(lg,'m')
  %infile example: traitfile = ['dat/' proj '/virus_traits'];
  load(infile);
  filename_infectionTreeData = strcat(infile, '_tree_', num2str(n_seqs));
  count = length(dat_VirusesArray(:,1));
  births = dat_VirusesArray(:,2);
  deaths = dat_VirusesArray(:,3);
  parent = dat_VirusesArray(:,4);
  infectionK = dat_VirusesArray(:,5);
  binding = dat_VirusesArray(:,7);
  bindingFinal = dat_VirusesArray(:,8); 
  antigenicTot = dat_VirusesArray(:,11);
end

if strcmp(lg,'s')
  %[births, deaths, parent] = GetInfectionTree(infile);
  count = length(dat_VirusesArray(:,1));
  births = dat_VirusesArray(:,2);
  deaths = dat_VirusesArray(:,3);
  parent = dat_VirusesArray(:,4);  
end
=======
%Reconstruct phylogenies
function void = main_generate_tree(smpno, proj, starttime, endtime)
p = path;
p = path(p,'lib/');
%clear all; close all;


if exist('proj','var') 
else
proj = '20130827_lar_n4_c08_N5';
end




n_seqs = str2num(smpno);

infile = ['dat/' proj '/dat_x_trans_tmp']; %N=0.5x10^6; 5yrs
traitfile = ['dat/' proj '/virus_traits'];
filename_infectionTreeData = strcat(infile, '_', num2str(n_seqs), '_indiv_infectionTree');

%[births, deaths, parent] = GetInfectionTree(infile);
load(infile);
load(traitfile);
count = length(dat_viruses(:,1));
births = dat_viruses(:,2);
deaths = dat_viruses(:,3);
parent = dat_viruses(:,4);
infectionK = dat_viruses(:,5);
binding = dat_VirusesArray(:,7);
bindingFinal = dat_VirusesArray(:,8); 
>>>>>>> ca2a392d1917ea7e07a21caa046c6ef60b579d07


%parfile = 'dat/params.mat';
%dat = load(parfile)
<<<<<<< HEAD
if exist('starttime','var') 
    starttime = str2num(starttime);
=======
global epi_params;
if exist('starttime','var') 
    starttime = str2num('starttime');
>>>>>>> ca2a392d1917ea7e07a21caa046c6ef60b579d07
else
    starttime = 30;
end
if exist('endtime','var') 
<<<<<<< HEAD
    endtime = str2num(endtime);
=======
    endtime = str2num('endtime');
>>>>>>> ca2a392d1917ea7e07a21caa046c6ef60b579d07
else
    endtime = max(deaths)-10;
end
epi_params.tRange_stoch(1,1)=starttime;
epi_params.tRange_stoch(1,2)=endtime; %final version 365*45 (1968-2013)



<<<<<<< HEAD
locs = find(parent == 0); 
parent(locs) = NaN;

if strcmp(lg,'c')
    save(filename_infectionTreeData, 'births', 'deaths', 'parent', 'infectionK', 'binding', 'bindingFinal','antigenicTot');
elseif strcmp(lg,'m')
    save(filename_infectionTreeData, 'births', 'deaths', 'parent', 'infectionK', 'binding', 'bindingFinal');
end

BuildTree_indiv_nexus(filename_infectionTreeData, n_seqs); % build nexus tree

% delete the temporary files
if exist([filename_infectionTreeData '.mat'], 'file') == 2
    delete([filename_infectionTreeData '.mat']);
    disp(['delete the temporary file: ' filename_infectionTreeData '.mat']);
end

end
=======

locs = find(parent == 0); 
parent(locs) = NaN;

save(filename_infectionTreeData, 'births', 'deaths', 'parent', 'infectionK', 'binding', 'bindingFinal');



%BuildTree_indiv_v2(filename_infectionTreeData, n_seqs);
BuildTree_indiv_nexus(filename_infectionTreeData, n_seqs);
>>>>>>> ca2a392d1917ea7e07a21caa046c6ef60b579d07
