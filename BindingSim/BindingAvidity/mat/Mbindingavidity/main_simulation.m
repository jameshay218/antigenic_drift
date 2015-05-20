function [] = main_simulation(cont)
%contact = num2str(cont);
contact = '15';
restoredefaultpath;
% define path
maindir = [pwd '/'];
locs = regexp(maindir,'main');
locs = regexp(maindir(1:locs),'[\\\/]');
rootdir = maindir(1:locs(end));
odedir = [rootdir 'binding_ode/'];
ibmsdir = [rootdir 'Mbindingavidity/'];
mainoutdir = [rootdir 'mainout/'];

clear x1;
global metadata;
kc = 0.09; %remove contact c from fitness function
ks = 0.0769;

%% Set metadata
%% Define project name
metadata.ode.proj = ['std_adapt_n4_c08']; 
metadata.ibms.proj = ['std_adapt_n4_c' contact];

%% Set metadata for ODE
metadata.ode.filename = 'params_std_reinfect_c08_meds.t2000';
metadata.ode.rho = 'transmission_meds.mat';

% Check the existence of the parameter files
if exist([odedir 'dat/' metadata.ode.filename])
  disp 'ODE parameters are read.'
end
if exist([odedir 'dat/' metadata.ode.rho])
  disp 'ODE transmission array are read.'
end

%% Define output folder
datepath = [datestr(now,10) datestr(now,5) datestr(now,7)];
out_dir = [mainoutdir 'ode/' datepath '/' metadata.ode.proj];
%out_dir = proj;
if(exist(out_dir)==7)
else
   mkdir(out_dir)
end
metadata.ode.out_dir = out_dir;

%%% Run ODE
cd(odedir);
if exist ([metadata.ode.out_dir '/dat_x_eq_ini.mat'])
  disp ('Output already exists. Skip ODE calculating.')
else
  disp ('Run ODE.'); 
  [x1] = run_without_vac(kc, ks, metadata);
end
cd(maindir);


%% Set metadata for IBMS
metadata.ibms.filename = ['params_std_reinfect_c' contact '_lar.dat'];
%%%^^^^^^^^^^^^^
%metadata.ibms.initFileODE = 'dat/20130826/20130827_lar_n4_c08_N5/dat_x_eq_ini.mat';
%metadata.ibms.initFileVirus = 'out/dscr/20130829_lar_n4_c08_kc05/virus_traits.mat';
%metadata.ibms.initFileSIR = 'out/dscr/20130829_lar_n4_c08_kc05/DataTLSIR.mat';
metadata.ibms.initFileODE = [metadata.ode.out_dir '/dat_x_eq_ini.mat'];
metadata.ibms.initFileVirus = '';
metadata.ibms.initFileSIR = '';
metadata.ibms.initFlag = 0; 
metadata.ibms.deltaV = 'dat/deltaVMatrix_kc05.mat';
metadata.ibms.rho = 'dat/transmission_meds.mat';

%% Check parameters
if exist([ibmsdir 'dat/' metadata.ibms.filename])
  disp 'IBMS parameters are read.'
else
  disp 'Error: failed to read IBMS parameters.';
end
if exist(metadata.ibms.initFileVirus) 
  if exist(metadata.ibms.initFileSIR) 
    metadata.ibms.initFlag = 2; %IBM
  end
else
  if exist(metadata.ibms.initFileODE)   
    metadata.ibms.initFlag = 1; %ODE
  end 
end

out_dir = [mainoutdir 'ibms/' datepath '/' metadata.ibms.proj];
%out_dir = proj;
if(exist(out_dir)==7)
else
   mkdir(out_dir)
end
metadata.ibms.out_dir = out_dir;


%%% Run IBM simulation
%% Set the parameter file
cd(ibmsdir);
for i=1:3
metadata.ibms.initFileVirus = [metadata.ibms.out_dir '/virus_traits.mat'];
if exist(metadata.ibms.initFileVirus)
  metadata.ibms.initFileVirus = [metadata.ibms.out_dir '/virus_traits.mat'];
  metadata.ibms.initFileSIR = [metadata.ibms.out_dir '/DataTLSIR.mat'];
  disp ('virus initFile exists, run simulation from previous result.');
  out_dir = [metadata.ibms.out_dir '_eq'];
  if(exist(out_dir)==7) 
  else
   mkdir(out_dir)
  end
  metadata.ibms.out_dir = out_dir;
  disp (['output directory:' metadata.ibms.out_dir]);
  metadata.ibms.initFlag = 2; %IBM
  run_ibms_simulation(metadata);  
  disp (['done! write output to ' metadata.ibms.out_dir]);
else
  disp ('run simulation with ODE as initial condition.');
  metadata.ibms.initFlag = 1; %ODE
  disp (['output directory:' metadata.ibms.out_dir]);
  run_ibms_simulation(metadata);
  disp (['done! write output to ' metadata.ibms.out_dir]);
end
end
cd(maindir);
end

function [x1] = run_without_vac(kc, ks, metadata)
        disp 'run ODE into equilibrium';
        firstrn = 100;        
        x1 = main_sirv_reinf_wan_seq(metadata.ode.filename,metadata.ode.proj,firstrn,ks,kc);
        metadata.ode.out_dir
        load([metadata.ode.out_dir '/dat_x_eq.tmp.mat']);
        x_eq1 = x_eq;
        x_eq = [];
        x_eq(1:100) = x_eq1(1:100);
        x_eq(101:200) = x_eq1(301:400);
        x_eq(201:300) = x_eq1(601:700);
        save([metadata.ode.out_dir '/dat_x_eq_ini.mat'], 'x_eq');
end

function [] = run_ibms_simulation(metadata)
	%% read Ntot from parameter file
	['dat/' metadata.ibms.filename];
	fid = fopen(['dat/' metadata.ibms.filename])
	p = textscan(fid,'%s %s');
	pa_array_headers = p{1};
	pa_array_values = p{2};
	N_tot = str2num(cell2mat(pa_array_values(strcmp(pa_array_headers,'N'))));
    disp (['total population size:' num2str(N_tot)]);
    %% Setup initial population from IBMS
     initViruses = [];
if metadata.ibms.initFlag == 2
    dat1 = open(metadata.ibms.initFileVirus);
    v = dat1.dat_VirusesArray;
    initViruses = v(find(v(:,3)==0),:);
    initViruses(:,1) = [1:length(initViruses(:,1))]; %update vid
    initViruses(:,2) = zeros(length(initViruses(:,1)),1); %update birth date
    initViruses(:,4) = zeros(length(initViruses(:,1)),1); %update parent
    dat2 = open(metadata.ibms.initFileSIR);
    x1 = dat2.dat_sir(end,2:end);
    S = x1(1:100);
    I = x1(101:200);
    R = x1(201:300);
    DataX = [S I R]; % Initial population for S I and R
    if sum(I)~= initViruses(end,1)
      disp 'error: virus number not consistent';
    end
    clear dat1;
    clear dat2;
    clear v;
    period = 365*45;
 end

%% Setup initial population from ODE
if metadata.ibms.initFlag == 1
S = [];%S0->S21
I = [];%I1->I100
init_I = 1000; %Inital infected individuals 
dat = open(metadata.ibms.initFileODE);
x1 = dat.x_eq(1,1:end);%retrieve the last row
ratio = N_tot./sum(x1);
S = round(x1(1:100).*ratio);
I = [init_I init_I init_I init_I init_I zeros(1,95)]; %Add Iini to I2 
%I = x1(101:200);
R = round(x1(201:300).*ratio);
S(1) = S(1)+N_tot-(sum(S)+sum(R))-init_I;
DataX = [S I R]; % Initial population for S I and R
period = 365*45;
end

%% Setup running time
last_time_point = 0;

%% Run simulation
[DataX VirusesArray] = tauleap_singlesir_ibm_matrix(period,DataX,0,0.5,initViruses); %change from 0.5->0.3
fx = find(VirusesArray(:,1)==0,1); % Because I prelocate a big array
if isempty(fx)
          fx = length(VirusesArray(:,1)) + 1;
end
last_time_point = VirusesArray(fx-1,2);
if last_time_point < period
  disp 'error: processes not finished';
end
disp 'processes completed';
end

