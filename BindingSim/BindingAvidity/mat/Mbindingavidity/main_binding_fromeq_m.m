%Simulate viral phylogeny from equilibrium
%Mar 23, 2013

global params;
global inputFile;
global Viruses;
global VirusArray;
global metadata;
metadata.ibms.deltaV = 'dat/deltaVMatrix_kc05.mat';
Viruses = [];



%% Set metadata
%% Define project name
metadata.ode.proj = ['std_adapt_n4_c07']; 
metadata.ibms.proj = ['std_adapt_n4_c07'];

%% Set metadata for ODE
metadata.ode.filename = 'params_std_reinfect_c08_meds.t2000';
metadata.ode.rho = 'transmission_meds.mat';

% Check the existence of the parameter files
if exist(['dat/' metadata.ode.filename])
  disp 'ODE parameters are read.'
end
if exist(['dat/' metadata.ode.rho])
  disp 'ODE transmission array are read.'
end

%% Define output folder
datepath = [datestr(now,10) datestr(now,5) datestr(now,7)];
out_dir = ['out' 'ode/' datepath '/' metadata.ode.proj];
%out_dir = proj;
if(exist(out_dir)==7)
else
   mkdir(out_dir)
end
metadata.ode.out_dir = out_dir;
metadata.ibms.out_dir = out_dir;
%metadata.ibms.filename = ['params_std_reinfect_med.dat'];
metadata.ibms.filename = ['params_std_reinfect_c07_lar.dat'];
metadata.ibms.initFlag = 2; %IBM initFlag=2: input inital viruses



% Read parameters
%params.filename = 'params_std_reinfect.t2000';
inputFile.filename = 'params_std_reinfect_med.dat';
%inputFile.filename = 'params_std_reinfect_small.dat';

% Setup initial population
S = [];%S0->S21
I = [];%I1->I100
init_I = 100; %Inital infected individuals 
%initFile = 'dat/dat_x_eq_ini.mat'; %from '../vaccination/out/20120503/std_vac_geo_kc0.18_ks0.075/dat_x_eq.tmp.mat
%initFile = 'dat/20130812/std_without_vac_c0.5_ks0.0769/dat_x_eq_ini.mat';
initFile = 'dat/20130815/std_adapt_n4_c07/dat_x_eq_ini.mat';
dat = open(initFile);
x1 = dat.x_eq(1,1:end);%retrieve the last row

%initFile = 'dat/20130323/DataTLSIR.mat'; 
%initFile = 'out/20130801/DataTLSIR.mat'; % small N
inputFile.initFile = initFile;
%dat = open(initFile);
%x1 = dat.dat_sir(end,2:end);%retrieve the last row

%modify on May 13, 2013
%To set N_tot > default 1000000
%Thesis version 1000000
%EEID2013 version 1300000
%Final thesis verison 2000000 <- Use this one

fid = fopen(['dat/' inputFile.filename])
p = textscan(fid,'%s %s');
pa_array_headers = p{1};
pa_array_values = p{2};
N_tot = str2num(cell2mat(pa_array_values(strcmp(pa_array_headers,'N')))); 

ratio = N_tot./sum(x1);
S = round(x1(1:100).*ratio);
I = [init_I zeros(1,99)];
%I = x1(101:200);
R = round(x1(201:300).*ratio);
S(1) = S(1)+N_tot-(sum(S)+sum(R))-init_I;
DataX = [S I R]; % Initial population for S I and R

%period = 365*29; %Should change to 45 for 1968-2013
%period = 365*30; %Change to 45 for 1968-2013
period = 365*2;
%period = 365*100; %100 yrs to reach equilibrium for small N
last_time_point = 0;
%while last_time_point <period
%tauleap_singlesir_ibm_matrix(EndTime, DataX, Steps, V_bding)


metadata.ibms.initFileVirus = ['dat/20130906/fb_n4_c07_eq_eq/virus_traits.mat'];
metadata.ibms.initFileSIR = ['dat/20130906/fb_n4_c07_eq_eq/DataTLSIR.mat'];
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



[DataX VirusesArray] = tauleap_singlesir_ibm_matrix(period,DataX,0,0.4, initViruses); %change from 0.5->0.3
fx = find(VirusesArray(:,1)==0);
last_time_point = VirusesArray(fx(1)-1,2);
if last_time_point < period
  disp 'error: processes not finished';
end
%end

disp 'processes completed';
