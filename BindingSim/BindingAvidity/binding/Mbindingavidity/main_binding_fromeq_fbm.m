%Simulate viral phylogeny from equilibrium
%Mar 23, 2013

global params;
global inputFile;
global Viruses;
global VirusArray;
Viruses = [];

% Read parameters
%params.filename = 'params_std_reinfect.t2000';
inputFile.filename = 'params_std_reinfect_lar.dat';

% Setup initial population
S = [];%S0->S21
I = [];%I1->I100
init_I = 100; %Inital infected individuals 
initFile = 'dat/20130827/std_fixb_c05/dat_x_eq_ini.mat';
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
period = 365*1; %Change to 45 for 1968-2013
%period = 365*0.5;
%period = 365*100; %100 yrs to reach equilibrium for small N
last_time_point = 0;
%while last_time_point <period
%tauleap_singlesir_ibm_matrix(EndTime, DataX, Steps, V_bding)
[DataX VirusesArray] = tauleap_singlesir_ibm_matrix_fixedb(period,DataX,0,0.4,0.4); %change from 0.5->0.3
fx = find(VirusesArray(:,1)==0);
last_time_point = VirusesArray(fx(1)-1,2);
if last_time_point < period
  disp 'error: processes not finished';
end
%end

disp 'processes completed';
