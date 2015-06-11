%% Declare variables
global params;
global VirusArray;
params.filename = 'params.dat';
%% Set file names
metadata.ibms.proj = ['std'];
metadata.ibms.deltaV = 'dat/deltaVMatrix_kc05.mat';
metadata.ibms.parameterFile = ['dat/' params.filename];
metadata.ibms.initVirusFlag = true; 
%% Define output folder
datepath = [datestr(now,10) datestr(now,5) datestr(now,7)];
out_dir = ['out/' 'ibms/' datepath '/' metadata.ibms.proj];
if(exist(out_dir)==7)
else
   mkdir(out_dir)
end
metadata.ibms.out_dir = out_dir;
%% Initialize Parameters
params = InitParameters(metadata.ibms.parameterFile);
['read parameter file from ' metadata.ibms.parameterFile];
params.out_dir = metadata.ibms.out_dir;
%% Initilize States
%Create Susceiptible groups from S(k=0) to S(k=N_Infect-1)
S  = zeros(1,params.N_Infect);
I = zeros(1,params.N_Infect);
R = zeros(1,params.N_Infect);
N = params.N;
N_Infect = params.N_Infect;
seed= params.seed;
k = 2;
S  = [(N./k).*ones(1,k) zeros(1,N_Infect-k)];
S(1) = S(1) - seed;
I = [seed zeros(1,N_Infect-1)];
R = zeros(1,N_Infect);
DataX = [S I R]; %X0
%% Run tauleap algorithm
time_step = 1;
total_days = 365*0.5;
[DataX Viruses] = tauleap_singlesir_ibm(total_days,DataX,time_step,metadata);
%% Plotting simulation log
cal_meanbding(metadata);
disp 'plot simulation log';
plot_simulation('out/ibms/20150526/std')
disp 'processes completed';