%% Declare variables
global params;
global Viruses;
global VirusArray;
Viruses = [];
params.filename = 'params.dat';
%% Set file names
metadata.ibms.proj = ['std'];
metadata.ibms.deltaV = 'dat/deltaVMatrix_kc05.mat';
metadata.ibms.filename = ['dat/params.dat'];
metadata.ibms.initViruses = true; 
%% Define output folder
datepath = [datestr(now,10) datestr(now,5) datestr(now,7)];
out_dir = ['out/' 'ibms/' datepath '/' metadata.ibms.proj];
if(exist(out_dir)==7)
else
   mkdir(out_dir)
end
metadata.ibms.out_dir = out_dir;
%% Initialize Parameters
filename = metadata.ibms.filename;
params = InitParameters(filename);
['read parameter file from ' filename];
params.out_dir = metadata.ibms.out_dir;
%% Initilize Statatus
%Create Susceiptible groups from S0 to S21
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
DataX = [S I R];
%% Run tauleap algorithm
[DataX Viruses] = tauleap_singlesir_ibm_matrix(365*0.5,DataX,1,metadata);
cal_meanbding(metadata);
disp 'plot simulation log';
plot_simulation('out/ibms/20150526/std')
disp 'processes completed';