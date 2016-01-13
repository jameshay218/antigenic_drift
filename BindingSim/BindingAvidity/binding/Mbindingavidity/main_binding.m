
global params;
global Viruses;
global VirusArray;
Viruses = [];
params.filename = 'params_std_reinfect.t2000';

S = [];%S0->S21
I = [];%I1->I100

%setup initial population
%initFile = 'dat/dat_x_eq_ini.mat'; %from '../vaccination/out/20120503/std_vac_geo_kc0.18_ks0.075/dat_x_eq.tmp.mat
%params.initFile = initFile;
%dat = open(initFile);
%x1 = dat.x_eq(1,:);%retrieve the first row


%Create Susceiptible groups from S0 to S21
%S  = [99990 zeros(1,99)];
S  = [10000.*ones(1,10) zeros(1,90)];
S(1) = S(1) - 10;
%Create Infected individuals groups from I1 to I20
%for viral binding avidity 1-5.
I = [10 zeros(1,99)];
R = zeros(1,100);


DataX = [S I R];
%tauleap_singlesir(5000,DataX);
%[DataX Viruses] = tauleap_singlesir_ibm_fixedb(100,DataX);


[DataX Viruses] = tauleap_singlesir_ibm_matrix(1000,DataX);
disp 'processes completed';