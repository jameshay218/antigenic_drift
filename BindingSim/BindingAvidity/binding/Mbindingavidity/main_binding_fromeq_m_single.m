%main function to simulate a single outbreak

function [] = main_binding_fromeq_m_single()
global params;
global metadata;
global initViruses;
Viruses = [];

%% Read parameters
%params.filename = 'params_std_reinfect.t2000';
metadata.ibms.filename = 'params_std_reinfect_small.new.dat';
params = InitParameters(['dat/' metadata.ibms.filename]);
metadata.ibms.initFileVirus = 'out/20151122/lar_n4_c15/virus_traits.mat';
metadata.ibms.initFileSIR = 'out/20151122/lar_n4_c15/DataTLSIR.mat';

metadata.ibms.initFlag = 0; 
metadata.ibms.deltaV = 'dat/deltaVMatrix_kc05.mat';
metadata.ibms.antigen = 1;
metadata.ibms.binding = 1;

%% Check parameters
if exist(metadata.ibms.initFileVirus) 
  if exist(metadata.ibms.initFileSIR) 
    metadata.ibms.initFlag = 2; %IBM
  end
else
  if exist(metadata.ibms.initFileODE)   
    metadata.ibms.initFlag = 1; %ODE
  end 
end

%% read Ntot from parameter file
fid = fopen(['dat/' metadata.ibms.filename])
p = textscan(fid,'%s %s');
pa_array_headers = p{1};
pa_array_values = p{2};
N_tot = str2num(cell2mat(pa_array_values(strcmp(pa_array_headers,'N')))); 

%% Define project name
proj = cell2mat(pa_array_values(strcmp(pa_array_headers,'proj'))); 
metadata.ibms.proj = proj;
        
%% Define output folder
out_dir = ['out/' datestr(now,10) datestr(now,5) datestr(now,7) '/' metadata.ibms.proj]
%out_dir = proj;
if(exist(out_dir)==7)
else
   mkdir(out_dir)
end
metadata.ibms.out_dir = out_dir;

%% Setup initial population from IBMS
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
ratio1 = round((sum(I)+sum(S)+ sum(R))./params.N);
S = round(S./ratio1);
I = round(I./ratio1);
R = round(R./ratio1);


%-- move half of R into S
%Rmove = round(R*(1/4));
%Rleft = R - Rmove;
%S = S + [0 R(1:end-1)];
%initViruses(sum(I)+1:end,:) = [];
%if sum(I)~= initViruses(end,1)
%  disp 'error: virus number not consistent';
%end
newinitViruses = [];
for k=1:100
%  I(k) = length(find(initViruses(:,5)==k));
totl=I(k);
idx=find(initViruses(:,5)==k);
newinitViruses=[newinitViruses;initViruses(idx(1:totl),:)];
end
initViruses=[];
newinitViruses(:,1)=1:length(newinitViruses(:,1));
initViruses=newinitViruses;

%totl=length(initViruses(:,1));
%ridx=round(sort(1+(totl-1)*rand(1,sum(I))));
%newinitViruses=initViruses(ridx,:);

%-- Boosting
%R = Rleft;
R_tmp = zeros(1,100);
for k=1:50
AbB = params.AbB;
Boost = poissrnd(AbB,1, sum(R(k)));
for i=1:length(Boost)
R_tmp(k+Boost(i)) = R_tmp(k+Boost(i))+1;
end
end
R = R_tmp; %update recovered individual boosting

DataX = [S I R]; % Initial population for S I and R
clear dat1;
clear dat2;
end

%% Setup initial population from ODE
if metadata.ibms.initFlag == 1
S = [];%S0->S21
I = [];%I1->I100
init_I = 100; %Inital infected individuals 
dat = open(metadata.ibms.initFileODE);
x1 = dat.x_eq(1,1:end);%retrieve the last row
ratio = N_tot./sum(x1);
S = round(x1(1:100).*ratio);
I = [init_I zeros(1,99)];
%I = x1(101:200);
R = round(x1(201:300).*ratio);
S(1) = S(1)+N_tot-(sum(S)+sum(R))-init_I;
DataX = [S I R]; % Initial population for S I and R
end

%% Setup running time
%period = 365*29; %Should change to 45 for 1968-2013
%period = 365*10; %Change to 45 for 1968-2013
period = 365*3;
last_time_point = 0;
%while last_time_point <period

%% Run simulation
[DataY VirusesArray] = tauleap_singlesir_ibm_matrix_antigen(period,DataX,1, 0, initViruses); %change from 0.5->0.3
%[DataY VirusesArray] = tauleap_singlesir_ibm_matrix(period,DataX,1, 0, initViruses); %change from 0.5->0.3
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