function [DataX Viruses] = tauleap_singlesir_ibm_matrix(Steps, DataX)
%Tauleap method for single strain SIR reinfection model
%Transmission rate b is constant
%When any variable becomes negative, devides the time_step by 2
%Input: [Steps, DataX]
%   Steps: total number of steps
%   DataX: initial condition
%Output: [TimeStep DataX]
%   DataX: data status
%   Viruses: Virus strains  
%Use predetermined matrix of deltaV
%Written on Aug 9, 2012


%List of Events:
%Birth:         uN,             S=S+1  
%Transmission:  beta(S)(I),     I=I+1, S=S-1
%Recovery:      g(I), I=I-1,    S=S+1
%Death:         -uS,-uI,        S=S-1, I=I-1
%Wanning:       wR,             R=R-1, S=S+1

%%%------------------------------------------------------------------------
%%%--------------------Declare global variables----------------------------
%%%------------------------------------------------------------------------
global params;      %Store parameters
global Viruses;     %Store virus stains in an array of structures
global VirusArray;  %Store currently existed virus strains
maxvid = [0];       %Tracking the maximum number of virus ID
rn = Steps;            %Save every 20 cycles

%Load deltaV Matrix
dat = open('dat/deltaVMatrix.mat');
deltaVMatrix = dat.deltaVMatrix;

%Metadata: column names of Viruses
VirusCol = struct('vid',1,'birth',2,'death',3,'parent',4,'infectionK',5,'beta',6,'initialV',7,'currentV',8);

%Initialize Parameters
filename = params.filename;
params = InitParameters(['dat/' filename]);
['read parameter file from ' 'dat/' filename]
%Time unit
time_step = 0.5;
iterations = Steps;

%Rate matrix
Rate_Matrix = [];
Rate_Matrix_Str = {}; %store the events names
Data_Events_rate = [];

%number of Sk, Ik, and Rk
S = [];%%Create Susceptible groups from S0-S99
I = [];%%Create Infected groups from I1-I100
R = [];%%Create Recovered groups from R1-R100
Beta = [];
S = DataX(1,1:params.N_Infect)';   
I = DataX(1,params.N_Infect+1:params.N_Infect*2)';  %20x5 I groups
R = DataX(1,params.N_Infect*2+1:params.N_Infect*3)'; %20 R groups
Sk = [0:params.N_Infect-1];

%Initialize virus strains
VirusArray = []; %current viral strains
%for i=1:10
Iini = sum(I);
for i=1:Iini
k=1; %Still working. Beta is a function of Sk and V
%v=0.1;
v=0.5;
if ~isempty(Viruses) && isstruct(Viruses)
   beta = get_beta(Sk, v, params.p, params.r, params.a, params.b, params.c);
   Viruses(i) = createVirus(i, 0, 0, 0, 1, beta, v, v);
else % initialize
   beta = get_beta(Sk, v, params.p, params.r, params.a, params.b, params.c);
   Viruses = createVirus(i, 0, 0, 0, 1, beta, v, v); %(vid, birth, death, sourceVirus, infectionK, beta, initialV, currentV) 
end
beta = get_beta(Sk, v, params.p, params.r, params.a, params.b, params.c);
VirusArray(i,:) = [i 0 0 0 1 beta(1) v v];%vid, birth, death, source, infectionK, beta, initialV, currentV %%%%%%%% Can't put array or structure into an element in an array 
maxvid = maxvid + 1;
end

%%%------------------------------------------------------------------------
%%%---------------------Simulation of the events---------------------------
%%%------------------------------------------------------------------------
PrevTime = 0;
CurrTime = 1;

%First column -> time
DataS = [CurrTime S'];
DataI = [0 I'];
DataR = [0 R'];
Data = [DataS DataI(2:params.N_Infect+1) DataR(2:params.N_Infect+1)];

%Convert to Beta matrix
%Beta(k,vi)  k=0-99, vi={v1-vn} n=number of infected individuals  
for i=1:length(VirusArray(:,VirusCol.vid))
Beta(:,i) = Viruses(VirusArray(i,VirusCol.vid)).beta';
end


for i1=2:iterations
 CurrTime
 time_step1 = Iterate(params);
 
 %update DataS, DataI and DataR matrix
 %DataS(end,1) = DataS(end-1,1)+time_step;
 DataS(end+1,1) = CurrTime;
 DataS(end,2:params.N_Infect+1) = S';
 %DataI(end,1) = DataI(end-1,1)+time_step;
 DataI(end+1,1) = CurrTime;
 DataI(end,2:params.N_Infect+1) = I';
 %DataR(end,1) = DataR(end-1,1)+time_step;
 DataR(end+1,1) = CurrTime;
 DataR(end,2:params.N_Infect+1) = R';
 %Data(end,1) = DataS(end-1,1)+time_step;
 Data(end+1,1) = CurrTime;
 Data(end,2:params.N_Infect*3+1) = [S' I' R'];
 
 if rem(i1,rn) == 0
    output_viruses(Viruses);
    Viruses(find([Viruses.death]~=0)) = []; %clear viruses after saving them
    output_sir(Data);
    disp "export"
    Data = [];
    DataS = [];
    DataI = [];
    DataR = [];
 end
end
return;
%%%------------------------End of Simulation-------------------------------




%%%------------------------------------------------------------------------
%%%---------------------Iteration of tao leap Algorithm--------------------
%%%------------------------------------------------------------------------
function [time_step1]=Iterate(Parameters)
time_step1=time_step;

%Updating the binding avidity dynamics
%Event of binding avidity changes
disp 'updating binding avidity';
vpop = length(VirusArray(:,VirusCol.vid))

if vpop < 2
   output_viruses(Viruses);
    Viruses(find([Viruses.death]~=0)) = []; %clear viruses after saving them
    output_sir(Data);
    disp "export"
end

%vidlist = VirusArray(:,VirusCol.vid);
%rowno_vid = find([Viruses.vid]==vidlist);
%%%Try to see whether there is a way to make it faster
vcurr = [];
rowno_list = [];
%disp '1';
for i=1:length(VirusArray(:,VirusCol.vid))
    vid = VirusArray(i,VirusCol.vid);
    rowno_vid = vid; %before clearing the Viruses structure, vid is same as rowno
    
    %Don't clear the Viruses, will produce errors
    %rowno_vid = find([Viruses.vid]==vid);%after clearing the Viruses
    %structure, vid need to be find
    %if isempty(rowno_vid)
    %  rowno_vid
    %end
    rowno_list(i) = rowno_vid;
    virus = Viruses(rowno_vid);
    vini = virus.currentV;
    k = virus.infectionK;
    if isempty(k) 
      disp "error: check viruses array";  %%%%%%show error message if there are empty values in records
    end
    %vcurr = getVChange_ode(vini,k,time_step1);
    vcurr(i) = vini + getDeltaV(vini,k); %% This causes the speed slow
    %updating the transmission rate beta    
    %%Viruses(rowno_vid).beta = beta;
    %Viruses(rowno_vid).beta = 0.5*ones(1,100);
    %%Viruses(rowno_vid).currentV = vcurr;
end
%beta(k,v), k=0-N_Infect-1,v=number of viral strains
%disp '2';
beta = get_beta_list(params.N_Infect-1,vcurr, params.p, params.r, params.a, params.b, params.c);
betamat = beta;
for i=1:length(rowno_list)
  Viruses(rowno_list(i)).beta = beta(:,i)';
  Viruses(rowno_list(i)).currentV = vcurr(i);
end
%disp '3';

%Rate of each epidemiological events
Rate_Birth = params.mu*params.N;
Rate_Death_S = params.mu*S; 
Rate_Death_I = params.mu*I;
Rate_Death_R = params.mu*R;
Rate_Recovery = params.gamma*I;
Rate_Wanning = params.wan*R;

I_total = sum(I);
%si = S*I_total;  



%%Store all Rate into a single rate matrix
%Rate_Matrix = [Rate_Matrix; Rate_Birth];
%Rate_Matrix = [Rate_Matrix; Rate_Death_S];
%Rate_Matrix = [Rate_Matrix; Rate_Death_I];
%Rate_Matrix = [Rate_Matrix; Rate_Recovery];
%Rate_Matrix = [Rate_Matrix; Rate_Infection];

%%%Number of events happened
%%%%Birth
dBirth = poissrnd(Rate_Birth*time_step);
%%%Death
%Death of Susceptible
dDeath_S = poissrnd(Rate_Death_S*time_step);
%Death of Infecteds
%Death of Recovered
dDeath_R = poissrnd(Rate_Death_R*time_step);

%disp '4';
%Death of Infecteds 
dDeath_I = zeros(params.N_Infect,1);
I_indv = binornd(1,params.mu*time_step, length(VirusArray(:,1)), 1);
de_id = find(I_indv==1);
for i=1:length(de_id)
        k = VirusArray(de_id(i),VirusCol.infectionK);
        dDeath_I(k) = dDeath_I(k)+1;
end
%Updating Death from infection on VirusArray
de_id = deathOfInfectds(CurrTime, de_id);

%disp '5';
%%%Recovery
dRecover = zeros(params.N_Infect,1);
R_indv = binornd(1,params.gamma*time_step, length(VirusArray(:,1)), 1);
rm_id = find(R_indv==1);
disp '5-1';
for i=1:length(rm_id)
        k = VirusArray(rm_id(i),VirusCol.infectionK);
        dRecover(k) = dRecover(k)+1;
end
%disp '5-2';
%Updating recovery (death time) for each viral strain
rm_id = recovery(CurrTime, rm_id);

%disp '6'; %%%codes below run very slow 
%%%Infection
%disp 'calculating infection';
Beta = betamat; %Just need to retrieve the beta matrix from the codes before
%Beta = getBeta();
%disp '7';
%for i=1:length(Beta(1,:))
Rate_Infection = repmat(S,[1 length(Beta(1,:))]).*Beta/params.N;
%end
dInfection_indv = poissrnd(Rate_Infection*time_step);
dInfection = sum(dInfection_indv,2);
dInfection_tot = sum(dInfection);

%%%Wanning
dWanning = poissrnd(Rate_Wanning*time_step);
dWanning_S = [0; dWanning(1:params.N_Infect-1)];

%disp '8';
%%Check whether any negative values from poisson distribution 
if isNegative(dBirth, dDeath_S, dInfection_tot, dInfection, dDeath_I, dRecover, dDeath_R, dWanning, dWanning_S) == 1
        time_step = time_step/2;
        return; %stop updating and return to start if there is negative values 
end

%disp '9';
%%Updating the status
%disp 'updating the status';
PrevTime = CurrTime;
CurrTime = PrevTime + time_step; %Tracking the current time
S(1)=S(1)+dBirth; %1) Birth
S=S-dDeath_S;     %2) Death from susceptible individuals
I=I-dDeath_I;     %3) Death from infected individuals


R=R-dDeath_R;     %4) Death from recovered individuals
I=I+dInfection;   %5) New Infectious individuals
%I(1,96:100) = I(1,96:100)+dInfection_List(1,101:105); %S20 reinfected to I20
S=S-dInfection;   %6) Loss of susceptible individuals from Infection
I=I-dRecover;     %7) Loss of infected individuals from Recovery
R=R+dRecover;     %8) New recovered individuals
R=R-dWanning;     %9) Loss if recovered individuals after immune wanning 
S=S+dWanning_S;   %10)Susceptible individuals from immune waning


%Updating Death from infection
%VirusArray(de_id,:) = []; %alreadt removes records in ln174


SourceViruses = VirusArray(:, VirusCol.vid);

%Updating transmission for each viral strain  
disp 'updating the transmission history';
for k=1:params.N_Infect %each k
    %NewVirusesID = find(dInfection_indv(k,:)>0);
    %SourceViruses = VirusArray(NewVirusesID, VirusCol.vid);
    %SourceViruses = VirusArray(:, VirusCol.vid);
    %NewViruses = dInfection_indv(k,find(dInfection_indv(k,:)>0));
    if sum(dInfection_indv(k,:))==0 %If there is no infection for Sk, skip this k
        continue;
    end
    v_size = length(Viruses);
    for ns=1:length(SourceViruses) %each source viruses     
        for nv=1:dInfection_indv(k,ns) %new viruses from the same source virus
         %vid = v_size+1;% this only current before clear Viruses array
         vid = maxvid+1;
         %(find([Viruses.vid]==vid)
         
         rowno_vid = SourceViruses(ns);
         vcurr = Viruses(rowno_vid).currentV; %before clearing viruses structure, rowno = vid
         %vcurr = Viruses(find([Viruses.vid]==SourceViruses(ns))).currentV;
         
         %vcurr = Viruses(SourceViruses(ns)).currentV %retrieve source viruses binding V
         %vcurr = Viruses(VirusArray(NewViruses(nv),VirusCol.vid)).currentV;
         Viruses(v_size+1)=createVirus(vid, CurrTime, 0, SourceViruses(ns), k, zeros(1,100), vcurr, vcurr); %createVirus( vid, birth, death, sourceVirus, infectionK, beta, initialV, currentV )
         VirusArray(length(VirusArray(:,1))+1,[VirusCol.vid VirusCol.infectionK VirusCol.beta]) = [vid,k,0]; %insert virus id into VirusArray (vid, infecitonK, beta)
         maxvid = maxvid+1; %keep updating maxvid
         v_size = length(Viruses);
        end
    end
end
%disp '10';
time_step1 = time_step;
end


function [beta] = getBeta()
    beta = [];
    %Convert to Beta matrix
    %%% Need testing
    for i=1:length(VirusArray(:,VirusCol.vid))
        vid = VirusArray(i,VirusCol.vid);
        beta(:,i) = Viruses(find([Viruses.vid]==vid)).beta';
    end 
end


function [de_id] = deathOfInfectds(CTime, de_id)
    for i=1:length(de_id)
        VirusArray(de_id(i),VirusCol.death) = PrevTime + time_step;
        %Viruses(VirusArray(de_id(i),1)).death = PrevTime + time_step;
        vid = VirusArray(de_id(i),VirusCol.vid);
        Viruses(find([Viruses.vid]==vid)).death = PrevTime + time_step;
    end
    VirusArray(de_id,:) = [];
end

function [rm_id] = recovery(CTime, rm_id)
    %R = binornd(1,0.5, length(VirusArray), 1)
    %rm_id = find(R==1);
    for i=1:length(rm_id)
        %i
        %rm_id(i)
        rowno_vid = VirusArray(rm_id(i),VirusCol.vid);
        %length(VirusArray)
        %Viruses(find([Viruses.vid]==vid))
        %if length(find([Viruses.vid]==rowno_vid)) > 1
        %   disp 'something is wrong';
        %end
        Viruses(rowno_vid).death = CTime;
        %Viruses(find([Viruses.vid]==vid)).death = CTime;
    end
    VirusArray(rm_id,:) = [];
end


%Check whether any negative values in the array
function [pos]=isNegative(dBirth, dDeath_S, dInfection_tot, dInfection, dDeath_I, dRecover, dDeath_R, dWanning, dWanning_S)
  pos = -1; 
  if(checkNegative(S(1)+dBirth)==1)
    pos = 1;
    return;
  end
  dS = S-dDeath_S-dInfection+dWanning_S;
  %if(checkNegative(S-dDeath_S-dInfection+dWanning_S)==1)
  if(checkNegative(dS)==1)
    pos = 1;
    return;
  end

  %if(checkNegative(S-dInfection)==1)
  %  pos = 1;
  %  return;
  %end

  dI = I-dDeath_I+dInfection-dRecover;
  if(checkNegative(dI)==1)
     disp('change timestep to avoid negative value of I');      
     pos = 1;
     return;
  end
  
  dR = R-dDeath_R+dRecover-dWanning;
  if(checkNegative(dR)==1)
     disp('change timestep to avoid negative value of R');      
     pos = 1;
     return;
  end
  
  function [pos]=checkNegative(x)
  pos = -1;
  sz = size(x);
  xlist = reshape(x,1,sz(1)*sz(2));
  for i=1:length(xlist)
      if(xlist(i)<0)
          pos=1;
          return;
      end
  end
  end
end

function output_sir(Data)
    % create output directory
    out_dir = ['out/' datestr(now,10) datestr(now,5) datestr(now,7)];
    %out_dir = params.out_dir
    if(exist(out_dir)==7)
    else
        mkdir(out_dir)
    end
    if exist([out_dir '/DataTLSIR.mat'], 'file')
        dat = open([out_dir '/DataTLSIR.mat']);
        dat_sir = dat.dat_sir;
        for i=1:length(Data(:,1))
            dat_sir = [dat_sir; Data(i,:)];
        end
        save([out_dir '/DataTLSIR.mat'],'dat_sir');
        clear dat_sir;
    else
        dat_sir = [];
        for i=1:length(Data(:,1))
            dat_sir(i,:) = Data(i,:);
        end
        save([out_dir '/DataTLSIR.mat'],'dat_sir');
        disp('done');
        clear dat_sir;
    end  
end


function output_viruses(Viruses)
    % create output directory
    out_dir = ['out/' datestr(now,10) datestr(now,5) datestr(now,7)];
    %out_dir = params.out_dir
    if(exist(out_dir)==7)
    else
        mkdir(out_dir)
    end
    if exist([out_dir '/dat_x_trans.tmp.mat'], 'file')
        dat = open([out_dir '/dat_x_trans.tmp.mat']);
        dat_viruses = dat.dat_viruses;
        for i=1:length(Viruses)
            dat_viruses(Viruses(i).vid,:) = [Viruses(i).vid, Viruses(i).birth, Viruses(i).death, Viruses(i).parent, Viruses(i).infectionK];
        end
        save([out_dir '/dat_x_trans.tmp.mat'],'dat_viruses');
        clear dat_viruses;
    else
        dat_viruses = [];
        for i=1:length(Viruses)
            dat_viruses(Viruses(i).vid,:) = [Viruses(i).vid, Viruses(i).birth, Viruses(i).death, Viruses(i).parent, Viruses(i).infectionK];
        end
        save([out_dir '/dat_x_trans.tmp.mat'],'dat_viruses');
        disp('done');
        clear dat_viruses;
    end  
end

%get estimate deltaV from a matrix
function [deltaV] = getDeltaV(vini,k)
   vlist = 0:1:60;
   vini_r=round(vini*10)+1;
   if isempty(find(vlist==vini_r))
     vini_r
     disp 'error during calculating deltaV';
   end
   deltaV=deltaVMatrix(k,find(vlist==vini_r));
end

%Plot population of different viruses 
function plot_binding_avidity(x1,t)
    figure;
    hold;
    meancharge_before = calculate_mean_proportion(x1);
    plot(t, meancharge_before);
end    

function plot_proportion(x)
    i = x(length(x),22:121); %retrieve data at the last date
    b = (reshape(i, 5, 20))';
    c = sum(b);
    ib = c./sum(c);
    figure;
    plot(ib(1,:),'go-', 'LineWidth', 1);
end 

function plot_si(x)
figure;
subplot(3,2,1);
%figure;
plot(x(1:length(x),1),x(1:length(x),2:6));
title('Number of Susceptible Individuals');
hleg1 = legend('S0','S1','S2','S3','S4');

subplot(3,2,2);
plot(x(1:length(x),1),x(1:length(x),7:11));
hleg2 = legend('S5','S6','S7','S8','S9');

%figure;
subplot(3,2,3);
plot(x(1:length(x),1),x(1:length(x),12:16));
hleg3 = legend('S10','S11','S12','S13','S14');
subplot(3,2,4);
plot(x(1:length(x),1),x(1:length(x),17:21));
hleg4 = legend('S15','S16','S17','S18','S19');

subplot(3,2,5);
plot(x(1:length(x),1),x(1:length(x),22:22));
hleg5 = legend('S20');

subplot(3,2,6);
p=x(x(1:length(x),1),1:length(x),2:22);
plot(sum(p'));
hleg6 = legend('Total (S)');
end

end






