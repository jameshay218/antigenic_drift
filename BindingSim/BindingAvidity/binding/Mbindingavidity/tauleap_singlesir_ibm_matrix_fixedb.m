function [DataX VirusesArray] = tauleap_singlesir_ibm_matrix_fixedb(EndTime, DataX, Steps, V_bding, initViruses)
%Tauleap method for single strain SIR reinfection model
%Transmission rate b is is stored in Viruses.beta
%When any variable becomes negative, devides the time_step by 2
%Input: [EndTime, DataX, Steps]
%   EndTime: time period
%   DataX: initial condition
%   Steps: total number of steps
%Output: [TimeStep DataX]
%   DataX: data status
%   CurrentViruses: Current active viruses (replaces VirusArray)
%   Viruses: {structure} Virus strains, contains beta array
%   VirusesArray: store traits for all active and hitory viruses 
%Clean version for bugs free
%Use predetermined matrix of deltaV
%Written on Mar 22, 2013
%Update death and recovery at the SAME time after update transmission events
%Otherwise few data will be lost
%v7, deletete Viruses strucure

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
global metadata;    %Store input files names
global inputFile;   %Store input files names
global Viruses;     %Store virus stains in an array of structures of 'vid','birth','death','parent','infectionK',5,'beta'
global VirusesArray;    %Store all historical viruses in array
global CurrentViruses;  %Store currently existed virus strains in array
maxvid = [0];       %Tracking the maximum number of virus ID

%VirusesArray = []; Need to prelocate the memory, otherwise very slow
VirusesArray = zeros(50E6,8); %what is the max capacity and max population size 20E6 too large
%May 11, 2013 11E6
CurrentViruses = [];

%Load deltaV Matrix
%dat = open('dat/deltaVMatrix.mat');
dat = open(metadata.ibms.deltaV);
deltaVMatrix = dat.deltaVMatrix;

%Metadata: column names of Viruses
VirusCol = struct('vid',1,'birth',2,'death',3,'parent',4,'infectionK',5,'beta',6,'initialV',7,'currentV',8);

%Initialize Parameters
filename = metadata.ibms.filename;
%filename = inputFile.filename;
params = InitParameters(['dat/' filename]);
['read parameter file from ' 'dat/' filename];
params.out_dir = metadata.ibms.out_dir;

%Load fixed rho
%dat1 = open('dat/Transmission.mat'); %replaced by metadata
dat1 = open(metadata.ibms.rho);
rho = dat1.rho;
params.rho = rho;

%Time unit
time_step = 1; % Change to 1day to simulate 45 years. Jul 16,2013
%iterations = Steps;

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
CurrentViruses = []; %current viral strains
if metadata.ibms.initFlag == 1
%for i=1:10
Iini = sum(I,2);
for i=1:length(Iini)
    for j=1:Iini(i)
    k=i; %Still working. Beta is a function of Sk and V
    %v=0.1;
    if exist('V_bding')
        v = V_bding;
    else
        v=0.5;
    end

    vid = maxvid+1;
    if ~isempty(Viruses) && isstruct(Viruses)
        %beta = get_beta(Sk, v, params.p, params.r, params.a, params.b, params.c, params.nv);
        %beta = b0;
        beta = params.c*params.rho(1);
        Viruses(vid) = createVirus2(vid, 0, 0, 0, k, beta, v, v);
    else % initialize
        %beta = get_beta(Sk, v, params.p, params.r, params.a, params.b, params.c, params.nv);
        %beta = b0;
        beta = params.c*params.rho(1);
        Viruses = createVirus2(vid, 0, 0, 0, k, beta, v, v); %(vid, birth, death, sourceVirus, infectionK, beta, initialV, currentV) 
    end

    %beta = get_beta(Sk, v, params.p, params.r, params.a, params.b, params.c, params.nv);
    %beta = b0;
    beta = params.c*params.rho(1); % transmission rate beta for S0
    birth = 0;
    death = 0;
    parent = 0;
    infectionK = i;
    beta = beta(1);
    initialV = v;
    currentV = v;
    CurrentViruses(vid,:) = [vid, birth, death, parent, infectionK, beta, initialV, currentV];
    VirusesArray(vid,:) = [vid, birth, death, parent, infectionK, beta, initialV, currentV];
    maxvid = maxvid + 1;
    end
end
elseif metadata.ibms.initFlag == 2
    CurrentViruses = initViruses;
    siz = length(initViruses(:,1)); %%%% This keeps the array size
    VirusesArray(1:siz,:) = initViruses;
    maxvid = siz;
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

Viruses = [];

while 1
 CurrTime;
 time_step1 = Iterate(params);
 
 vpopsize = length(CurrentViruses(:,VirusCol.vid));

 if vpopsize < 2
    %output_viruses(Viruses);
    %Viruses(find([Viruses.death]~=0)) = []; %clear viruses after saving them
    %output_sir(Data);
    disp('Viruses extinction');
    return;
 end

 %update DataS, DataI and DataR matrix
 DataS(end+1,1) = CurrTime;
 DataS(end,2:params.N_Infect+1) = S';
 DataI(end+1,1) = CurrTime;
 DataI(end,2:params.N_Infect+1) = I';
 DataR(end+1,1) = CurrTime;
 DataR(end,2:params.N_Infect+1) = R';
 Data(end+1,1) = CurrTime;
 Data(end,2:params.N_Infect*3+1) = [S' I' R'];
 
if CurrTime>1 & rem(CurrTime,20) == 0
 CurrTime
 vpopsize
 %PopTotal = sum(S)+sum(I)+sum(R)
end
  if CurrTime >= EndTime
    output_viruses2(VirusesArray);
    output_sir(Data);
    output_traits(CurrTime);
    disp "export"
    Data = [];
    DataS = [];
    DataI = [];
    DataR = [];
    return;
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
vpop = length(CurrentViruses(:,VirusCol.vid));



%if vpop < 2
%    output_viruses(Viruses);
%    %%Viruses(find([Viruses.death]~=0)) = []; %clear viruses after saving them
%    output_sir(Data);
%    error('No viruses anymore')
%end


vcurr = [];
rowno_list = [];


%Updating binding avidity V and beta
vid_list = CurrentViruses(:,VirusCol.vid);
vini_list = VirusesArray(vid_list,VirusCol.currentV);
k_list = VirusesArray(vid_list,VirusCol.infectionK);
vcurr_list = vini_list + getDeltaV(vini_list, k_list).*time_step;

%if (CurrTime<180)
%    beta = get_beta_list(params.N_Infect-1,vcurr_list', params.p, params.r, params.a, params.b, 0.45, params.nv);
%else
%end
    %beta1 = get_beta_list(params.N_Infect-1,vcurr_list', params.p, params.r, params.a, params.b, params.c, params.nv);
    %beta = b0*ones(params.N_Infect,length(vcurr_list));
    beta = repmat(params.c*params.rho',[1 length(vcurr_list)]);
betamat = beta;
rowno_list = vid_list;
beta_c = num2cell(beta(:,1:length(vid_list))'); %same #row of currentviruses
vcurr_c = num2cell(vcurr_list);
VirusesArray(vid_list(:),VirusCol.beta) = beta(1,:); %save in vid, only store beta for naive
VirusesArray(vid_list(:),VirusCol.currentV) = vcurr_list(:,1);

CurrentViruses(:,VirusCol.beta) = beta(1,:); %save in vid, only store beta for naive
CurrentViruses(:,VirusCol.currentV) = vcurr_list(:,1);
%for i=1:length(rowno_list)
%  Viruses(i).beta = beta(:,i)'; %same #row of currentviruses
%  Viruses(i).currentV = vcurr_list(i,1);
%  VirusesArray(rowno_list(i),VirusCol.beta) = beta(1,i); %save in vid, only store beta for naive
%  VirusesArray(rowno_list(i),VirusCol.currentV) = vcurr_list(i,1);
%end

%Rate of each epidemiological events
Ntot = sum(S)+sum(I)+sum(R);
%Rate_Birth = params.mu*params.N;
Rate_Birth = params.mu*Ntot;
Rate_Death_S = params.mu*S; 
Rate_Death_I = params.mu*I;
Rate_Death_R = params.mu*R;
Rate_Recovery = params.gamma*I;
Rate_Wanning = params.wan*R;
I_total = sum(I);

%%%Number of events happened
%%%Birth
dBirth = poissrnd(Rate_Birth*time_step);

%%%Death
%Death of Susceptible
dDeath_S = poissrnd(Rate_Death_S*time_step);
%Death of Recovered
dDeath_R = poissrnd(Rate_Death_R*time_step);
%Death of Infecteds 
dDeath_I = zeros(params.N_Infect,1);

%%I_indv = poissrnd(Rate_Death_I*time_step);
%%I_indv = binornd(1,params.mu*time_step, length(CurrentViruses(:,1)), 1);
%%If death and recovered happens for the same virus, what next?
%%use multinomial

prob = [params.mu*time_step,params.gamma*time_step,1-params.mu*time_step-params.gamma*time_step];
I_loss = mnrnd(1,prob,length(CurrentViruses(:,1)));
I_indv = I_loss(:,1);
R_indv = I_loss(:,2);

de_id = find(I_indv==1); % is this right
d_list = CurrentViruses(de_id, VirusCol.infectionK); %return a list of k
d_uni_list = unique(d_list);
for i=1:length(d_uni_list)
    k = d_uni_list(i);
    dDeath_I(k) = dDeath_I(k)+length(d_list(d_list==k));
end

%%%Recovery
dRecover = zeros(params.N_Infect,1);
%R_indv = poissrnd(Rate_Recovery*time_step);
%R_indv = binornd(1,params.gamma*time_step, length(CurrentViruses(:,1)), 1); %recovered_ID, determines which infected individuals are going to be recovered
rm_id = find(R_indv==1);
k_list = CurrentViruses(rm_id, VirusCol.infectionK); %return a list of k
k_uni_list = unique(k_list);
for i=1:length(k_uni_list)
    k = k_uni_list(i);
    dRecover(k) = dRecover(k)+length(k_list(k_list==k));
end


%%%Infection
Beta = betamat; %Just need to retrieve the beta matrix from the codes before
Rate_Infection = repmat(S,[1 length(Beta(1,:))]).*Beta/params.N;
dInfection_indv = poissrnd(Rate_Infection*time_step);
dInfection = sum(dInfection_indv,2);
dInfection_tot = sum(dInfection);

%%%Wanning
dWanning = poissrnd(Rate_Wanning*time_step);
dWanning_S = [0; dWanning(1:params.N_Infect-1)];

%%Check whether any negative values from poisson distribution 
if isNegative(dBirth, dDeath_S, dInfection_tot, dInfection, dDeath_I, dRecover, dDeath_R, dWanning, dWanning_S) == 1
        %time_step = time_step/2;
        return; %stop updating and return to start if there is negative values 
end


%%Updating the status
PrevTime = CurrTime;
CurrTime = PrevTime + time_step; %Tracking the current time
%prevI=sum(I);
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
if CurrTime>1 & rem(CurrTime,10) == 0
 CurrI=sum(I)
end






SourceViruses = CurrentViruses(:, VirusCol.vid);

%Updating transmission for each viral strain  
number_new_viruses = sum(dInfection_indv); %sum of new viruses for each current active virus
ns_list = find(number_new_viruses>0); %the index of the viruses who produce offsprings

    %only search for viruses who produce offsprings
    for i =1:length(ns_list) %each source viruses     
        ns = ns_list(i);
        %if sum(dInfection_indv(:,ns))==0 %If there is no infection for this virus, skip it
        %    continue;
        %end
        vcurr = CurrentViruses(ns,VirusCol.currentV);
        number_of_infected_by_virus = sum(dInfection_indv(:,ns));
        k_list1 = find(dInfection_indv(:,ns)>0);
        vid_array = [];
        k_array = [];
        %for k=1:params.N_Infect %each k
        for i=1:length(k_list1) %each k
            k = k_list1(i);
            number_of_infected = dInfection_indv(k,ns);
            vid_array(end+1:end+number_of_infected) = maxvid+1:maxvid+number_of_infected;
            k_array(end+1:end+number_of_infected) = k.*ones(1,number_of_infected);
            maxvid = maxvid+number_of_infected;
        end
        vid_j = vid_array;
        birth_j = CurrTime;
        death_j = 0;
        parent_j = SourceViruses(ns);
        infectionK_j = k_array;
        beta_j = zeros(1,100);
        initialV_j = vcurr; 
        currentV_j = vcurr; 

    
        if number_of_infected_by_virus>1
          number_of_infected_by_virus;
        end
        %Insert into CurrentViruses
        CurrentViruses(end+1:end+number_of_infected_by_virus,[VirusCol.vid VirusCol.infectionK VirusCol.beta]) = [vid_j',infectionK_j',zeros(1,number_of_infected_by_virus)']; %insert virus id into CurrentViruses (vid, infecitonK, beta)
        %Insert into VirusesArray
        VirusesArray(vid_j,:) = [vid_j', birth_j.*ones(1,number_of_infected_by_virus)', death_j.*ones(1,number_of_infected_by_virus)', parent_j.*ones(1,number_of_infected_by_virus)', infectionK_j', beta_j(1).*ones(1,number_of_infected_by_virus)', initialV_j.*ones(1,number_of_infected_by_virus)', currentV_j.*ones(1,number_of_infected_by_virus)']; %insert virus id into CurrentViruses (vid, infecitonK, beta)
    end
     

%Updating Death from infection on CurrentViruses
de_id1 = deathOfInfectds(PrevTime, de_id); %death_ID, determines which infected individuals are going to be died

%Updating recovery (death time) for each viral strain
rm_id1 = recovery(PrevTime, rm_id);

removeViruses(de_id, rm_id);
CurrentViruses;
%disp '10';
time_step1 = time_step;

end
%%---end of iteration---###

function [de_id] = deathOfInfectds(CTime, de_id)
    %for i=1:length(de_id)
    %    CurrentViruses(de_id(i),VirusCol.death) = CTime;
    %    vid = CurrentViruses(de_id(i),VirusCol.vid);
    %    VirusesArray(find(VirusesArray(:,VirusCol.vid)==vid),VirusCol.death) = CTime;
    %end
    CurrentViruses(de_id,VirusCol.death) = CTime;
    vid_list = CurrentViruses(de_id(:),VirusCol.vid);
    fx=find(ismember(VirusesArray(:,1),vid_list)==1);
    VirusesArray(fx,VirusCol.death) = CTime;
    %CurrentViruses(de_id,:) = []; %shouldn't remove now. Should remove with reovery events
    %Viruses(de_id) = [];
end

function [rm_id] = recovery(CTime, rm_id)
    %for i=1:length(rm_id)
    %    vid = CurrentViruses(rm_id(i),VirusCol.vid); %convert to vid
    %    VirusesArray(vid,VirusCol.death) = CTime;
    %end
    
    CurrentViruses(rm_id,VirusCol.death) = CTime;
    vid_list = CurrentViruses(rm_id(:),VirusCol.vid);
    fx=find(ismember(VirusesArray(:,1),vid_list)==1);
    VirusesArray(fx,VirusCol.death) = CTime;
    %CurrentViruses(rm_id,:) = [];
    %Viruses(rm_id) = [];
end

function [] = removeViruses(de_id, rm_id)
    remove_set = union(de_id,rm_id);
    CurrentViruses(remove_set,:) = []; %shouldn't remove now. Should remove with reovery events
end

%Check whether any negative values in the array
function [pos]=isNegative(dBirth, dDeath_S, dInfection_tot, dInfection, dDeath_I, dRecover, dDeath_R, dWanning, dWanning_S)
  pos = -1; 
  if(checkNegative(S(1)+dBirth)==1)
    pos = 1;
    return;
  end
  
  dS = S-dDeath_S-dInfection+dWanning_S;
  if(checkNegative(dS)==1)
    pos = 1;
    return;
  end

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
    out_dir = ['out/' datestr(now,10) datestr(now,5) datestr(now,7) '_' params.proj];
    out_dir = params.out_dir
    if(exist(out_dir)==7)
    else
        mkdir(out_dir)
    end
        dat_sir = [];
        dat_sir(:,:) = Data(:,:);
        save([out_dir '/DataTLSIR.mat'],'dat_sir','params');
        disp('done');
        clear dat_sir;
end

function output_traits(CurrTime)
    % create output directory
    out_dir = ['out/' datestr(now,10) datestr(now,5) datestr(now,7) '_' params.proj];
    out_dir = params.out_dir
    if(exist(out_dir)==7)
    else
        mkdir(out_dir)
    end
    fx = find(VirusesArray(:,VirusCol.vid)==0);
    dat_VirusesArray = VirusesArray(1:fx(1)-1,:);
    save([out_dir '/virus_traits.mat'],'dat_VirusesArray');
end


function output_viruses2(VsArray)
    % create output directory
    out_dir = ['out/' datestr(now,10) datestr(now,5) datestr(now,7) '_' params.proj];
    out_dir = params.out_dir
    if(exist(out_dir)==7)
    else
        mkdir(out_dir)
    end
        dat_viruses = [];
        fx = find(VsArray(:,VirusCol.vid)==0);
        dat_viruses = [VsArray(1:fx(1)-1,VirusCol.vid), VsArray(1:fx(1)-1,VirusCol.birth), VsArray(1:fx(1)-1,VirusCol.death), VsArray(1:fx(1)-1,VirusCol.parent), VsArray(1:fx(1)-1,VirusCol.infectionK)];
        save([out_dir '/dat_x_trans_tmp.mat'],'dat_viruses');
        disp('done');
        clear dat_viruses;
end

%get estimate deltaV from a matrix
function [deltaV] = getDeltaV(vini,k)
   vlist = 0:1:500;
   vini_r=round(vini*100)+1;
   for i=1:length(vini_r)
   if isempty(find(vlist==vini_r(i)))
     vini_r(i)
     disp 'error during calculating deltaV';
   end
   end
   deltaV=diag(deltaVMatrix(k,vini_r));
end

%Plot population of different viruses 
function plot_binding_avidity(x1,t)
    figure;
    hold;
    meancharge_before = calculate_mean_proportion(x1);
    plot(t, meancharge_before);
end    
end






%%
