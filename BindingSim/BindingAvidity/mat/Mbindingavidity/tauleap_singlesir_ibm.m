function [DataXt VirusesArray] = tauleap_singlesir_ibm(EndTime, DataX, time_step, meta)
%Tauleap method for single strain SIR reinfection model simulation
%This function tracks the evolution of receptor binding of viruses strains 
%in a population described by an SIR reinfection model 
%Refer to section 2A in the paper (RoySoc2013)
%Input: [EndTime, DataX, time_step, meta]
%   EndTime: integer representing the final time point and therefore time
%   period of the simulation
%   DataX: initial conditions for the SIR compartments
%   time_step: total number of steps within given time period
%   meta: stores input file names
%Output: (DataXt VirusesArray)
%   DataX: array containing time series states of SIR compartment (Xt)
%   VirusesArray: array containing traits for all active and historical viruses

%Clean version for bugs free
%Written on Mar 22, 2013
%Updated Jul 01, 2015 (JH)


%List of Events:
%Birth:         uN,             S=S+1  
%Transmission:  beta(S)(I),     I=I+1, S=S-1
%Recovery:      g(I), I=I-1,    S=S+1
%Death:         -uS,-uI,        S=S-1, I=I-1
%Wanning:       wR,             R=R-1, S=S+1

%%%------------------------------------------------------------------------
%%%--------------------StepA.1. Declare global variables---------------------
%%%------------------------------------------------------------------------
global params;      %Store parameters
global metadata;    %Store input files names
global VirusesArray;    %Store all historical viruses in an array
global CurrentViruses;  %Store current virus strains in array
maxvid = [0];       %Tracking the maximum virus ID
metadata = meta;
VirusesArray = zeros(10E6,8); %preallocate the memory for the virus array

%Load deltaV Matrix (for deltaV, please refer to eq2.8 in RoySci)
dat = open(metadata.ibms.deltaV);
%VirusCol: column names of Viruses
VirusCol = struct('vid',1,'birth',2,'death',3,'parent',4,'infectionK',5,'beta',6,'initialV',7,'currentV',8);


%%%------------------------------------------------------------------------
%%%--------------------StepA.2. Declare local variables----------------------
%%%------------------------------------------------------------------------

%number of Sk, Ik, and Rk
S = [];%%Create Susceptible groups from S(0)-S(k-1) where k = N_Infects (max number of prev. infections)
I = [];%%Create Infected groups from I(1)-I(k)
R = [];%%Create Recovered groups from R(1)-R(k)
S = DataX(1,1:params.N_Infect)';   
I = DataX(1,params.N_Infect+1:params.N_Infect*2)';  %N_Infects immune groups
R = DataX(1,params.N_Infect*2+1:params.N_Infect*3)'; 
Sk = [0:params.N_Infect-1];


%%%------------------------------------------------------------------------
%%---------------------StepA.3. Initialize virus strains---------------------
%%%------------------------------------------------------------------------

%Initialize virus strains
CurrentViruses = []; %current viral strains
if metadata.ibms.initVirusFlag == true

for i=1:length(I) %For each infected compartment
  for j=1:I(i) %For each infected individual in this compartment, create one new virus
    if exist('V_bding')
        v = V_bding;
    else
        v=0.5;
    end
    vid = maxvid+1;
    
    %calculate transmission rate based on f(k,V), g(V) and within host R0
    beta = get_beta(Sk, v, params.p, params.r, params.a, params.b, params.c, params.nv);
    
    %traits for initial viruses
    birth = 0;
    death = 0;
    parent = 0;
    infectionK = i;
    beta = beta(i); 
    initialV = v;
    currentV = v;
    
    %store initial virus traits into the virus array
    CurrentViruses(vid,:) = [vid, birth, death, parent, infectionK, beta, initialV, currentV];
    VirusesArray(vid,:) = [vid, birth, death, parent, infectionK, beta, initialV, currentV];
    maxvid = maxvid + 1;
  end
end



%%%------------------------------------------------------------------------
%%----------------------StepB. Simulation of the events--------------------
%%%------------------------------------------------------------------------
PrevTime = 0;
CurrTime = 1;
%First column -> time
DataS = [CurrTime S'];
DataI = [CurrTime I'];
DataR = [CurrTime R'];
DataXt = [DataS DataI(2:params.N_Infect+1) DataR(2:params.N_Infect+1)];

while 1
    CurrTime;
    % Single iteration of tau leap algorithm; updating state values
    % after one time step
    time_step1 = Iterate(params);
    
    % Number of current viruses
    vpopsize = length(CurrentViruses(:,VirusCol.vid));

    %stop simulation when viruses go to extinction 
    if vpopsize < 2
        % Save final states and virus array
        output_virus_transmission(VirusesArray);
        output_sir(DataXt);
        output_traits(CurrTime);
        VirusesArray(find(VirusesArray(:,VirusCol.vid)==0),:) = [];
        disp('Viruses extinction');
        disp "export"
        return;
    end
  
    %update states Xt using DataS, DataI and DataR matrix
    DataS(end+1,1) = CurrTime;
    DataS(end,2:params.N_Infect+1) = S';
    DataI(end+1,1) = CurrTime;
    DataI(end,2:params.N_Infect+1) = I';
    DataR(end+1,1) = CurrTime;
    DataR(end,2:params.N_Infect+1) = R';
    DataXt(end+1,1) = CurrTime;
    DataXt(end,2:params.N_Infect*3+1) = [S' I' R'];
    
    % Print log every 20 time steps
    if CurrTime>1 & rem(CurrTime,20) == 0
        CurrTime
        vpopsize
    end
 
    %At end of simulation time frame, stop simulation and output results
    if CurrTime >= EndTime
        % Save states and virus array
        output_virus_transmission(VirusesArray);
        output_sir(DataXt);
        output_traits(CurrTime);
        VirusesArray(find(VirusesArray(:,VirusCol.vid)==0),:) = [];
        disp "export"
        
        % Clear temporary compartment states
        DataXt = [];
        DataS = [];
        DataI = [];
        DataR = [];
        return;
    end
end
end
%%%------------------------End of Simulation-------------------------------




%%%------------------------------------------------------------------------
%%----------------------Iteration of tao leap Algorithm--------------------
%%%------------------------------------------------------------------------

function [time_step1] = Iterate(Parameters)
time_step1=time_step;

%% Update binding avidity V
vcurr = [];                                                %binding v of current viruses
vid_list = CurrentViruses(:,VirusCol.vid);                 %viruses id
vini_list = VirusesArray(vid_list,VirusCol.currentV);      %virus binding
k_list = VirusesArray(vid_list,VirusCol.infectionK);       %immune status
%calculate virus binding after one timestep using precalculated matrix
%1)use precalculated array %Currently disabled
%vcurr_list = vini_list + getDeltaV(vini_list, k_list).*time_step;  
%2)use ode simulation
vcurr_list = vini_list + getDeltaV_ode(vini_list,k_list,time_step1);%use
VirusesArray(vid_list(:),VirusCol.currentV) = vcurr_list(:,1);            
CurrentViruses(:,VirusCol.currentV) = vcurr_list(:,1);

%% Define rate of each epidemiological events
Ntot = sum(S)+sum(I)+sum(R);
Rate_Birth = params.mu*Ntot;
Rate_Death_S = params.mu*S; 
Rate_Death_I = params.mu*I;
Rate_Death_R = params.mu*R;
Rate_Recovery = params.gamma*I;
Rate_Wanning = params.wan*R;
I_total = sum(I);

%% Calculate number of events happened
%% Birth
dBirth = poissrnd(Rate_Birth*time_step);
%% Death
%Death of Susceptible
dDeath_S = poissrnd(Rate_Death_S*time_step);
%Death of Recovered
dDeath_R = poissrnd(Rate_Death_R*time_step);
%Death of Infecteds 
dDeath_I = zeros(params.N_Infect,1);

%use multinomial to prevent both 'death' and 'recovery' events occur for the same virus
%Update death and recovery at the SAME time after update transmission events
prob = [params.mu*time_step,params.gamma*time_step,1-params.mu*time_step-params.gamma*time_step];
I_loss = mnrnd(1,prob,length(CurrentViruses(:,1)));
I_indv = I_loss(:,1);       %death
R_indv = I_loss(:,2);       %recovery
de_id = find(I_indv==1);    %index of viruses among host died   
d_list = CurrentViruses(de_id, VirusCol.infectionK); %return a list of k
d_uni_list = unique(d_list);
for i=1:length(d_uni_list)
    k = d_uni_list(i);
    dDeath_I(k) = dDeath_I(k)+length(d_list(d_list==k));
end
%% Recovery
dRecover = zeros(params.N_Infect,1);
rm_id = find(R_indv==1);
k_list = CurrentViruses(rm_id, VirusCol.infectionK); %return a list of k
k_uni_list = unique(k_list);
for i=1:length(k_uni_list)
    k = k_uni_list(i);
    dRecover(k) = dRecover(k)+length(k_list(k_list==k));
end
%% Infection
beta = get_beta_list(params.N_Infect-1,vcurr_list', params.p, params.r, params.a, params.b, params.c, params.nv);
Beta = beta; %Just need to retrieve the beta matrix from the codes before
Rate_Infection = repmat(S,[1 length(Beta(1,:))]).*Beta/params.N;
dInfection_indv = poissrnd(Rate_Infection*time_step);
dInfection = sum(dInfection_indv,2); 
dInfection_tot = sum(dInfection);
%% Wanning
dWanning = poissrnd(Rate_Wanning*time_step);
dWanning_S = [0; dWanning(1:params.N_Infect-1)];

% Check whether any negative values would be generated before updating the event 
if isNegative(dBirth, dDeath_S, dInfection_tot, dInfection, dDeath_I, dRecover, dDeath_R, dWanning, dWanning_S) == 1
        %time_step = time_step/2;
        return; %stop updating and return to start if there is negative values 
end

%% Updating the epidemiological status Xt
PrevTime = CurrTime;
CurrTime = PrevTime + time_step; %Tracking the current time
S(1)=S(1)+dBirth; %1) Birth
S=S-dDeath_S;     %2) Death from susceptible individuals
I=I-dDeath_I;     %3) Death from infected individuals
R=R-dDeath_R;     %4) Death from recovered individuals
I=I+dInfection;   %5) New Infectious individuals
S=S-dInfection;   %6) Loss of susceptible individuals from Infection
I=I-dRecover;     %7) Loss of infected individuals from Recovery
R=R+dRecover;     %8) New recovered individuals
R=R-dWanning;     %9) Loss if recovered individuals after immune wanning 
S=S+dWanning_S;   %10)Susceptible individuals from immune waning

%% Updating transmission tree for each viral strain in CurrentViruses 
SourceViruses = CurrentViruses(:, VirusCol.vid);        %the source viruses
number_new_viruses = sum(dInfection_indv);              %[Sk x Vi] number of the new viruses offsprings from each source viruses
ns_list = find(number_new_viruses>0);                   %the index of the source viruses who produce offsprings>0
%only search for viruses who produce offsprings
for i =1:length(ns_list) %run for each source viruses (i)    
        ns = ns_list(i); %source viruses index
        vcurr = CurrentViruses(ns,VirusCol.currentV);   %source binding avidity
        number_of_infected_by_virus = sum(dInfection_indv(:,ns));
        k_list1 = find(dInfection_indv(:,ns)>0);        %get the list of [Sk] that are infected
        vid_array = [];
        k_array = [];
        beta_array = [];
        for i=1:length(k_list1)     %run for each Sk infected by virus (i)
            number_of_inf_Sk = dInfection_indv(k_list1(i),ns); %#of Sk
            vid_array(end+1:end+number_of_inf_Sk) = maxvid+1:maxvid+number_of_inf_Sk; %vid of newly infected
            k_array(end+1:end+number_of_inf_Sk) = k_list1(i).*ones(1,number_of_inf_Sk);%k array for newly infected
            beta_array(end+1:end+number_of_inf_Sk) = Beta(k_list1(i),ns).*ones(1,number_of_inf_Sk);%susceptibility of each new patient
            maxvid = maxvid+number_of_inf_Sk;                                           %update maxvid
        end
        %same as the columns of the CurrentViruses 
        vid_j = vid_array;
        birth_j = CurrTime;
        death_j = 0;
        parent_j = SourceViruses(ns);
        infectionK_j = k_array;
        beta_j = beta_array;
        initialV_j = vcurr; 
        currentV_j = vcurr;     
        if number_of_infected_by_virus>1
          number_of_infected_by_virus;
        end
        %Insert offspring viruses into CurrentViruses 
        CurrentViruses(end+1:end+number_of_infected_by_virus,:) = [vid_j', birth_j.*ones(1,number_of_infected_by_virus)', death_j.*ones(1,number_of_infected_by_virus)', parent_j.*ones(1,number_of_infected_by_virus)', infectionK_j', beta_j', initialV_j.*ones(1,number_of_infected_by_virus)', currentV_j.*ones(1,number_of_infected_by_virus)'];
        %Insert into VirusesArray
        VirusesArray(vid_j,:) = [vid_j', birth_j.*ones(1,number_of_infected_by_virus)', death_j.*ones(1,number_of_infected_by_virus)', parent_j.*ones(1,number_of_infected_by_virus)', infectionK_j', beta_j', initialV_j.*ones(1,number_of_infected_by_virus)', currentV_j.*ones(1,number_of_infected_by_virus)']; %insert virus id into CurrentViruses (vid, infecitonK, beta)
end
     
%% Removing viruses from CurrentViruses after recovery or death
de_id1 = deathOfInfectds(PrevTime, de_id); %death_ID, determines which infected individuals are going to be died
rm_id1 = recovery(PrevTime, rm_id);        %recovered_ID, determines which infected individuals are going to be recovered
%Clear the viruses
removeViruses(de_id, rm_id);
time_step1 = time_step;
end
%%---end of iteration---###

%%%------------------------------------------------------------------------
%%--------------------------Auxiliary functions----------------------------
%% Sub functions to remove viruses
function [de_id] = deathOfInfectds(CTime, de_id)
    CurrentViruses(de_id,VirusCol.death) = CTime;
    vid_list = CurrentViruses(de_id(:),VirusCol.vid);
    fx=find(ismember(VirusesArray(:,1),vid_list)==1);
    VirusesArray(fx,VirusCol.death) = CTime;
end

function [rm_id] = recovery(CTime, rm_id)
    CurrentViruses(rm_id,VirusCol.death) = CTime;
    vid_list = CurrentViruses(rm_id(:),VirusCol.vid);
    fx=find(ismember(VirusesArray(:,1),vid_list)==1);
    VirusesArray(fx,VirusCol.death) = CTime;
end

function [] = removeViruses(de_id, rm_id)
    remove_set = union(de_id,rm_id);
    CurrentViruses(remove_set,:) = [];
end

%% Check whether any negative values in the array
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

%% Sub functions for saving output status Xt
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

%% Sub functions for saving output Virus traits
function output_traits(CurrTime)
    % create output directory
    out_dir = ['out/' datestr(now,10) datestr(now,5) datestr(now,7) '_' params.proj];
    out_dir = params.out_dir
    if(exist(out_dir)==7)
    else
        mkdir(out_dir)
    end
    fx = find(VirusesArray(:,VirusCol.vid)==0,1); % Don't retrieve VirusCol.vid = 0.
    if isempty(fx)
          fx = length(VirusesArray(:,1)) + 1;
    end
    dat_VirusesArray = VirusesArray(1:fx-1,:);
    save([out_dir '/virus_traits.mat'],'dat_VirusesArray');
end

%% Sub functions for saving output Virus Transmission Arrays
function output_virus_transmission(VsArray)
    % create output directory
    out_dir = ['out/' datestr(now,10) datestr(now,5) datestr(now,7) '_' params.proj];
    out_dir = params.out_dir
    if(exist(out_dir)==7)
    else
        mkdir(out_dir)
    end
        dat_viruses = [];
        fx = find(VsArray(:,VirusCol.vid)==0,1);
        if isempty(fx)
          fx = length(VsArray(:,1)) + 1;
        end
        dat_viruses = [VsArray(1:fx-1,VirusCol.vid), VsArray(1:fx-1,VirusCol.birth), VsArray(1:fx-1,VirusCol.death), VsArray(1:fx-1,VirusCol.parent), VsArray(1:fx-1,VirusCol.infectionK)];
        
        save([out_dir '/dat_x_trans_tmp.mat'],'dat_viruses');
        disp('done');
        clear dat_viruses;
end

%% Other Sub functions might be used in the future
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




