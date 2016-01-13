%New infection fitness model developed in 2015
%Subplot1 f(k,V)
%Subplot2 g(V) = exp(-a*V.^b)
%Subplot3 Beta(k,V)
%Subplot4 dBeta(k,V)/dV
%Author: Hsiang-Yu Yuan
%4 May 2015


%Binding avidity range
V = 0:0.005:2;

%Transmission parameters
%%p=2,r=2

p = 4; %original value is 4
%r = 1;
r = 70;
b = 3;
a = 0.7;
c = 2; % contact rate
nv = 4; % average copies number of each virion

%Epidemiological parameters
mu_val = 1/(70*365.25); %lifespan = 70y
gamma_val=1/5; %infecious period should change to 5d


Trans_array = []; %Array of probability of evasion by immune system
Trans_Pr_array = [];

N_Reinfect = 32; 

figure;
%colorm = {[1 1 0];[1 0 1];[0 1 1];[1 0 0];[0 1 0];[0 0 1];[0 0 0]}; 
%
%for i=1:N_Reinfect;
maximm_Reinfect = N_Reinfect;
for i=1:maximm_Reinfect;
colorm(i) = {((maximm_Reinfect-i)./(maximm_Reinfect-1))*[0 0 1]+((i-1)./(maximm_Reinfect-1))*[1 0 0]}; 
end
if maximm_Reinfect < N_Reinfect
for i=maximm_Reinfect+1:N_Reinfect
  colorm(i) = colorm(maximm_Reinfect);
end
end

D = 0; %antigenic distance by boosting
delta = 0; %antigenic distance by drift
for k = 0:N_Reinfect-1
%for k=1:20:21;   
   P_Ab = exp(-p*(V+1));
   P_s0 = 1-P_Ab;
   %%%P_Trans: f(k,V)
   %%%P_Trans_Pr: f'(k,V)
   

   immK = r*k-delta;
   if immK < 0
     immK = 0;
   end
   P_Trans = (P_s0).^(immK);
   if k>=1
       P_Trans_Pr = immK*p.*((1-P_Ab).^(immK-1)).*(P_Ab);
   elseif k == 0
       P_Trans_Pr = zeros(1,length(V));
   end
   
   %subplot(2,2,1); plot(V,P_Trans,'color',cell2mat(colorm(k+1))); hold on;
   Trans_array(k+1,:) = P_Trans;
   Trans_Pr_array(k+1,:) = P_Trans_Pr;
end
%subplot(2,2,1); xlabel('Binding avidity V'); ylabel('Probability of immune escape');

Trans_array;

%%%P_Rep: g(V)
%%%P_Rep_Pr: g'(V)
P_Rep = exp(-a*V.^b);
P_Rep_Pr = -a*b*exp(-a*V.^b).*(V.^(b-1));

%subplot(2,2,2); plot(V, P_Rep,'color',[0 0 0]); hold on;
%subplot(2,2,2); xlabel('Binding avidity V'); ylabel('Probability of succesfull replication')




for i=1:length(Trans_Pr_array(:,1))
%R0 = f(k,V)g(V)n
%rho = 1 - 1/R0
%beta = c x rho
R0_Trans = Trans_array(i,:).*P_Rep.*nv;
Rho_Trans = 1 - R0_Trans.^-1;
Rho_Trans(find(Rho_Trans<0))=0;
Rho_Trans_array(i,:) = Rho_Trans;
%%Calculate the derivative of phi
B_Pr = Trans_Pr_array(i,:).*P_Rep + Trans_array(i,:).*P_Rep_Pr;
B_Pr_array(i,:) = B_Pr;

%B_Trans = c.*Rho_Trans;

subplot(1,2,2);
plot(V,Rho_Trans,'color',cell2mat(colorm(i)),'LineWidth',1.2); hold on;
%subplot(2,2,4); plot(V,B_Pr,'color',cell2mat(colorm(i))); hold on;


B_Trans = Trans_array(i,:).*P_Rep;
%R0_Trans = c.*Trans_array(i,:).*P_Rep/(mu_val+gamma_val);
subplot(1,2,1); plot(V,B_Trans,'color',cell2mat(colorm(i)),'LineWidth',1.2); hold on;
end

%%Calculate the max rho
rho0 = max(Rho_Trans_array(1,:));
rho1 = max(Rho_Trans_array(2,:));
rho2 = max(Rho_Trans_array(3,:));

subplot(1,2,1);
xlabel('binding avidity'); ylabel('probability of replication within a host (R_i_n)');
subplot(1,2,2);
xlabel('binding avidity'); ylabel('probability of infection between the hosts');
line([0 2], [rho0 rho0], 'LineWidth',1.2, 'LineStyle','--');
line([0 2], [rho1 rho1], 'LineWidth',1.2, 'LineStyle','--');
line([0 2], [rho2 rho2], 'LineWidth',1.2, 'LineStyle','--');
%subplot(2,2,4); xlabel('Binding avidity V'); ylabel('d\beta/dV');
   
figure;
hold on;
for i=1:length(Trans_Pr_array(:,1))
    plot(V,B_Pr_array(i,:),'color',cell2mat(colorm(i)));
end

