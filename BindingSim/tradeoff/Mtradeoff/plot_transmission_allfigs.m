%Subplot1 f(k,V)
%Subplot2 g(V) = exp(-a*V.^b)
%Subplot3 Beta(k,V)
%Subplot4 dBeta(k,V)/dV
%Author: Hsiang-Yu Yuan
%1st version: Jul 16, 2012
%Rho = f x g
%2nd verison: Aug 01, 2013
%Rho = 1 - R0^-1
%R0 = f x g x n

%Binding avidity range
V = 0:0.01:6;

%Transmission parameters
p = 2;
r = 1;
b = 3;
a = 0.7;
c = 0.5; % contact rate
nv = 8; % average copies number of each virion

%Epidemiological parameters
mu_val = 1/(70*365.25); %lifespan = 70y
gamma_val=1/5; %infecious period should change to 5d


Trans_array = []; %Array of probability of evasion by immune system
Trans_Pr_array = [];


N_Reinfect = 20; 
%colorm = {[1 1 0];[1 0 1];[0 1 1];[1 0 0];[0 1 0];[0 0 1];[0 0 0]}; 
%
for i=1:N_Reinfect;
colorm(i) = {((N_Reinfect-i)./(N_Reinfect-1))*[0 0 1]+((i-1)./(N_Reinfect-1))*[1 0 0]}; 
end
for k = 0:N_Reinfect-1
%for k=1:20:21;   
   P_Ab = exp(-p*(V+1));
   
   %%%P_Trans: f(k,V)
   %%%P_Trans_Pr: f'(k,V)
   P_Trans = (1-P_Ab).^(r*k); 
   if k>=1
       P_Trans_Pr = r*k*p.*((1-P_Ab).^(r*k-1)).*(P_Ab);
   elseif k == 0
       P_Trans_Pr = zeros(1,length(V));
   end
   subplot(2,2,1); plot(V,P_Trans,'color',cell2mat(colorm(k+1))); hold on;
   Trans_array(k+1,:) = P_Trans;
   Trans_Pr_array(k+1,:) = P_Trans_Pr;
end
subplot(2,2,1); xlabel('Binding avidity V'); ylabel('Probability of immune escape');

Trans_array;

%%%P_Rep: g(V)
%%%P_Rep_Pr: g'(V)
P_Rep = exp(-a*V.^b);
P_Rep_Pr = -a*b*exp(-a*V.^b).*(V.^(b-1));
subplot(2,2,2); plot(V, P_Rep,'color',[0 0 0]); hold on;
%subplot(3,2,4); plot(B, P_Rep_Pr); hold on;
subplot(2,2,2); xlabel('Binding avidity V'); ylabel('Probability of succesfull replication')
B_Pr_array = [];

for i=1:length(Trans_Pr_array(:,1))
%R0 = f(k,V)g(V)n
%rho = 1 - 1/R0
%beta = c x rho
%R0_Trans = Trans_array(i,:).*P_Rep.*nv;
%Rho_Trans = 1 - R0_Trans.^-1;
%Rho_Trans(find(Rho_Trans<0))=0;
%B_Trans = c.*Rho_Trans;
%subplot(2,2,3); plot(V,Rho_Trans,'color',cell2mat(colorm(i))); hold on;

B_Trans = c.*Trans_array(i,:).*P_Rep;
%R0_Trans = c.*Trans_array(i,:).*P_Rep/(mu_val+gamma_val);
subplot(2,2,3); plot(V,B_Trans,'color',cell2mat(colorm(i))); hold on;
%B_Pr: f'(k,V)g(V) + f(k,V)g'(V);
B_Pr = Trans_Pr_array(i,:).*P_Rep + Trans_array(i,:).*P_Rep_Pr;
B_Pr_array(i,:) = B_Pr;

subplot(2,2,4); plot(V,B_Pr,'color',cell2mat(colorm(i))); hold on;
end
subplot(2,2,3); xlabel('Binding avidity V'); ylabel('Infection probability given a contact')
subplot(2,2,4); xlabel('Binding avidity V'); ylabel('d\beta/dV');
line([0; 3],[0; 0]);
B_Pr_array;
   
   
R0_S0= c.*Trans_array(1,:).*P_Rep/(mu_val+gamma_val);
R0_S0(1,1)

