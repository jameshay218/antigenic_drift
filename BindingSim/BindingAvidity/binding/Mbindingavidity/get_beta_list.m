%Aug 09, 2012
%The program calculate the changes of transmission per unit change in binding avidity 
%input: k, the times of infection
%       v, the binding avidity
%output: dB(k,v)/dV

function  [val] = get_beta_list(sk, v, p, r, a, b, c, n0)
   %p, r: parameters for transmission benefit
   %a, b: parameters for replication cost
   %c: contact rate
   %no = average number of RNA copies per each virion 
   if ~exist('p')
    p = 2;
   end
   if ~exist('r')
    r = 1;
   end
   %r = 2;
   if ~exist('a')
     a = 0.7;
   end
   if ~exist('b')
     b = 3;
   end
   if ~exist('c')
     c = 0.5; %contact rate
   end
   for k=0:sk
     P_Ab = exp(-p*(v+1));        % Probability to be recognized by antibodies
    %P_Trans = (1-P_Ab).^(r*k);   % Probability to escape immunity
     P_Trans = (1-P_Ab).^(r*k);   % Probability to escape immunity
     P_Rep = exp(-a*v.^b);        % Probability to replicate
     R0 = (P_Trans.*P_Rep).*n0;
     rho = 1 - R0.^-1;
     loc_r = find(rho<=0);
     rho(loc_r) = 0;
     %val(k+1,1:length(v)) = c.*(P_Trans.*P_Rep);
     val(k+1,1:length(v)) = c.*rho;
   end   
end
   

