%Nov 15, 2011
%The program calculate the changes of transmission per unit change in binding avidity 
%input: k, the times of infection
%       v, the binding avidity
%output: dB(k,v)/dV

function  val = dBeta_dV(k, v, p, r, a, b, c)
   %p, r: parameters for transmission benefit
   %a, b: parameters for replication cost
   %c: contact rate
   
   P_Ab = exp(-p*(v+1));        % Probability to be recognized by antibodies
   P_Trans = (1-P_Ab).^(r*k);   % Probability to escape immunity
   P_Trans_Pr = r*k*p.*((1-P_Ab).^(r*k-1)).*(P_Ab);
   P_Rep = exp(-a*v.^b);        % Probability to replicate
   P_Rep_Pr = -a*b*exp(-a*v.^b).*(v.^(b-1));
   %val = c.*(P_Trans_Pr.*P_Rep + P_Trans.*P_Rep_Pr);
   val = P_Trans_Pr.*P_Rep + P_Trans.*P_Rep_Pr; 
end
   

