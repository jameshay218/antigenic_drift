function [ X1 ] = getNewImmuneStatus( k, mu, p, v, r )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


   P_Ab = exp(-p*(v+1));         % Probability to be recognized by antibodies
   P_Esc_Boost = 1-P_Ab;               % Probability of immune escape per each boosting
   
   if size(k,2)==1
        P_Esc = (P_Esc_Boost).^(r*k);   % Probability to escape immunity for given k
        X1= round((log(P_Esc + (1-P_Esc).*(1-exp(-mu)))./log(P_Esc_Boost))./r);
   else
   %Probability to escape immunity for all given k
   for k=0:99
    P_Esc(:,k+1) = (P_Esc_Boost).^(r*k);
   end
   log_esc_boost_array = repmat(log(P_Esc_Boost),[1 100]);
   X1= round((log(P_Esc + (1-P_Esc).*(1-repmat(exp(-mu),[1 100])))./log_esc_boost_array)./r);
   end
   %X1 = floor((log(P_Esc + (1-P_Esc).*(1-exp(-mu)).*0)./log(P_Esc_Boost))./r)+1;
   %X1 = floor((log(P_Esc)./log(P_Esc_Boost))./r);
   X1(find(X1<1))=1; 
end

