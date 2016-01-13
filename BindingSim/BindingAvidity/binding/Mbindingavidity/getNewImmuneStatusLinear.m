function [ X1 ] = getNewImmuneStatusLinear( k, mu)
%immune status adjusted by mutation
%   Detailed explanation goes here
   if length(k) == length(mu)
        X1 = round(k - mu); 
        X1(find(X1<1))=1; 
   else
        matK = repmat(k,length(mu));
        matM = repmat(mu',1,length(k));
   end
end

