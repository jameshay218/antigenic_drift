function [ X1_2D ] = getNewImmuneStatusLinear2D( k, mu)
%immune status adjusted by mutation
%   Detailed explanation goes here

        matK = repmat(k,length(mu),1);
        matM = repmat(mu,1,length(k));
        X1_2D = matK - matM;
end

