
p = 1;
k = 1;
g = @(V,p,k,a,b) exp(-(V+1))*(p*k*1+a*b*V^(b-1))-a*b*V.^(b-1);

fzero(g,0,[],1,1,0.7,3)