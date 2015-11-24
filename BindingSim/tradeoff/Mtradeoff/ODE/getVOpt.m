function v2 = getVOpt(num)

%p=2;
p=1;
r=1;
a=0.7;
b=3;
c=0.5;

i=0:num
h=exp(p)*(((b-1)/p).^(b-1))
k2= (i/h).^(1/(b-1));
v1 = lambertw(k2);
v2 = v1*((b-1)/p);

end