function [ output_args ] = plot_prob_transmission( input_args )
%PROB_TRANSMISSION Summary of this function goes here
%   Detailed explanation goes here

b = 2;
a = 0.005;


d2 = 0:1:12;
K2 = 2.^d2;
%K1 = 0:0.01:1;
d1 = 0:0.5:30;
K1 = 2.^d1;
d3 = [15 20];
K3 = 2.^d3;
c = 10;
c3 = 1;
for i=1:length(K2)
  f(i,:) = 1./(1+c*(K2(i)./K1)); %fitness is defined by contact rate 
end

figure;
subplot(2,2,1);
hold;
for i1=1:length(f(:,1))
 plot(f(i1,:),'b');
end

subplot(2,2,2);
hold;

for i1=1:length(K3)
        for i=1:length(K1)
        g(i1,i) = 1./(1+c3*(K1(i)./K3(i1))); %fitness is defined by contact rate 
        %if i1 == 1
        %    plot(g(i,:));
        %else
        %    plot(g(i,:),'r');
        %end
    end
end

for i=1:length(g(:,1))
    if i==1
            plot(g(i,:));
        else
            plot(g(i,:),'r');
    end
end
%P_Rep = exp(-a.*d1.^b);
%P_Rep = 
subplot(2,2,3);
hold;
for i1=1:length(f(:,1))
 plot(f(i1,:).*g(1,:),'b');
end

for i=1:length(K2)
df(i,:) = (c.*K2(i).*(f(i,:).^2)).*(2.^(-d1));
end

for i1=1:length(K3)
dg(i1,:) = -(c3./K3(i1)).*(g(i1,:).^2).*(2.^d1);
end

subplot(2,2,4);
hold;
for i=1:length(K2)
dfg = df(i,:).*g(1,:)+dg(1,:).*f(i,:);
plot(dfg);
end

%for i1=1:length(f(:,1))
% plot(f(i1,:).*g(:,2)','r');
%end

%r_v=0.8
%for k=1:20
%  f = r_v^k
%  plot(k,f,'r.');
%end

end

