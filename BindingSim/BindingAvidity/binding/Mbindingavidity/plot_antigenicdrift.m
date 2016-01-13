function [ ] = plot_antigenicdrift( dat_VirusesArray )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%-- Plot antigenic drift by days
figure;
day = dat_VirusesArray(:,2);
plot(day, dat_VirusesArray(:,end),'.'); %drift
hold on;


vid = find(dat_VirusesArray(:,9)>4.5);
x1 = dat_VirusesArray(vid,2);
y1 = dat_VirusesArray(vid,11);
vid2 = dat_VirusesArray(vid,4);
rmid = find(vid2==0);
vid2(rmid) = [];
x1(rmid) = [];
y1(rmid) = [];
x2 = dat_VirusesArray(vid2,2);
y2 = dat_VirusesArray(vid2,11);
for i = 1:length(x2)
	plot([x2(i) x1(i)],[y2(i) y1(i)],'g-');
end

%%-- ByProduct
vid = dat_VirusesArray(1000:end,1);
total = length(vid);
for i = 1: 3000
vid = round(total*rand(1,1));
x1 = dat_VirusesArray(vid,2);
y1 = dat_VirusesArray(vid,10);
vid2 = dat_VirusesArray(vid,4);
if vid2==0
  continue;
end
x2 = dat_VirusesArray(vid2,2);
y2 = dat_VirusesArray(vid2,11);
y1 = y1 + y2;
bypro = abs(y2-y1);
if bypro>0.01
plot([x2 x1],[y2 y1],'r-','LineWidth',2);
end
if bypro<=0.01
plot([x2 x1],[y2 y1],'c-','LineWidth',2);
end
end


%%-- ByProduct
%virusTOT = dat_VirusesArray;
%vid = find(dat_VirusesArray(:,9)>5);
%virusATG = dat_VirusesArray(vid,:);

%for i=1:length(virusATG(:,1))
%x2 = virusATG(i,2); 
%y2 = virusATG(i,11);
%vidchildren = find(virusTOT(:,4)==virusATG(i,1));
%for j=1:length(vidchildren)
%  x1 = virusTOT(vidchildren(j),2);
%  y1 = virusTOT(vidchildren(j),11);
%  plot([x2 x1],[y2 y1],'r-');
%end
%end
figure;
hold on
virus = dat_VirusesArray;
for d = 1:5:900
        vid = find(virus(:,2)>d & virus(:,2)<d+10);
	immuneJ = virus(vid,5)-(virus(vid,11)-virus(vid,10)-virus(vid,9)); 
        immuneJ(find(immuneJ<0)) = 0;
        plot(d+10,mean(immuneJ),'k.');
        plot(d+10,mean(virus(vid,7)),'.'); 
end



end

