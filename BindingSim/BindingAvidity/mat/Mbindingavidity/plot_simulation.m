function [] = plot_simulation_log_i(input_dir)
% plot the log of viruses simulation


proj = '';

x1 = [];
    dat = load([input_dir '/DataTLSIR.mat']);
    dat1 = load([input_dir '/virus_traits.mat']);
    x1 = dat.dat_sir;

    end_date = length(x1(:,1));
ranges2 = [0 end_date] %starting date to end date
sample_year = end_date-30;

params = dat.params;
N = sum(x1(1,2:end),2);
s = x1(:,2:1+params.N_Infect);
i = x1(:,1+params.N_Infect+1:1+params.N_Infect*2);
sum_s = sum(s,2);
sum_s_fraction = sum_s./N;
sum_i = sum(i,2);
sum_i_fraction = sum_i./N;
TF = find(x1(:,1)>=sample_year & x1(:,1)<sample_year+1);
sample_s = sum(s(TF,:),1);
sample_i = sum(i(TF,:),1);
VirusesArray = dat1.dat_VirusesArray;
TF1 = find(VirusesArray(:,2)<=sample_year & VirusesArray(:,3)>sample_year);
sample_v = VirusesArray(TF1,7); %binding avidity at transmission time
%sample_v = VirusesArray(TF1,8); %binding avidity before final clearance
meanBinding = [];
if exist(input_dir)
    dat1 = load([input_dir '/mean_binding.mat']);
    meanBinding = dat1.meanBinding;
else
    dat1 = load(['out/dscr/' proj '/mean_binding.mat']);
    meanBinding = dat1.meanBinding;
end
lastday = x1(end,1);
%dat1.meanBinding
%
h = figure;
%The fraction of Sk/N
subplot(2,2,1);
%plot(x1(:,1),sum_s_fraction);hold on;
%plot(x1(:,1),sum_i_fraction*50,'r');
%xlabel('time(yrs)');
%ylabel('S_t_o_t/N');
%set(gca,'XTick',0:365*5:365*yrs_e);
%set(gca,'XTickLabel',{'0','5','10','15','20','25','30','35','40','45'});
%xlim([0 365*yrs_e]);

[AX,H1,H2] = plotyy(x1(:,1),sum_i_fraction,x1(:,1),sum_s_fraction,'plot');
xlabel('time(days)');
set(get(AX(1),'Ylabel'),'String','I_t_o_t/N'); 
%set(get(AX(2),'Ylabel'),'String','I_t_o_t/N');
ha=findobj(gcf,'type','axes');
%set(ha(1),'ylim',[0 60],'ytick',[20 40]);
set(ha(1),'ylim',[0 1],'ytick',[0 1]);

%set(AX,'XTick',0:365*5:365*yrs_e);
set(ha(1),'xlim',[0 ranges2(2)],'xtick',[0 ranges2(2)]);
set(ha(2),'xlim',[0 ranges2(2)],'xtick',[0 ranges2(2)]);
set(AX,'XTickLabel',{'0', ranges2(2)});


%The histrogram of Sk at year
subplot(2,2,2);
%bar(sample_s./sum(sample_s));
bar(sample_i./sum(sample_i));
xlabel('no. previous infection k');
ylabel('I_K/I_t_o_l');
set(gca,'XTick',[1:10:21]);
set(gca,'XTickLabel',{'0','10','20'});
xlim([0.5 25]);


%Mean binding avidity
subplot(2,2,3);
plot(meanBinding);
xlabel('time(days)');
ylabel('Average binding avidity (V)');
set(gca,'XTick',[0 ranges2(2)]);
set(gca,'XTickLabel',{'0', ranges2(2)});
xlim([0 ranges2(2)]);


%%he histogram of V at year
subplot(2,2,4);
hist(sample_v, 17);
xlabel('Binding avidity');
ylabel('Frequency');
%The number of Ik
%subplot(3,2,3);
%plot(x1(:,1));
%The histrogram of Ik at year
%subplot(3,2,4);
%hist(sample_i);
hgsave(h,[input_dir '/Figure6.fig']);
print(h,'-dpng',[input_dir '/Figure6.png']);
end