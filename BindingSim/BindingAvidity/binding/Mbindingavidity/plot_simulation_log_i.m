function [] = plot_simulation_log_i(input_dir, sy)
% plot the log of viruses simulation
addpath(genpath('..'));
%addpath(genpath('../../lib'));
%ranges = [
%          1993.601 1994.6; 
%          1994.601 1995.6;
%          ]
%ranges2 = [1993.601 2005.6];
%ranges2 = [0 35*365];   
%sample_year = 28.*365;
yrs_e = 45;
%yrs_e = 30;
ranges2 = [0 yrs_e*365];   
sample_year = sy.*365;
%sample_year = 1.*365;
%dat = load('out/20130308_2/20130314/DataTLSIR.mat'); %(thesis version)
%dat1 = load('out/20130308_2/20130314/virus_traits.mat'); %(thesis version)
%proj = '20130325';%Thesis version
%proj = '20130513';%EEID meeting
%proj = '20130630';%Testing
%proj = '20130718';%Testing
%proj = '20130801';
%dat = load(['out/' proj '/DataTLSIR.mat']);
%dat1 = load(['out/' proj '/virus_traits.mat']);

% from dscr_out
%proj = '20130812_med';
%proj = '20130814_med_fixb';
%proj = '20130815_med_n8';
%proj = '20130816_med_n4_c09';
%proj = '20130816_med_n8_c08';
%proj = '20130817_lar_n4_c07';
%proj = '20130819_lar_fixb28';
%proj = '20130827_lar_n4_c08_N5';
%proj = '20130828_lar_n4_c08_kc01';
%proj = '20130827_lar_fixb_c05';
%proj = '20130828_lar_fixb_c05_kc01';
proj = 'std_adapt_n4_c15_wan3_meds';
x1 = [];
if exist(input_dir)
    dat = load([input_dir '/DataTLSIR.mat']);
    dat1 = load([input_dir '/virus_traits.mat']);
    x1 = dat.dat_sir;
else
    dat = load(['out/' proj '/DataTLSIR.mat']);
    dat1 = load(['out/' proj '/virus_traits.mat']);
    x1 = dat.dat_sir;
end
N = sum(x1(1,2:301),2);
s = x1(:,2:101);
i = x1(:,102:201);
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
    dat1 = load(['out/' proj '/mean_binding.mat']);
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
xlabel('time(yrs)');
set(get(AX(1),'Ylabel'),'String','I_t_o_t/N'); 
%set(get(AX(2),'Ylabel'),'String','I_t_o_t/N');
ha=findobj(gcf,'type','axes');
%set(ha(1),'ylim',[0 60],'ytick',[20 40]);
set(ha(1),'ylim',[0 1],'ytick',[0 1]);

%set(AX,'XTick',0:365*5:365*yrs_e);
set(ha(1),'xlim',[0 365*yrs_e],'xtick',[0 365*yrs_e]);
set(ha(2),'xlim',[0 365*yrs_e],'xtick',[0 365*yrs_e]);
set(AX,'XTickLabel',{'0','45'});
%set(AX,'XTickLabel',{'0','5','10','15','20','25','30','35','40','45'});
%set(AX,'XTickLabel',{'0','5','10','15','20','25','30'});
%xlim([0 365*yrs_e]);
%text(100,100,'a');
%mTextBox = uicontrol('style','text')
%set(mTextBox,'String','a',BackgroundColor',[1 1 1]);
%set(mTextBox,'Position',[100 100 1 1]);

%The histrogram of Sk at year
subplot(2,2,2);
%bar(sample_s./sum(sample_s));
bar(sample_i./sum(sample_i));
xlabel('no. previous infection k');
ylabel('I_K/I_t_o_l');
set(gca,'XTick',[1:10:21]);
set(gca,'XTickLabel',{'0','10','20'});
xlim([0.5 25]);
%text(200,100,'b');
%mTextBox = uicontrol('style','text')
%set(mTextBox,'String','b',BackgroundColor',[1 1 1]);
%set(mTextBox,'Position',[200 100 1 1]);

%Mean binding avidity
subplot(2,2,3);
plot(meanBinding);
xlabel('time(yrs)');
ylabel('Average binding avidity (V)');
set(gca,'XTick',0:365*10:365*yrs_e-1);
set(gca,'XTickLabel',{'0','10','20','30','40','45'});
%set(gca,'XTickLabel',{'0','5','10','15','20','25','30','35','40','45'});
xlim([0 365*yrs_e]);
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
save([input_dir '/log.mat'], 'x1', 'sum_i_fraction','x1','sum_s_fraction','sample_i','sample_v')
end