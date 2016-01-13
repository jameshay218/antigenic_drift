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
%yrs_e = 45;
yrs_e = 2;
ranges2 = [0 yrs_e*365];   
%sample_year = 40.*365;
sample_year = 1.5.*365;
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

%optimum
%proj = '20130827_lar_n4_c08_N5';
proj = '20151120/sml';

%kc = 0.1
%proj = '20130828_lar_n4_c08_kc01';
%proj = '20130828_lar_fixb_c05_kc01';

dat = load(['out/' proj '/DataTLSIR.mat']);
dat1 = load(['out/' proj '/virus_traits.mat']);
x1 = dat.dat_sir;
N = sum(x1(1,2:601),2);
s = x1(:,2:201);
i = x1(:,102:401);
sum_s = sum(s,2);
sum_s_fraction = sum_s./N;
sum_i = sum(i,2);
sum_i_fraction = sum_i./N;
TF = find(x1(:,1)>=sample_year & x1(:,1)<sample_year+1);
sample_s = sum(s(TF,:),1);
sample_i = sum(i(TF,:),1);
VirusesArray = dat1.dat_VirusesArray;
TF1 = find(VirusesArray(:,2)<=sample_year & VirusesArray(:,3)>sample_year);
sample_v = VirusesArray(TF1,7); %initial V
%dat1 = load(['out/' proj '/mean_binding.mat']);
%dat1 = load(['dscr_out/' proj '/mean_binding.mat']);
%meanBinding = dat1.meanBinding;
%lastday = x1(end,1);

figure;
%The fraction of Sk/N
subplot(2,2,1);
[AX,H1,H2] = plotyy(x1(:,1),sum_s_fraction,x1(:,1),sum_i_fraction,'plot');
xlabel('time(yrs)');
set(get(AX(1),'Ylabel'),'String','S_t_o_t/N'); 
%set(get(AX(2),'Ylabel'),'String','I_t_o_t/N');
set(AX,'XTick',0:100:2*365);


%The histrogram of Sk at year
subplot(2,2,2);
%bar(sample_s./sum(sample_s));
bar(sample_i./sum(sample_i));
xlabel('no. previous infection k');
ylabel('S_K/S_t_o_l');
set(gca,'XTick',[1:20:61]);
set(gca,'XTickLabel',{'0','20','40','60'});
xlim([0.5 61]);
%text(200,100,'b');
%mTextBox = uicontrol('style','text')
%set(mTextBox,'String','b',BackgroundColor',[1 1 1]);
%set(mTextBox,'Position',[200 100 1 1]);
