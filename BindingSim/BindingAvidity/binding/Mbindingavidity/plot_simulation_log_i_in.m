function [] = plot_simulation_log_i(input_dir)
% plot the log of viruses simulation

    dat = load([input_dir '/mean_binding_in.mat']);
    meanBinding = dat.meanBinding;
    x1 = dat.x1;
    sum_s_fraction = dat.sum_s_fraction;
    sum_i_fraction = dat.sum_i_fraction;
    sample_i = dat.sample_i;
    sample_v = dat.sample_v;
    lastday = x1(end,1);

    yrs_e = 90;
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
set(AX,'XTickLabel',{'0','90'});
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
set(gca,'XTick',0:365*10:365*yrs_e);
set(gca,'XTickLabel',{'0','10','20','30','40','50','60','70','80','90'});
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