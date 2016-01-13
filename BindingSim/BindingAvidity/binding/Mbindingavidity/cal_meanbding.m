%output = 'out/20130325/'; %Thesis version
%output = 'out/20130513/'; %EEID version
%output = 'out/20130630/'; %Testing
%output = 'out/20130717/'; %Testing N=2x10^6
%output = 'out/20130718/'; %Large sample size N=3x10^6. Use this as
%representative result.
%output = 'out/20130801/';
%output = 'out/dscr_out/20130829_lar_n4_c08_kc05/';
%output = 'out/dscr_out/20130829_lar_fixb_c08_kc05/';
global metadata;
output = [metadata.ibms.out_dir '/'];
load([output 'virus_traits.mat']);
VirusesArray = dat_VirusesArray;
disp (['calculate mean binding avidity from ' output]);
for d=1:round(VirusesArray(end,2))
  TF = find(VirusesArray(:,2)<=d-1 & VirusesArray(:,3)>d-1);
  meanBinding(d,1) = mean(VirusesArray(TF,7));
  meanBindingFinal(d,1) = mean(VirusesArray(TF,8));
end
save([output 'mean_binding.mat'], 'meanBinding', 'meanBindingFinal');
