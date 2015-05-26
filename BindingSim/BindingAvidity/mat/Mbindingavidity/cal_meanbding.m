function [] = cal_meanbding(metadata)
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
end