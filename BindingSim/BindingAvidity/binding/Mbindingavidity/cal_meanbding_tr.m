%calculate mean binding at transmission time
%output = 'out/dscr_out/20130829_lar_fixb_c08_kc05/';
global metadata;
output = [metadata.ibms.out_dir '/'];
load([output 'virus_traits.mat']);
VirusesArray = dat_VirusesArray;
interval = 60;
disp (['calculate mean binding avidity from ' output]);
mean_sampled_binding = [];
global params;
params.kc = 0.5;
params.p = 4; 
params.r = 70;
params.a = 0.7;
params.b = 3;
params.c = 0.5;

for d=1:interval:VirusesArray(end,2)
  d
  %if d == round(VirusesArray(end,2))
  %      break;
  %end
  TF = find(VirusesArray(:,2)<=d & VirusesArray(:,3)>d); %current active viruses

  indiv_sampled_infectionK = VirusesArray(TF,5);
  indiv_sampled_binding_s = VirusesArray(TF,7);
  indiv_sampled_infecteddays = d-VirusesArray(TF,2);
  
  %meanBinding(d,1) = mean(VirusesArray(TF,7));
  %meanBindingFinal(d,1) = mean(VirusesArray(TF,8));
  n_tot_samples = length(TF);
  indiv_sampled_binding = [];
  for i = 1:n_tot_samples
    %indiv_sampled_binding(i) = 0.5; %Just testing
    indiv_sampled_binding(i) = getVChange_ode(indiv_sampled_binding_s(i),indiv_sampled_infectionK(i)-1,indiv_sampled_infecteddays(i)); %getVChange_ode(initialV,infectionK,time_step)
  end
  mean_sampled_binding(end+1,1) = mean(indiv_sampled_binding);
end

save([output 'mean_binding_sampled.mat'], 'mean_sampled_binding');

