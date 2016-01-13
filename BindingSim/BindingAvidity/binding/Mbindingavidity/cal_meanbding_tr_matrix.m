function [] = cal_meanbding_tr_matrix()
%calculate mean binding at transmission time
%output = 'out/dscr_out/20130829_lar_fixb_c08_kc05/';
global metadata;
output = [metadata.ibms.out_dir '/'];
dat = load([output 'virus_traits.mat']);
VirusesArray = dat.dat_VirusesArray;
dat = load('dat/deltaVMatrix_kc05.mat');
deltaVMatrix = dat.deltaVMatrix;
interval = 30;
disp (['calculate mean binding avidity from ' output]);
mean_sampled_binding = [];

for d=1:interval:VirusesArray(end,2)
  d
  %if d == round(VirusesArray(end,2))
  %      break;
  %end
  TF = find(VirusesArray(:,2)<=d & VirusesArray(:,3)>d); %current active viruses

  indiv_sampled_infectionK = VirusesArray(TF,5);
  indiv_sampled_binding_s = VirusesArray(TF,7);
  indiv_sampled_infecteddays = d-VirusesArray(TF,2);
  vcurr_list = indiv_sampled_binding_s; 
  %meanBinding(d,1) = mean(VirusesArray(TF,7));
  %meanBindingFinal(d,1) = mean(VirusesArray(TF,8));
  n_tot_samples = length(TF);
  indiv_sampled_binding = [];
  time_step = zeros(n_tot_samples,1);
  for i = 0:max(indiv_sampled_infecteddays)
    %indiv_sampled_binding(i) = 0.5; %Just testing
    %indiv_sampled_binding(i) = getVChange_ode(indiv_sampled_binding_s(i),indiv_sampled_infectionK(i)-1,indiv_sampled_infecteddays(i)); %getVChange_ode(initialV,infectionK,time_step)
    vini_list = vcurr_list;
    k_list = indiv_sampled_infectionK;
    remain_time_step = indiv_sampled_infecteddays - i;
    time_step = zeros(n_tot_samples,1);
    time_step(find(remain_time_step>0),1)=1;
    vcurr_list = vini_list + getDeltaV(vini_list, k_list).*time_step;
  end
  indiv_sampled_binding = vcurr_list;
  mean_sampled_binding(end+1,1) = mean(indiv_sampled_binding);
end

save([output 'mean_binding_sampled(matrix).mat'], 'mean_sampled_binding');

%get estimate deltaV from a matrix
function [deltaV] = getDeltaV(vini,k)
   vlist = 0:1:500;
   vini_r=round(vini*100)+1;
   %for i=1:length(vini_r)
   %if isempty(find(vlist==vini_r(i)))
   %  vini_r(i)
   %  disp 'error during calculating deltaV';
   %end
   %end
   deltaV=diag(deltaVMatrix(k,vini_r));
end
end