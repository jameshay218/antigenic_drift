
function [] = merge_meanbding(input_dir1, input_dir2, target_dir)
input1 = input_dir1;
input2 = input_dir2;

dat = load([input1 '/mean_binding_sampled(matrix).mat']);
meanBindingTrans1 = dat.mean_sampled_binding;
dat = load([input2 '/mean_binding_sampled(matrix).mat']);
meanBindingTrans2 = dat.mean_sampled_binding;
meanBindingTrans = [meanBindingTrans1; meanBindingTrans2];

dat = load([input1 '/mean_binding.mat']);
meanBindingFinal1 = dat.meanBindingFinal;
dat = load([input2 '/mean_binding.mat']);
meanBindingFinal2 = dat.meanBindingFinal;
meanBindingFinal = [meanBindingFinal1; meanBindingFinal2];


dat = load([input1 '/mean_binding.mat']);
meanBinding1 = dat.meanBinding;
dat = load([input2 '/mean_binding.mat']);
meanBinding2 = dat.meanBinding;
meanBinding = [meanBinding1; meanBinding2];



dat = load([input1 '/log.mat']);
x1_1 = dat.x1;
sum_i_fraction_1 = dat.sum_i_fraction;
sum_s_fraction_1 = dat.sum_s_fraction;

dat = load([input2 '/log.mat']);
x1_2 = dat.x1;
x1_2(:,1) = x1_2(:,1) + x1_1(end,1);
sum_i_fraction_2 = dat.sum_i_fraction;
sum_s_fraction_2 = dat.sum_s_fraction;
sample_i = dat.sample_i;
sample_v = dat.sample_v;
sample_vfinal = dat.sample_vfinal;
x1 = [x1_1; x1_2];
sum_i_fraction = [sum_i_fraction_1; sum_i_fraction_2];
sum_s_fraction = [sum_s_fraction_1; sum_s_fraction_2];
save([target_dir '/mean_binding_merged.mat'],'meanBindingTrans','meanBindingFinal','meanBinding', 'x1','sum_i_fraction','sum_s_fraction','sample_i','sample_v','sample_vfinal');
end
%load([output 'virus_traits.mat']);
%VirusesArray = dat_VirusesArray;
%disp (['calculate mean binding avidity from ' output]);
%for d=1:round(VirusesArray(end,2))
%  TF = find(VirusesArray(:,2)<=d-1 & VirusesArray(:,3)>d-1);
%  meanBinding(d,1) = mean(VirusesArray(TF,7));
%  meanBindingFinal(d,1) = mean(VirusesArray(TF,8));
%end
%save([output 'mean_binding.mat'], 'meanBinding', 'meanBindingFinal');
