
function [] = merge_meanbding(input_dir1, input_dir2, input_dir3, target_dir)
input1 = input_dir1;
input2 = input_dir2;
input3 = input_dir3;

dat = load([input1 '/mean_binding_sampled(matrix).mat']);
meanBindingSample1 = dat.mean_sampled_binding;
dat = load([input2 '/mean_binding_sampled(matrix).mat']);
meanBindingSample2 = dat.mean_sampled_binding;
dat = load([input3 '/mean_binding_sampled(matrix).mat']);
meanBindingSample3 = dat.mean_sampled_binding;
meanBindingSample = [meanBindingSample1; meanBindingSample2; meanBindingSample3];

%dat = load([input1 '/mean_binding.mat']);
%meanBindingFinal1 = dat.meanBindingFinal;
%dat = load([input2 '/mean_binding.mat']);
%meanBindingFinal2 = dat.meanBindingFinal;
%dat = load([input3 '/mean_binding.mat']);
%meanBindingFinal3 = dat.meanBindingFinal;
%meanBindingFinal = [meanBindingFinal1; meanBindingFinal2; meanBindingFinal3];


%dat = load([input1 '/mean_binding.mat']);
%meanBinding1 = dat.meanBinding;
%dat = load([input2 '/mean_binding.mat']);
%meanBinding2 = dat.meanBinding;
%dat = load([input3 '/mean_binding.mat']);
%meanBinding3 = dat.meanBinding;
%meanBinding = [meanBinding1; meanBinding2; meanBinding3];



dat = load([input1 '/log.mat']);
x1_1 = dat.x1;
sum_i_fraction_1 = dat.sum_i_fraction;
sum_s_fraction_1 = dat.sum_s_fraction;

dat = load([input2 '/log.mat']);
x1_2 = dat.x1;
x1_2(:,1) = x1_2(:,1) + x1_1(end,1);
sum_i_fraction_2 = dat.sum_i_fraction;
sum_s_fraction_2 = dat.sum_s_fraction;

dat = load([input3 '/log.mat']);
x1_3 = dat.x1;
x1_3(:,1) = x1_3(:,1) + x1_2(end,1);
sum_i_fraction_3 = dat.sum_i_fraction;
sum_s_fraction_3 = dat.sum_s_fraction;

sample_i = dat.sample_i;
sample_v = dat.sample_v;
sample_vfinal = dat.sample_vfinal;

x1 = [x1_1; x1_2; x1_3];
sum_i_fraction = [sum_i_fraction_1; sum_i_fraction_2; sum_i_fraction_3];
sum_s_fraction = [sum_s_fraction_1; sum_s_fraction_2; sum_s_fraction_3];
save([target_dir '/mean_binding_merged.mat'],'meanBindingSample', 'x1','sum_i_fraction','sum_s_fraction','sample_i','sample_v','sample_vfinal');
%save([target_dir '/mean_binding_merged.mat'],'meanBindingTrans','meanBindingFinal','meanBinding', 'x1','sum_i_fraction','sum_s_fraction','sample_i','sample_v','sample_vfinal');
end

