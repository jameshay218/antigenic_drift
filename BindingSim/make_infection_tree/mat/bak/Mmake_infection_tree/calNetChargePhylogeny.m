function [ netcharge infectK] = calNetChargePhylogeny( treefilename, chargefilename, tip )
% Calculate Netcharge for viral phylogeny
% The function returns mean netcharge and infectK (immune status) for trunk nodes and non-trunk nodes
% Jul 19, 2013
dat = load(['dat/' treefilename]);
pair = dat.dat_viruses(:,4); % default use pair from whole set
dat = load(['dat/' chargefilename]);
charge_txt = dat.dat_VirusesArray(:,[8 5]);

%tip = 281;
n_tips = [];
n_total = length(pair);
n_tips = (n_total+1)./2;

%tip = pair(end);
%lin_start = round(n_total./2); % Lineage should coalesce to other lineages no later than lin_start
lin_start = 100; % Retrieve only one lineage
anc_set = [];
tip = pair(end:-1:end-100);
anci(1)= {getParentList(pair, tip(1))}; % First lineage
anc_set = union(anc_set,anci{1});
for i=2:length(tip)
anci(i)= {getParentList(pair, tip(i))}; %633 ABJ53482_CY016995_2006/04/06_11_17,
lcai = find(ismember(anci{i},anc_set)==1,1);
if anci{i}(lcai)<lin_start
  anc_set = union(anc_set,anci{i});
  length(anc_set)
end
%anci(find(anci==tip)) = [];
end

%remove root node
anc_set(find(anc_set==0)) = [];

cutoff = 4500000;
nodes = [cutoff:n_total];

% remove non lineage strains

% find the non trunk
nt = nodes(~ismember(nodes, anc_set));

% find the trunk nodes
tk = anc_set;
loc_t = find(tk<=cutoff);
tk(loc_t) = [];
netcharge.total_mean = mean(charge_txt(nodes,1));
netcharge.total_var = var(charge_txt(nodes,1));
netcharge.total_sd = (var(charge_txt(nodes,1))).^0.5;
netcharge.trunk_mean = mean(charge_txt(tk,1));
netcharge.trunk_var = var(charge_txt(tk,1));
netcharge.trunk_sd = (var(charge_txt(tk,1))).^0.5;
netcharge.nontrunk_mean = mean(charge_txt(nt,1));
netcharge.nontrunk_var = var(charge_txt(nt,1));
netcharge.nontrunk_sd = (var(charge_txt(nt,1))).^0.5;
%netcharge.internal_mean = mean(charge_txt(in,1));
%netcharge.internal_var = var(charge_txt(in,1));
%netcharge.external_mean = mean(charge_txt(ex,1));
%netcharge.external_var = var(charge_txt(ex,1));

infectK.total_mean = mean(charge_txt(nodes,2));
infectK.total_var = var(charge_txt(nodes,2));
infectK.total_sd = (var(charge_txt(nodes,2))).^0.5;
infectK.trunk_mean = mean(charge_txt(tk,2));
infectK.trunk_var = var(charge_txt(tk,2));
infectK.trunk_sd = (var(charge_txt(tk,2))).^0.5;
infectK.nontrunk_mean = mean(charge_txt(nt,2));
infectK.nontrunk_var = var(charge_txt(nt,2));
infectK.nontrunk_sd = (var(charge_txt(nt,2))).^0.5;
%infectK.internal_mean = mean(charge_txt(in,2));
%infectK.internal_var = var(charge_txt(in,2));
%infectK.external_mean = mean(charge_txt(ex,2));
%infectK.external_var = var(charge_txt(ex,2));

% Draw histogram
 hist(charge_txt(nodes,1),60);
 hold on;
 line([netcharge.trunk_mean netcharge.trunk_mean], [0 2*10^6]);
 line([netcharge.nontrunk_mean netcharge.nontrunk_mean], [0 2*10^6]);
end




