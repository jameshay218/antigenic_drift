function [ netcharge infectK] = calNetCharge( pairfilename, tip )
% The function return the netcharge and infectK(immune status) for a viral
% lineage that ends at a given tip
% Jul 19, 2013

dat = load(['dat/' pairfilename]);
pair = dat.new_pair;
charge_txt = dat.charge_txt;

%tip = 281;
n_tips = [];
n_total = length(pair);
n_tips = (n_total+1)./2;

[anc]= getParentList(pair, tip); %633 ABJ53482_CY016995_2006/04/06_11_17,
anc(find(anc==tip)) = [];

nodes = [1:n_total];
rm_id = [];
for i=1:length(nodes)
  if(pair(nodes(i))==0)
    rm_id(end+1) = nodes(i);
  end
end
nodes(rm_id) = [];
% remove non lineage strains

% find the non trunk
nt = nodes(~ismember(nodes, anc))

% find the external tips
ex = [1:n_tips];
ex = ex(ismember(ex,nodes));

% find the internal tips of non-trunk
in = intersect([n_tips+1:n_total], nt);

% find the trunk nodes
tk = anc;
tk(find(tk==0))=[];
netcharge.trunk_mean = mean(charge_txt(tk,1));
netcharge.trunk_var = var(charge_txt(tk,1));
netcharge.nontrunk_mean = mean(charge_txt(nt,1));
netcharge.nontrunk_var = var(charge_txt(nt,1));
netcharge.internal_mean = mean(charge_txt(in,1));
netcharge.internal_var = var(charge_txt(in,1));
netcharge.external_mean = mean(charge_txt(ex,1));
netcharge.external_var = var(charge_txt(ex,1));

infectK.trunk_mean = mean(charge_txt(tk,2));
infectK.trunk_var = var(charge_txt(tk,2));
infectK.nontrunk_mean = mean(charge_txt(nt,2));
infectK.nontrunk_var = var(charge_txt(nt,2));
infectK.internal_mean = mean(charge_txt(in,2));
infectK.internal_var = var(charge_txt(in,2));
infectK.external_mean = mean(charge_txt(ex,2));
infectK.external_var = var(charge_txt(ex,2));
end




