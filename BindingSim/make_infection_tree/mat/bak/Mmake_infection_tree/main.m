function [] = main(tip)
addpath(genpath('../plot_trait_change'));

%proj = '20130812_med/';
proj = '20130814_med_fixb/';
%proj = '20130816_med_n4_c09/';

[tree taxa] = annotate_tree(['dat/' proj 'indiv_genealogy_300.branch.tree']);
[t, n, b, nm] = read_newick_seq(['dat/' proj 'indiv_genealogy_300.tree']);
v = zeros(length(b),1)
for i=1:length(taxa)
    if i<=300
        id = str2num(taxa(1,i).name);
    else
        id = i;
    end
    str = taxa(1,i).annotation;
    str1 = regexp(str, 'Bindingscr.set=\"\d.\d+\"', 'match')
    if ~isempty(char(str1))
        v(id,1)=str2num(cell2mat(regexp(char(str1), '\d.\d+', 'match')));
    else
        v(id,1)=0;
    end
end
pairs = string2pairs(['dat/' proj 'indiv_genealogy_300.topo.tree'],300);
disp pairs

dist_root = [];
n_total = length(pairs);
for i=1:length(pairs)
    if i == 300
      disp test
    end
    [lin]= getParentList(pairs, i);
    distance = 0;
    distance = sum(b(lin(1:end-1)));
    %for j=1:length(lin)
    %    if lin(j) == 0
    %        break;
    %    end
    %    distance = distance + b(lin(j)); 
    %end
    dist_root(i,:) = distance;
end

save(['dat/' proj 'indiv_genealogy_300_taxa.mat'], 't', 'n', 'b', 'nm', 'taxa', 'pairs', 'dist_root', 'v');
disp distance
v(587)=0.4;
v(598)=0.4;
v(599)=0.4;
figure; hold on;
for tip=300:-1:299
    [lin]= getParentList(pairs, tip);
    plot(dist_root(lin(1:end-1)),v(lin(1:end-1)),'.r-');
end
plot(dist_root(1:300),v(1:300),'rx');
end